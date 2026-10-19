#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
oligopool_primer_design.py - 为 oligopool 设计两侧通用引物结合序列

oligo 结构: [引物F结合序列] + payload + [引物R结合序列(RC)]
本脚本针对 pool 全部序列(含反向互补链)筛选安全随机引物:
  1. GC 含量在 [gc_min, gc_max]
  2. 无 >= max_homopolymer 的同碱基连续(如 AAAAA)
  3. 无 >= max_gc_run 的 G/C 连续
  4. 3′ 端 tail_k-mer 不出现在 pool 中     -- 防 3′ 错配起扩(最关键)
  5. 全长任意 inner_k-mer 不出现在 pool 中 -- 防长片段同源结合
  6. 引物间/自身 3′ 端 dimer_k-mer 互补     -- 防 primer-dimer(新增, --dimer-k 0 关闭)

注: 原版同时检查 11-mer 和 12-mer; 12-mer 命中必然伴随其子 11-mer 命中,
    12-mer 检查是纯冗余, 本版合并为单一 inner_k。
    tail_k / inner_k 由程序按 pool 大小自动选择(可行的最小 k), 不接受指定。

多文件模式 (--fasta 可传多个文件, --mode 选择):
    separate  模式1: 每个文件各自单独一张芯片合成, 依次独立设计,
                    每个文件输出一组引物(自动按文件名命名输出)
    shared    模式2: 多个文库共用一张芯片合成, 每文库分配一对专属引物;
                    所有引物对"全部文库的联合 pool + 其他文库侧翼(含RC)"联合筛选,
                    文库间扩增干扰最小, 用各自引物对即可分离单一文库;
                    结果输出为单个文件

用法示例:
    python3 oligopool_primer_design.py --fasta combined_slice.fa
    python3 oligopool_primer_design.py --fasta xx.fa -n 40 --threads 16
    python3 oligopool_primer_design.py --fasta xx.fa --check safe_primers.txt   # 校验已有引物
    python3 oligopool_primer_design.py --fasta A.fa B.fa                       # 模式1: 依次独立设计
    python3 oligopool_primer_design.py --fasta A.fa B.fa --mode shared         # 模式2: 共芯片联合设计
"""

import argparse
import os
import random
import sys
import time
import multiprocessing as mp

# ========================= 常量 =========================
BASE2BIT = {"A": 0, "C": 1, "G": 2, "T": 3}
RC_TAB = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

# numpy 可选加速(有则 k-mer 构建快一个数量级, 无则纯 python 回退)
try:
    import numpy as np
    _LUT = np.full(256, 255, np.uint8)
    for _b, _v in BASE2BIT.items():
        _LUT[ord(_b)] = _v
except ImportError:
    np = None

# fork 模式下 worker 通过写时复制继承这些全局变量, 避免大集合反复 pickle
_G = {"sets": None, "cfg": None, "rng": None}


def rc_seq(s):
    return s.translate(RC_TAB)[::-1]


def encode_kmer(seq):
    x = 0
    for b in seq:
        x = (x << 2) | BASE2BIT[b]
    return x


# ========================= FASTA =========================
def read_fasta(path):
    seqs, name, buf = [], None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    seqs.append("".join(buf))
                name, buf = line[1:], []
            else:
                buf.append(line)
    if name is not None:
        seqs.append("".join(buf))
    return seqs


# ========================= k-mer 集合 =========================
def build_kmer_sets_np(seqs, ks):
    """numpy 快速路径: 所有序列用 N 拼接(防跨序列假 k-mer), 滚动位编码一次成型"""
    big = "N".join(seqs).upper().encode("ascii")
    codes = _LUT[np.frombuffer(big, np.uint8)]
    n = len(codes)
    c64 = codes.astype(np.uint64)
    bad = np.zeros(n + 1, np.int64)          # bad[i] = 前 i 个碱基中非 ACGT 数
    np.cumsum(codes == 255, out=bad[1:])
    sets = {}
    for k in ks:
        m = n - k + 1
        x = np.zeros(m, np.uint64)
        for j in range(k):
            x = x * np.uint64(4) + c64[j:j + m]
        valid = (bad[k:k + m] - bad[:m]) == 0
        sets[k] = set(np.unique(x[valid]).tolist())
    return sets


def build_kmer_sets_py(seqs, ks):
    """纯 python 回退: 滚动编码, 遇非 ACGT 碱基重置窗口"""
    ks = sorted(ks)
    masks = {k: (1 << (2 * k)) - 1 for k in ks}
    sets = {k: set() for k in ks}
    get = BASE2BIT.get
    for seq in seqs:
        xs = dict.fromkeys(ks, 0)
        run = 0
        for b in seq.upper():
            v = get(b)
            if v is None:
                xs = dict.fromkeys(ks, 0)
                run = 0
                continue
            run += 1
            for k in ks:
                xs[k] = ((xs[k] << 2) | v) & masks[k]
                if run >= k:
                    sets[k].add(xs[k])
    return sets


def build_kmer_sets(seqs, ks):
    if np is not None:
        return build_kmer_sets_np(seqs, ks)
    return build_kmer_sets_py(seqs, ks)


# ========================= 引物评估 =========================
def max_runs(p):
    """返回 (最长同碱基连续, 最长 G/C 连续)"""
    maxh = cur = 1
    maxg = gcrun = (1 if p[0] in "GC" else 0)
    for i in range(1, len(p)):
        cur = cur + 1 if p[i] == p[i - 1] else 1
        if cur > maxh:
            maxh = cur
        if p[i] in "GC":
            gcrun = gcrun + 1 if p[i - 1] in "GC" else 1
            if gcrun > maxg:
                maxg = gcrun
        else:
            gcrun = 0
    return maxh, maxg


def fail_reasons(p, cfg, sets, detail=False):
    """返回不满足的条件列表; detail=True 时附带匹配窗口信息(用于 --check)"""
    reasons = []
    L = len(p)
    gc = (p.count("G") + p.count("C")) / L
    if gc < cfg["gc_min"] or gc > cfg["gc_max"]:
        reasons.append("GC %.3f 超出 [%.2f, %.2f]" % (gc, cfg["gc_min"], cfg["gc_max"]))
    maxh, maxg = max_runs(p)
    if maxh >= cfg["max_homopolymer"]:
        reasons.append("同碱基连续 x%d (>= %d)" % (maxh, cfg["max_homopolymer"]))
    if maxg >= cfg["max_gc_run"]:
        reasons.append("G/C 连续 x%d (>= %d)" % (maxg, cfg["max_gc_run"]))
    clamp = cfg["gc_clamp"]
    if clamp and sum(c in "GC" for c in p[-5:]) < clamp:
        reasons.append("3′ 端 5nt 内 G/C 数 < %d (gc-clamp)" % clamp)
    tk = cfg["tail_k"]
    tail = p[-tk:]
    if encode_kmer(tail) in sets[tk]:
        reasons.append("3′ 端 %d-mer %s 存在于 pool [3′ 错配起扩风险]" % (tk, tail))
    ik = cfg["inner_k"]
    hits = [p[i:i + ik] for i in range(L - ik + 1)
            if encode_kmer(p[i:i + ik]) in sets[ik]]
    if hits:
        msg = "%d-mer 命中 pool: %s" % (ik, ",".join(hits[:3]))
        if len(hits) > 3:
            msg += " ...共 %d/%d 窗口" % (len(hits), L - ik + 1)
        reasons.append(msg)
    return reasons


def fail_reasons_both(p, cfg, sets, detail=False):
    """同时检查正向与反向互补两个形式.
    3′ 侧翼位点的反向互补才是实际使用的反向 PCR 引物, 其 3′ 端
    对应位点序列的 5′ 端头部 -- 两个方向都必须安全."""
    r = fail_reasons(p, cfg, sets, detail)
    r2 = fail_reasons(rc_seq(p), cfg, sets, detail)
    if r2:
        r = r + ["[反向互补形式, 即实际反向引物] " + x for x in r2]
    return r


def dimer_ok(cand, accepted, dimer_k):
    """3′ 端互补检查: cand 尾部 rc 出现在自身或任一已接受引物中(及反向)即拒绝"""
    if dimer_k <= 0:
        return True
    tail_rc = rc_seq(cand[-dimer_k:])
    if tail_rc in cand:
        return False
    for a in accepted:
        if tail_rc in a or rc_seq(a[-dimer_k:]) in cand:
            return False
    return True


# ========================= worker =========================
def _worker(task):
    # 任务 = (任务号, 长度); 随机流按 (seed, 任务号, 长度) 派生:
    # 1) fork 出的 worker 不会共享同一条随机流(原版 bug)
    # 2) 与调度无关, 同一 --seed 多线程也可精确复现
    i, L = task
    cfg, sets = _G["cfg"], _G["sets"]
    rng = random.Random(((cfg["seed"] or 0) << 32) ^ (i << 7) ^ L)
    for _attempt in range(cfg["max_attempts"]):
        p = "".join(rng.choices("ACGT", k=L))
        if not fail_reasons_both(p, cfg, sets):
            return p
    return None


def sample_presence(sets, k, rng, n_sample=2000):
    """抽样估计随机 k-mer 落在 pool 中的比例, 用于可行性预检"""
    s = sets[k]
    hit = 0
    for _ in range(n_sample):
        if encode_kmer("".join(rng.choices("ACGT", k=k))) in s:
            hit += 1
    return hit / float(n_sample)


def pilot_hits(sets, tk, ik, primer_len=20, gc_min=0.40, gc_max=0.60,
               maxh=5, maxg=5, n_try=150000, seed=888777):
    """预试验: 抽 n_try 条随机候选引物, 数出能通过全部安全条件的条数.
    通过率衡量当前 (tk, ik) 是否可行, 比理论估计可靠(基因组 k-mer 分布高度偏态)."""
    rng = random.Random(seed)
    choices = rng.choices
    lo = gc_min * primer_len
    hi = gc_max * primer_len

    def one_ok(s):
        g = s.count("G") + s.count("C")
        if g < lo or g > hi:
            return False
        h, gr = max_runs(s)
        if h >= maxh or gr >= maxg:
            return False
        if encode_kmer(s[-tk:]) in sets[tk]:
            return False
        for i in range(len(s) - ik + 1):
            if encode_kmer(s[i:i + ik]) in sets[ik]:
                return False
        return True

    hit = 0
    for _ in range(n_try):
        p = "".join(choices("ACGT", k=primer_len))
        if one_ok(p) and one_ok(rc_seq(p)):
            hit += 1
    return hit


def resolve_k_and_build(pool_seqs, primer_len=20, verbose=True, n_try=150000):
    """自动确定 (tail_k, inner_k) 并构建 k-mer 集合, 返回 (tk, ik, sets).
    从最严格的 (10, 11) 开始, 用预试验实测通过率, 不足则逐步 +1,
    选到可行的最小 k (最安全). k 为内部参数, 不对外暴露."""
    t0 = time.time()
    tk, ik = 10, 11
    k_max = min(21, primer_len)
    cache = {}
    while True:
        for kk in (tk, ik):
            if kk not in cache:
                cache[kk] = build_kmer_sets(pool_seqs, [kk])[kk]
        sets = {kk: cache[kk] for kk in (tk, ik)}
        hits = pilot_hits(sets, tk, ik, primer_len=primer_len, n_try=n_try)
        if verbose:
            print("  自动尝试 tail-k=%d / inner-k=%d: 预试验 %d 条随机候选通过 %d 条" % (
                tk, ik, n_try, hits))
        if hits >= 1:
            if verbose:
                print("  选定 tail-k=%d / inner-k=%d (可行的最小 k, 用时 %.1f s)" % (
                    tk, ik, time.time() - t0))
            return tk, ik, sets
        if ik >= k_max:
            if verbose:
                print("  警告: k 到 %d/%d 仍无可行候选" % (tk, ik))
            return tk, ik, sets
        tk, ik = tk + 1, ik + 1


# ========================= 生成主流程 =========================
def generate_primers(sets, cfg):
    # 按长度分组生成: cfg["primer_lens"] 为 [L] 或 [F长, R长];
    # seen/dimer 池贯穿全部组, 跨组任取两条也互相兼容
    lens = cfg["primer_lens"]
    n = cfg["n_primers"]
    labels = ["F", "R"] if len(lens) > 1 else ["P"]
    groups = []
    seen, dim_pool = set(), []
    t0 = time.time()
    print("")
    print("[3/4] 生成引物 (n=%d/组, 长度 %s, threads=%d) ..." % (
        n, "/".join(str(x) for x in lens), cfg["threads"]))

    def accept_stream(cands, bucket):
        got = 0
        for c in cands:
            if c is None or c in seen:
                continue
            seen.add(c)
            if not dimer_ok(c, dim_pool, cfg["dimer_k"]) or not dimer_ok(rc_seq(c), dim_pool, cfg["dimer_k"]):
                continue
            bucket.append(c)
            dim_pool.extend([c, rc_seq(c)])
            got += 1
            if len(bucket) >= n:
                break
        return got

    ctx = None
    if cfg["threads"] > 1:
        try:
            ctx = mp.get_context("fork")
        except ValueError:
            ctx = None
        if ctx is None:
            print("警告: 当前平台不支持 fork, 退化为单线程")
            cfg["threads"] = 1

    task_base = 0
    for gi, L in enumerate(lens):
        bucket = []
        if ctx is not None:
            _G["sets"], _G["cfg"] = sets, cfg
            with ctx.Pool(cfg["threads"]) as pool:
                for _round in range(100):
                    batch = max(64, (n - len(bucket)) * 3)
                    got = accept_stream(
                        pool.map(_worker, [(task_base + j, L) for j in range(batch)]),
                        bucket)
                    task_base += batch
                    if len(bucket) >= n or got == 0:
                        break
        else:
            _G["sets"], _G["cfg"] = sets, cfg
            while len(bucket) < n:
                c = _worker((task_base, L))
                task_base += 1
                if c is None:
                    break
                accept_stream([c], bucket)
        groups.append((labels[gi], bucket))
        if len(bucket) < n:
            sys.stderr.write("[错误] 长度 %d 只得到 %d/%d 条引物, 约束可能不可满足;" % (L, len(bucket), n))
            sys.stderr.write("  请尝试放宽 --gc-min/--gc-max 或增大 --max-attempts.")
            if any(g[1] for g in groups):
                _write_outputs(groups, cfg)
                sys.stderr.write("已将部分结果写出, 请调整参数后重试.")
            sys.exit(1)
        print("  长度 %d: 完成 %d 条 (累计用时 %.1f s)" % (L, len(bucket), time.time() - t0))
    return groups

def tm_wallace(p):
    return 2 * (p.count("A") + p.count("T")) + 4 * (p.count("G") + p.count("C"))


def _write_outputs(groups, cfg):
    # groups = [(标签, [引物])]; txt 兼容旧格式, tsv 带组前缀与明细
    with open(cfg["out"], "w") as f:
        for label, seqs in groups:
            for p in seqs:
                f.write(p + chr(10))
    tsv = cfg["out"] + ".tsv"
    with open(tsv, "w") as f:
        f.write("# %s%s" % (cfg["param_line"], chr(10)))
        f.write("name%ssequence%slen%sGC%sTm_Wallace%s" % (chr(9), chr(9), chr(9), chr(9), chr(10)))
        for label, seqs in groups:
            for i, p in enumerate(seqs, 1):
                f.write("%s%02d%s%s%s%d%s%.3f%s%d%s" % (
                    label, i, chr(9), p, chr(9), len(p), chr(9),
                    (p.count("G") + p.count("C")) / len(p), chr(9),
                    tm_wallace(p), chr(10)))
    print("")
    print("[4/4] 结果写入:")
    print("  %s  (每行一条, 兼容旧格式)" % cfg["out"])
    print("  %s  (带 GC/Tm 明细)" % tsv)


def print_report(groups):
    print("")
    print("========== 安全引物 ==========")
    print("%-6s %-24s %6s %6s" % ("name", "sequence", "GC%", "Tm"))
    for label, seqs in groups:
        for i, p in enumerate(seqs, 1):
            print("%s%02d    %-24s %5.1f%% %5dC" % (
                label, i, p, 100 * (p.count("G") + p.count("C")) / len(p),
                tm_wallace(p)))
    print("==============================")
    if len(groups) > 1:
        print("用法提示: 5′ 端从 F 组任取, 3′ 端从 R 组任取(R 的反向互补即为反向引物)")
    else:
        print("用法提示: 任取两条 P_i/P_j; 5′ 端加 P_i, 3′ 端加 P_j(其反向互补即为返端引物)")

# ========================= --check 模式 =========================
def check_primers(sets, cfg, path):
    primers = []
    seen = set()
    with open(path) as f:
        for line in f:
            # 按空白/制表符拆字段, 提取纯 ACGT 且 >=10nt 的 token:
            # 兼容纯列表 / safe_primers.txt.tsv / shared 模式多列 tsv
            for tok in line.split():
                s = tok.upper().strip(",;\"")
                if len(s) >= 10 and s not in seen and all(c in "ACGT" for c in s):
                    seen.add(s)
                    primers.append(s)
    if not primers:
        sys.stderr.write("[错误] %s 中未读到有效引物序列" % path)
        sys.exit(1)
    print("")
    print("校验 %d 条引物 (tail_k=%d, inner_k=%d):" % (
        len(primers), cfg["tail_k"], cfg["inner_k"]))
    print("")
    n_fail = 0
    for i, p in enumerate(primers, 1):
        reasons = fail_reasons_both(p, cfg, sets, detail=True)
        gc = 100 * (p.count("G") + p.count("C")) / len(p)
        if reasons:
            n_fail += 1
            print("P%02d %s GC=%.0f%%  [不通过]" % (i, p, gc))
            for r in reasons:
                print("     - %s" % r)
        else:
            print("P%02d %s GC=%.0f%%  [通过]" % (i, p, gc))
    print("")
    print("汇总: %d/%d 通过" % (len(primers) - n_fail, len(primers)))
    if cfg["dimer_k"] > 0:
        dk = cfg["dimer_k"]
        pairs = []
        for i in range(len(primers)):
            for j in range(i + 1, len(primers)):
                if rc_seq(primers[i]) == primers[j]:
                    continue  # 位点与其自身反向引物互补属设计内, 跳过
                if rc_seq(primers[i][-dk:]) in primers[j] or \
                        rc_seq(primers[j][-dk:]) in primers[i]:
                    pairs.append((i, j))
        self_bad = [i for i, p in enumerate(primers) if rc_seq(p[-dk:]) in p]
        if pairs or self_bad:
            print("引物间 3%s 端互补检查 (dimer_k=%d):" % (chr(8242), dk))
            for i in self_bad:
                print("  P%02d 自身 3%s 端可形成二聚体" % (i + 1, chr(8242)))
            for i, j in pairs:
                print("  P%02d 与 P%02d 3%s 端互补, 同管使用有二聚体风险" % (
                    i + 1, j + 1, chr(8242)))
            n_fail += 1
        else:
            print("引物间/自身 3%s 端互补检查: 无风险" % chr(8242))
    sys.exit(1 if n_fail else 0)


# ========================= 多文件流程 =========================
def flank_conflict(primer, flank_pool, tk, ik):
    """候选引物任一 k-mer 窗口出现在其他文库侧翼序列(含反向互补)中.
    共芯片时其他文库的侧翼会随其扩增富集, 扩增产物流失到别的反应即为
    高拷贝污染源, 因此跨文库侧翼与 payload 采用同样的 k-mer 隔离标准."""
    L = len(primer)
    for k in (tk, ik):
        if k > L:
            continue
        for i in range(L - k + 1):
            w = primer[i:i + k]
            for fs in flank_pool:
                if w in fs:
                    return True
    return False


def _gen_shared(sets, cfg, flank_pool, accepted, rng, L):
    """生成长度 L 的共芯片安全引物: 常规安全检查 + 跨文库侧翼检查 + 全体两两二聚体检查"""
    for _attempt in range(cfg["max_attempts"]):
        p = "".join(rng.choices("ACGT", k=L))
        if fail_reasons_both(p, cfg, sets):
            continue
        if flank_conflict(p, flank_pool, cfg["tail_k"], cfg["inner_k"]):
            continue
        if not dimer_ok(p, accepted, cfg["dimer_k"]) or                 not dimer_ok(rc_seq(p), accepted, cfg["dimer_k"]):
            continue
        return p
    return None


def _read_union(fasta_paths):
    """读取全部输入文件, 返回 (联合序列列表, 各文件名列表), 并打印每文件统计"""
    union, stems = [], []
    for fp in fasta_paths:
        seqs = [s.upper() for s in read_fasta(fp) if s]
        if not seqs:
            sys.stderr.write("[错误] %s 中没有序列" % fp)
            sys.exit(1)
        print("  %s: %d 条, %.2f Mb" % (
            os.path.basename(fp), len(seqs), sum(len(s) for s in seqs) / 1e6))
        union += seqs
        stems.append(os.path.splitext(os.path.basename(fp))[0])
    return union, stems


def process_single(fasta_path, args, cfg):
    """单文件流程(单文件输入与模式1共用); cfg["out"] 由调用方设置"""
    t0 = time.time()
    print("[1/4] 读取 %s ..." % fasta_path)
    seqs = [s.upper() for s in read_fasta(fasta_path) if s]
    if not seqs:
        sys.stderr.write("[错误] FASTA 中没有序列")
        sys.exit(1)
    n_bad = sum(1 for s in seqs for c in s if c not in "ACGTN")
    total_bp = sum(len(s) for s in seqs)
    print("  %d 条序列, 共 %.2f Mb; 非法碱基 %d 个(自动按断点处理)" % (
        len(seqs), total_bp / 1e6, n_bad))
    print("  加入反向互补链后总池 %.2f Mb" % (2 * total_bp / 1e6))

    print("[2/4] 自动确定 k 阈值并构建 k-mer 集合 (%s) ..." % (
        "numpy 加速" if np is not None else "纯 python"))
    pool_seqs = seqs + [rc_seq(s) for s in seqs]
    tk, ik, sets = resolve_k_and_build(pool_seqs, primer_len=cfg["primer_len"])
    cfg["tail_k"], cfg["inner_k"] = tk, ik
    for k in sorted({tk, ik}):
        print("  k=%d: %d 个不同 k-mer" % (k, len(sets[k])))
    print("  用时 %.1f s" % (time.time() - t0))

    rng0 = random.Random(args.seed if args.seed is not None else 20260929)
    for k in sorted({tk, ik}):
        print("  随机 %d-mer 落在 pool 的比例: %.2f%%" % (
            k, 100 * sample_presence(sets, k, rng0)))

    groups = generate_primers(sets, cfg)
    print_report(groups)
    _write_outputs(groups, cfg)
    print("")
    print("总用时 %.1f s" % (time.time() - t0))


def process_shared(fasta_paths, args, cfg):
    """模式2: 多文库共一张芯片. 每文库分配一对专属引物, 全部引物对
    联合 pool(全部文库 payload)+其他文库侧翼联合筛选, 文库间互相隔离."""
    t0 = time.time()
    print("[1/3] 读取 %d 个文库(共芯片):" % len(fasta_paths))
    union, stems = _read_union(fasta_paths)
    total_bp = sum(len(s) for s in union)
    print("  联合 pool: %d 条, %.2f Mb (双链 %.2f Mb)" % (
        len(union), total_bp / 1e6, 2 * total_bp / 1e6))

    print("[2/3] 对联合 pool 自动选 k 并构建 k-mer 集合 (%s) ..." % (
        "numpy 加速" if np is not None else "纯 python"))
    pool = union + [rc_seq(s) for s in union]
    tk, ik, sets = resolve_k_and_build(pool, primer_len=cfg["primer_len"])
    cfg["tail_k"], cfg["inner_k"] = tk, ik
    for k in sorted({tk, ik}):
        print("  k=%d: %d 个不同 k-mer" % (k, len(sets[k])))
    print("  用时 %.1f s" % (time.time() - t0))

    print("[3/3] 逐文库分配引物对 (对联合 pool + 已选侧翼联合筛选; 每文库固定一对, -n 不适用) ...")
    rng = random.Random(args.seed)
    accepted, flanks, pairs = [], [], []
    for stem in stems:
        F = _gen_shared(sets, cfg, flanks, accepted, rng, cfg["f_len"])
        if F is None:
            sys.stderr.write("[错误] 为文库 %s 生成 F 引物失败, 约束过紧;" % stem)
            sys.stderr.write("  可减少文库数/放宽 --gc-min/--gc-max/增大 --max-attempts,")
            sys.stderr.write("  或改用 --mode separate 分芯片设计.")
            sys.exit(1)
        accepted.extend([F, rc_seq(F)])
        flanks.extend([F, rc_seq(F)])
        R = _gen_shared(sets, cfg, flanks, accepted, rng, cfg["r_len"])
        if R is None:
            sys.stderr.write("[错误] 为文库 %s 生成 R 引物失败, 约束过紧;" % stem)
            sys.exit(1)
        accepted.extend([R, rc_seq(R)])
        flanks.extend([R, rc_seq(R)])
        pairs.append((stem, F, R))

    print("")
    print("========== 各文库引物对 ==========")
    for stem, F, R in pairs:
        print("%-18s F=%s R=%s (PCR 反向引物 %s)" % (stem, F, R, rc_seq(R)))
    print("==================================")

    with open(cfg["out"], "w") as f:
        f.write("# mode=shared 多文库共芯片; tail-k=%d inner-k=%d dimer-k=%d; seed=%s%s" % (
            tk, ik, cfg["dimer_k"], args.seed, chr(10)))
        f.write("# library%sF_primer(5侧翼,即PCR正向引物)%sR_site(3侧翼)%sR_pcr_primer%s" % (
            chr(9), chr(9), chr(9), chr(10)))
        for stem, F, R in pairs:
            f.write("%s%s%s%s%s%s%s%s" % (
                stem, chr(9), F, chr(9), R, chr(9), rc_seq(R), chr(10)))
    print("")
    print("[输出] %s" % cfg["out"])
    print("")
    print("用法提示: 每个文库用各自的 F/R 跑 build_oligos.py 出合成单;")
    print("  分离单一文库 = 用该文库的引物对 PCR; 全部引物已对联合 pool 与其他")
    print("  文库侧翼做过 k-mer 隔离, 文库间扩增干扰最小.")
    print("")
    print("总用时 %.1f s" % (time.time() - t0))


def _run_check(fasta_paths, args, cfg):
    """--check: 对全部输入文件校验已有引物(多文件 = 联合 pool, 即共芯片口径)"""
    t0 = time.time()
    print("[读取] %d 个输入文件:" % len(fasta_paths))
    union, _ = _read_union(fasta_paths)
    print("[构建] 对联合 pool 自动选 k 并构建 k-mer 集合 ...")
    pool = union + [rc_seq(s) for s in union]
    tk, ik, sets = resolve_k_and_build(pool, primer_len=cfg["primer_len"])
    cfg["tail_k"], cfg["inner_k"] = tk, ik
    print("  用时 %.1f s" % (time.time() - t0))
    if len(fasta_paths) > 1:
        print("  (多个输入文件已合并为联合 pool 校验, 即共芯片口径)")
    check_primers(sets, cfg, args.check)


# ========================= main =========================
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="oligopool 两侧通用引物结合序列设计",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--fasta", required=True, nargs="+", metavar="FASTA",
                    help="oligopool 序列(FASTA), 可传多个文件(配合 --mode)")
    ap.add_argument("--mode", choices=["separate", "shared"], default="separate",
                    help="多文件模式(仅多文件时生效): separate=模式1,每个文件独立一张芯片;"
                         " shared=模式2,多文库共一张芯片并隔离各文库")
    ap.add_argument("--out", default=None,
                    help="输出文件(单文件/shared 模式生效; 多文件 separate 模式自动按文件名输出)")
    ap.add_argument("-n", "--n-primers", type=int, default=20,
                    help="生成引物条数(单文件/separate 模式; shared 模式每文库固定一对)")
    ap.add_argument("--threads", type=int, default=20, help="并行进程数(1=单线程)")
    ap.add_argument("--primer-len", type=int, nargs="+", default=[20],
                    metavar=("F_LEN", "R_LEN"),
                    help="引物长度; 1 个值=F/R 同长, 2 个值=分别指定 F 和 R 长度")
    ap.add_argument("--gc-min", type=float, default=0.40, help="GC 下限")
    ap.add_argument("--gc-max", type=float, default=0.60, help="GC 上限")
    ap.add_argument("--max-homopolymer", type=int, default=5,
                    help="同碱基连续达到该值即拒绝")
    ap.add_argument("--max-gc-run", type=int, default=5,
                    help="G/C 连续达到该值即拒绝")
    ap.add_argument("--dimer-k", type=int, default=8,
                    help="引物间 3′ 端互补检查长度(0=关闭)")
    ap.add_argument("--gc-clamp", type=int, default=0,
                    help="要求 3′ 端 5nt 内至少 N 个 G/C(0=关闭)")
    ap.add_argument("--seed", type=int, default=None, help="随机种子(可复现)")
    ap.add_argument("--max-attempts", type=int, default=2000000,
                    help="每条引物最大尝试次数(防死循环)")
    ap.add_argument("--check", default=None, metavar="PRIMER_FILE",
                    help="校验模式: 检查已有引物(每行一条)对当前 pool 是否仍安全")
    args = ap.parse_args(argv)

    if not (1 <= len(args.primer_len) <= 2):
        ap.error("--primer-len 需要 1 或 2 个值")
    for v in args.primer_len:
        if not (4 <= v <= 60):
            ap.error("--primer-len 每个值应在 4-60 之间")
    f_len = args.primer_len[0]
    r_len = args.primer_len[-1]
    if not (0 < args.gc_min < args.gc_max < 1):
        ap.error("需满足 0 < gc-min < gc-max < 1")

    cfg = vars(args).copy()
    cfg["param_line"] = " ".join(argv if argv is not None else sys.argv[1:])
    cfg["f_len"], cfg["r_len"] = f_len, r_len
    cfg["primer_len"] = max(f_len, r_len)
    cfg["primer_lens"] = [f_len] if f_len == r_len else [f_len, r_len]
    fastas = args.fasta

    if args.check:
        _run_check(fastas, args, cfg)
        return

    if len(fastas) == 1:
        if cfg["out"] is None:
            cfg["out"] = os.path.join(
                os.path.dirname(os.path.abspath(fastas[0])), "safe_primers.txt")
        process_single(fastas[0], args, cfg)
        return

    if args.mode == "separate":
        print("模式1 (separate): %d 个文件各自独立设计(每个文件单独一张芯片), 依次处理" % len(fastas))
        print("")
        for i2, fp in enumerate(fastas, 1):
            out = os.path.join(os.path.dirname(os.path.abspath(fp)),
                               os.path.splitext(os.path.basename(fp))[0] + "_safe_primers.txt")
            print("===== [%d/%d] %s =====" % (i2, len(fastas), fp))
            cfg["out"] = out
            process_single(fp, args, cfg)
            print("")
    else:
        print("模式2 (shared): %d 个文库共一张芯片, 联合设计并隔离各文库" % len(fastas))
        if cfg["out"] is None:
            cfg["out"] = os.path.join(
                os.path.dirname(os.path.abspath(fastas[0])), "safe_primers_shared.tsv")
        process_shared(fastas, args, cfg)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
oligo.py - 为 oligopool 设计两侧通用引物结合序列

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

用法示例:
    python3 oligo.py --fasta combined_slice.fa
    python3 oligo.py --fasta xx.fa -n 40 --threads 16
    python3 oligo.py --fasta xx.fa --check safe_primers.txt   # 校验已有引物
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
def _worker(i):
    # 每个任务按 (seed, 任务号) 派生独立随机流:
    # 1) fork 出的 worker 不会共享同一条随机流(原版 bug)
    # 2) 与调度无关, 同一 --seed 多线程也可精确复现
    cfg, sets = _G["cfg"], _G["sets"]
    rng = random.Random(((cfg["seed"] or 0) << 32) ^ i)
    L = cfg["primer_len"]
    for _attempt in range(cfg["max_attempts"]):
        p = "".join(rng.choices("ACGT", k=L))
        if not fail_reasons(p, cfg, sets):
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
    hit = 0
    for _ in range(n_try):
        p = "".join(choices("ACGT", k=primer_len))
        g = p.count("G") + p.count("C")
        if g < lo or g > hi:
            continue
        h, gr = max_runs(p)
        if h >= maxh or gr >= maxg:
            continue
        if encode_kmer(p[-tk:]) in sets[tk]:
            continue
        ok = True
        for i in range(primer_len - ik + 1):
            if encode_kmer(p[i:i + ik]) in sets[ik]:
                ok = False
                break
        if ok:
            hit += 1
    return hit


def resolve_k_and_build(pool_seqs, primer_len=20, verbose=True, n_try=150000):
    """自动确定 (tail_k, inner_k) 并构建 k-mer 集合, 返回 (tk, ik, sets).
    从最严格的 (10, 11) 开始, 用预试验实测通过率, 不足则逐步 +1,
    选到可行的最小 k (最安全). k 为内部参数, 不对外暴露."""
    t0 = time.time()
    tk, ik = 10, 11
    k_max = min(21, primer_len)
    while True:
        sets = build_kmer_sets(pool_seqs, [tk, ik])
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
    n = cfg["n_primers"]
    accepted, seen = [], set()
    t0 = time.time()
    print("")
    print("[3/4] 生成引物 (n=%d, threads=%d) ..." % (n, cfg["threads"]))

    def accept_stream(cands):
        got = 0
        for c in cands:
            if c is None or c in seen:
                continue
            seen.add(c)
            if not dimer_ok(c, accepted, cfg["dimer_k"]):
                continue
            accepted.append(c)
            got += 1
            if len(accepted) >= n:
                break
        return got

    if cfg["threads"] > 1:
        try:
            ctx = mp.get_context("fork")
        except ValueError:
            ctx = None
        if ctx is None:
            print("警告: 当前平台不支持 fork, 退化为单线程")
            cfg["threads"] = 1
        else:
            _G["sets"], _G["cfg"] = sets, cfg
            task_base = 0
            with ctx.Pool(cfg["threads"]) as pool:
                for _round in range(100):
                    batch = max(64, (n - len(accepted)) * 3)
                    # pool.map 保持任务号顺序返回, 保证 --seed 下贪心接受顺序确定
                    got = accept_stream(
                        pool.map(_worker, range(task_base, task_base + batch)))
                    task_base += batch
                    if len(accepted) >= n or got == 0:
                        break
    if cfg["threads"] <= 1:
        _G["sets"], _G["cfg"] = sets, cfg
        counter = 0
        while len(accepted) < n:
            c = _worker(counter)
            counter += 1
            if c is None:
                break
            accept_stream([c])

    if len(accepted) < n:
        sys.stderr.write("")
        sys.stderr.write("[错误] 只得到 %d/%d 条引物, 约束可能不可满足:" % (len(accepted), n))
        sys.stderr.write("  请尝试放宽 --gc-min/--gc-max 或增大 --max-attempts.")
        if accepted:
            _write_outputs(accepted, cfg)
            sys.stderr.write("已将部分结果写出, 请调整参数后重试.")
        sys.exit(1)

    print("  完成: %d 条引物, 用时 %.1f s" % (len(accepted), time.time() - t0))
    return accepted


def tm_wallace(p):
    return 2 * (p.count("A") + p.count("T")) + 4 * (p.count("G") + p.count("C"))


def _write_outputs(accepted, cfg):
    with open(cfg["out"], "w") as f:
        for p in accepted:
            f.write(p + chr(10))
    tsv = cfg["out"] + ".tsv"
    with open(tsv, "w") as f:
        f.write("# %s%s" % (cfg["param_line"], chr(10)))
        f.write("name	sequence	len	GC	Tm_Wallace%s" % chr(10))
        for i, p in enumerate(accepted, 1):
            f.write("P%02d	%s	%d	%.3f	%d%s" % (
                i, p, len(p), (p.count("G") + p.count("C")) / len(p),
                tm_wallace(p), chr(10)))
    print("")
    print("[4/4] 结果写入:")
    print("  %s  (每行一条, 兼容旧格式)" % cfg["out"])
    print("  %s  (带 GC/Tm 明细)" % tsv)


def print_report(accepted):
    print("")
    print("========== 安全引物 ==========")
    print("%-6s %-22s %6s %6s" % ("name", "sequence", "GC%", "Tm"))
    for i, p in enumerate(accepted, 1):
        print("P%02d    %-22s %5.1f%% %5dC" % (
            i, p, 100 * (p.count("G") + p.count("C")) / len(p), tm_wallace(p)))
    print("==============================")
    print("用法提示: 任取两条 P_i/P_j; 5′ 端加 P_i, 3′ 端加 P_j(其反向互补即为返端引物)")


# ========================= --check 模式 =========================
def check_primers(sets, cfg, path):
    primers = []
    with open(path) as f:
        for line in f:
            s = "".join(line.split()).upper()
            if s and all(c in "ACGT" for c in s):
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
        reasons = fail_reasons(p, cfg, sets, detail=True)
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


# ========================= main =========================
def main():
    ap = argparse.ArgumentParser(
        description="oligopool 两侧通用引物结合序列设计",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--fasta", required=True, help="oligopool 序列(FASTA)")
    ap.add_argument("--out", default=None,
                    help="输出引物文件(默认在 fasta 同目录 safe_primers.txt)")
    ap.add_argument("-n", "--n-primers", type=int, default=20, help="生成引物条数")
    ap.add_argument("--threads", type=int, default=20, help="并行进程数(1=单线程)")
    ap.add_argument("--primer-len", type=int, default=20, help="引物长度")
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
    args = ap.parse_args()

    if not (4 <= args.primer_len <= 60):
        ap.error("--primer-len 应在 4-60 之间")
    if not (0 < args.gc_min < args.gc_max < 1):
        ap.error("需满足 0 < gc-min < gc-max < 1")

    cfg = vars(args).copy()
    cfg["param_line"] = " ".join(sys.argv[1:])
    if cfg["out"] is None:
        cfg["out"] = os.path.join(
            os.path.dirname(os.path.abspath(args.fasta)), "safe_primers.txt")

    t0 = time.time()
    print("[1/4] 读取 %s ..." % args.fasta)
    seqs = [s.upper() for s in read_fasta(args.fasta) if s]
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
    tk, ik, sets = resolve_k_and_build(pool_seqs, primer_len=args.primer_len)
    cfg["tail_k"], cfg["inner_k"] = tk, ik
    for k in sorted({tk, ik}):
        print("  k=%d: %d 个不同 k-mer" % (k, len(sets[k])))
    print("  用时 %.1f s" % (time.time() - t0))

    if args.check:
        check_primers(sets, cfg, args.check)
        return

    # 信息: 随机 k-mer 落在 pool 的比例(自动选 k 已保证可行性)
    rng0 = random.Random(args.seed if args.seed is not None else 20260929)
    for k in sorted({tk, ik}):
        print("  随机 %d-mer 落在 pool 的比例: %.2f%%" % (
            k, 100 * sample_presence(sets, k, rng0)))

    accepted = generate_primers(sets, cfg)
    print_report(accepted)
    _write_outputs(accepted, cfg)
    print("")
    print("总用时 %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main()

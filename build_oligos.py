#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_oligos.py - 接受 payload(oligo) FASTA 与两条引物序列, 生成可直接送合成的完整序列

oligo 结构:
    [primer-F] + payload + [primer-R 位点]

约定 (与 oligo.py 配套):
    --primer-f : 直接加在 5′ 端的序列, 同时也是 PCR 正向引物
    --primer-r : 直接加在 3′ 端的序列; PCR 反向引物是它的反向互补
                 (若手上拿的已经是反向引物本身, 加 --primer-r-is-rc 自动取反)

输出 (以 --out 为前缀):
    out.fa               完整序列(保留原始名称)
    out.csv              name,sequence 两列, 可直接提交合成厂商
    out.pcr_primers.txt  实际使用的 PCR 引物对 (F 与 rc(R))

默认在写出前先对两条引物做 pool 安全性校验(tail-k/inner-k 错配 + 二聚体),
不安全则拒绝输出 (--force 可强制写出, --no-verify 可跳过校验)。

用法:
    python3 build_oligos.py --fasta combined_slice.fa --primer-f AGCTCTGTCGCTACGAACGG --primer-r CATCACCAGCGTCCGTAAGC
"""

import argparse
import csv
import os
import sys
import time

from oligo import rc_seq, encode_kmer, resolve_k_and_build


def is_acgt(s):
    return all(c in "ACGT" for c in s)


def read_fasta_named(path):
    """返回 [(名称, 序列)]"""
    recs, name, buf = [], None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:], []
            else:
                buf.append(line)
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def pool_hits(primer, sets, k):
    """返回 primer 中出现在 pool k-mer 集合里的窗口"""
    return [primer[i:i + k] for i in range(len(primer) - k + 1)
            if encode_kmer(primer[i:i + k]) in sets[k]]


def verify_pair(f, r_site, sets, tail_k, inner_k, dimer_k):
    """校验两条引物对 pool 的安全性, 返回问题列表(空=通过)"""
    problems = []
    for label, p in (("primer-F", f), ("primer-R 位点", r_site)):
        tail = p[-tail_k:]
        if encode_kmer(tail) in sets[tail_k]:
            problems.append("%s: 3′ 端 %d-mer %s 存在于 pool [3′ 错配起扩风险]" % (label, tail_k, tail))
        ih = pool_hits(p, sets, inner_k)
        if ih:
            problems.append("%s: %d-mer 命中 pool 共 %d/%d 窗口 (%s ...)" % (
                label, inner_k, len(ih), len(p) - inner_k + 1, ",".join(ih[:3])))
    if dimer_k > 0:
        f_trc = rc_seq(f[-dimer_k:])
        r_trc = rc_seq(r_site[-dimer_k:])
        if f_trc in f or f_trc in r_site or r_trc in f or r_trc in r_site:
            problems.append("引物 3′ 端互补(自身或两条之间), 有 primer-dimer 风险 (dimer_k=%d)" % dimer_k)
    return problems


def main():
    ap = argparse.ArgumentParser(
        description="oligopool 加引物侧翼, 生成送合成完整序列",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--fasta", required=True, help="payload 序列(FASTA)")
    ap.add_argument("--primer-f", required=True, help="5′ 端引物结合序列(即 PCR 正向引物)")
    ap.add_argument("--primer-r", required=True,
                    help="3′ 端引物结合序列(PCR 反向引物为其反向互补)")
    ap.add_argument("--primer-r-is-rc", action="store_true",
                    help="给定的 --primer-r 已是反向引物本身, 自动取反向互补后再拼接")
    ap.add_argument("--out", default=None,
                    help="输出前缀(默认在 fasta 同目录 oligo_for_synthesis)")
    ap.add_argument("--no-verify", action="store_true", help="跳过引物安全性校验")
    ap.add_argument("--force", action="store_true", help="校验不通过也强制写出")
    ap.add_argument("--dimer-k", type=int, default=8, help="校验用二聚体检查长度(0=关)")
    ap.add_argument("--max-len", type=int, default=300,
                    help="完整 oligo 长度告警阈值(仅提示, 不截断)")
    args = ap.parse_args()

    pf = args.primer_f.upper().strip()
    pr = args.primer_r.upper().strip()
    for label, p in (("--primer-f", pf), ("--primer-r", pr)):
        if not is_acgt(p):
            ap.error("%s 含非 ACGT 字符: %s" % (label, p))
    if args.dimer_k and (len(pf) < args.dimer_k or len(pr) < args.dimer_k):
        ap.error("引物长度需 >= --dimer-k = %d" % args.dimer_k)
    r_site = rc_seq(pr) if args.primer_r_is_rc else pr
    if args.out is None:
        args.out = os.path.join(
            os.path.dirname(os.path.abspath(args.fasta)), "oligo_for_synthesis")

    t0 = time.time()
    print("[1/3] 读取 payload %s ..." % args.fasta)
    recs = read_fasta_named(args.fasta)
    if not recs:
        sys.stderr.write("[错误] FASTA 中没有序列")
        sys.exit(1)
    kept, dropped = [], []
    seen = set()
    dup = 0
    for name, seq in recs:
        seq = seq.upper()
        if not seq or not is_acgt(seq):
            dropped.append(name)
            continue
        if name in seen:
            dup += 1
        seen.add(name)
        kept.append((name, seq))
    if dropped:
        print("  跳过 %d 条(空序列或含非 ACGT 碱基, 不能送合成): %s%s" % (
            len(dropped), ",".join(dropped[:5]), " ..." if len(dropped) > 5 else ""))
    if dup:
        print("  警告: 有 %d 个重名序列, 名称未去重" % dup)
    if not kept:
        sys.stderr.write("[错误] 过滤后没有可用序列")
        sys.exit(1)
    print("  保留 %d 条" % len(kept))

    full = [(name, pf + seq + r_site) for name, seq in kept]
    lens = [len(s) for _, s in full]
    over = sum(1 for L in lens if L > args.max_len)
    print("  完整 oligo 长度: min=%d mean=%.1f max=%d (侧翼 +%d/+%d)" % (
        min(lens), sum(lens) / len(lens), max(lens), len(pf), len(r_site)))
    if over:
        print("  警告: %d 条超过 --max-len %d nt, 请确认合成厂商长度上限" % (over, args.max_len))

    if not args.no_verify:
        print("[2/3] 校验引物对当前 pool 的安全性 (k 自动选择) ...")
        pool = [s for _, s in kept] + [rc_seq(s) for _, s in kept]
        tk, ik, sets = resolve_k_and_build(pool, primer_len=len(pf))
        for k in sorted({tk, ik}):
            print("  k=%d: %d 个不同 k-mer" % (k, len(sets[k])))
        problems = verify_pair(pf, r_site, sets, tk, ik, args.dimer_k)
        if problems:
            sys.stderr.write(chr(10) + "[校验不通过]" + chr(10))
            for p in problems:
                sys.stderr.write("  - %s" % p + chr(10))
            if not args.force:
                sys.stderr.write(chr(10) + "已拒绝输出。确认无误可 --force 强制, 或 --no-verify 跳过校验。" + chr(10))
                sys.exit(1)
            sys.stderr.write(chr(10) + "[--force] 继续写出。" + chr(10))
        else:
            print("  两条引物均通过: 无 3′ 错配, 无内部 k-mer 命中, 无二聚体风险")
    else:
        print("[2/3] 跳过校验 (--no-verify)")

    print("[3/3] 写出 ...")
    with open(args.out + ".fa", "w") as f:
        for name, seq in full:
            f.write(">%s%s" % (name, chr(10)))
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + chr(10))
    with open(args.out + ".csv", "w", newline="") as f:
        w = csv.writer(f, quoting=csv.QUOTE_ALL)
        w.writerow(["name", "sequence"])
        for name, seq in full:
            w.writerow([name, seq])
    with open(args.out + ".pcr_primers.txt", "w") as f:
        f.write("# PCR primer pair for this pool%s" % chr(10))
        f.write("# forward : 5′ 侧翼序列本身%s" % chr(10))
        f.write("# reverse : 3′ 侧翼序列的反向互补%s" % chr(10))
        f.write("F	%s%s" % (pf, chr(10)))
        f.write("R	%s%s" % (rc_seq(r_site), chr(10)))

    print("  %s.fa" % args.out)
    print("  %s.csv" % args.out)
    print("  %s.pcr_primers.txt" % args.out)
    print("")
    print("PCR 引物对:")
    print("  F: %s" % pf)
    print("  R: %s (3′ 侧翼 %s 的反向互补)" % (rc_seq(r_site), r_site))
    print("")
    print("完成: %d 条, %d-%d nt, 总用时 %.1f s" % (
        len(full), min(lens), max(lens), time.time() - t0))


if __name__ == "__main__":
    main()

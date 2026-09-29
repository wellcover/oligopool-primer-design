# oligopool-primer-design

Design **safe universal primer binding sites** for oligonucleotide pool (oligopool) synthesis — for MPRA, CRISPR libraries, barcode pools, and any pooled oligo order where you want to amplify the whole pool with a single primer pair.

## The problem

Pooled oligo orders (Twist, IDT, etc.) are typically synthesized as:

```
[Primer F binding site] + payload (variable region) + [Primer R binding site (RC)]
```

so that the entire pool can be PCR-amplified with one primer pair. The two flanking sequences must not misprime anywhere in the pool — a single 3′-end match to some payload sequence causes off-target amplification and distorts pool representation.

This tool randomly samples candidate primers and keeps only those that are safe against **every sequence in your pool, on both strands**.

## Safety criteria

A candidate primer is accepted only if:

1. GC content within `[gc-min, gc-max]` (default 40–60%)
2. No homopolymer run ≥ `max-homopolymer` (default 5)
3. No G/C-only run ≥ `max-gc-run` (default 5)
4. Its 3′ terminal `tail-k`-mer (default 10) does **not** occur anywhere in the pool (both strands) — the critical anti-mispriming filter
5. No `inner-k`-mer (default 11) of the full primer occurs anywhere in the pool
6. Optional but on by default: no 3′-end complementarity between any two accepted primers or with itself (`dimer-k`, default 8) — prevents primer-dimers when you pick a pair

The pool's reverse complements are always included automatically, so a primer that is safe here is safe in both orientations.

> Note: the classic "check 11-mer and 12-mer separately" scheme is redundant — any 12-mer match implies an 11-mer match. This tool uses a single `--inner-k`, default 11, equivalent to the stricter of the two.

## Requirements

- Python 3.6+ (stdlib only)
- `numpy` **optional** — makes k-mer set construction ~10× faster on large pools (pure-Python fallback is automatic)

## Quick start

```bash
# generate 20 candidate primers with default criteria
python3 oligo.py --fasta my_pool.fa

# pick a pair, e.g. P03 as forward and P11 as reverse
#   5′ flank  = P03
#   3′ flank  = P11          (its reverse complement is the actual reverse primer)

# larger pool? raise the k-mer thresholds
python3 oligo.py --fasta my_pool.fa -n 40 --threads 16 --tail-k 12 --inner-k 13

# reproducible runs (exact same output, single- or multi-threaded)
python3 oligo.py --fasta my_pool.fa --seed 42

# re-validate primers you designed earlier against an UPDATED pool
python3 oligo.py --fasta my_pool_v2.fa --check safe_primers.txt
```

Output: one primer per line (plain `safe_primers.txt`, pipe-friendly) plus a `.tsv` with GC and Wallace Tm for each primer. `--check` exits non-zero if any primer fails, so it can run in CI/snakemake.

## Parameters

| Flag | Default | Meaning |
|---|---|---|
| `--fasta` | (required) | oligopool sequences, FASTA |
| `--out` | fasta dir + `safe_primers.txt` | output file |
| `-n, --n-primers` | 20 | how many primers to generate |
| `--threads` | 20 | worker processes (1 = single-thread) |
| `--primer-len` | 20 | primer length |
| `--gc-min / --gc-max` | 0.40 / 0.60 | GC window |
| `--max-homopolymer` | 5 | reject if any single-base run ≥ this |
| `--max-gc-run` | 5 | reject if any G/C-only run ≥ this |
| `--tail-k` | 10 | 3′ terminal k-mer that must be absent from pool |
| `--inner-k` | 11 | internal k-mer that must be absent from pool |
| `--dimer-k` | 8 | pairwise 3′ complementarity check (0 = off) |
| `--gc-clamp` | 0 | require ≥ N G/C in last 5 nt (0 = off) |
| `--seed` | none | random seed, fully reproducible |
| `--max-attempts` | 2000000 | per-primer attempt cap (no infinite loops) |
| `--check FILE` | — | validate existing primers against the pool |

## Choosing `--tail-k` / `--inner-k` for your pool size

The tool samples random k-mers and reports what fraction fall inside your pool. Rule of thumb:

| Pool size (both strands) | suggested `--tail-k / --inner-k` |
|---|---|
| < 0.1 Mb | 10 / 11 |
| ~1 Mb | 11 / 12 |
| 5–15 Mb | 10 / 11 usually still works (GC-balanced k-mers are the rare ones); if generation fails, use 12 / 13 |
| > 50 Mb | 13 / 14 |

If the constraint set is unsatisfiable, the tool tells you instead of looping forever.

## Performance

42,000 oligos (6.7 Mb pool, 13.4 Mb with reverse complements): **~7 s total** with numpy — k-mer sets built by vectorized rolling 2-bit encoding, shared with worker processes via fork copy-on-write, primers found by parallel rejection sampling.

## 中文说明

给 oligopool（MPRA 文库、CRISPR 文库、barcode 池等）设计两侧通用引物结合位点。oligo 结构为 `[引物F结合序列] + payload + [引物R结合序列(RC)]`，用一对引物即可扩增整个池。

脚本随机采样候选引物，只保留满足以下条件者：GC 40–60%、无同碱基连续≥5、无 G/C 连续≥5、3′ 端 10-mer 与全长 11-mer 均不出现在 pool 任何序列（含反向互补链）中、引物间 3′ 端无互补（防二聚体，`--dimer-k 0` 可关）。

用法：

```bash
python3 oligo.py --fasta my_pool.fa                        # 默认参数生成 20 条
python3 oligo.py --fasta my_pool.fa --tail-k 12 --inner-k 13   # 大 pool 提高阈值
python3 oligo.py --fasta my_pool_v2.fa --check safe_primers.txt  # pool 更新后校验旧引物
```

结果中任取两条 P_i / P_j：5′ 端加 P_i，3′ 端加 P_j（其反向互补即为反向引物）。同一 `--seed` 下结果可精确复现。仅依赖 Python 标准库，装了 numpy 更快。

## License

MIT — see [LICENSE](LICENSE).

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
4. Its 3′ terminal `tail-k`-mer does **not** occur anywhere in the pool (both strands) — the critical anti-mispriming filter
5. No `inner-k`-mer of the full primer occurs anywhere in the pool

Both k values are auto-selected per pool (smallest feasible pair, see below); they are not user-settable.
6. Optional but on by default: no 3′-end complementarity between any two accepted primers or with itself (`dimer-k`, default 8) — prevents primer-dimers when you pick a pair

The pool's reverse complements are always included automatically, so a primer that is safe here is safe in both orientations.

> Note: the classic "check 11-mer and 12-mer separately" scheme is redundant — any 12-mer match implies an 11-mer match. This tool uses a single `--inner-k`, default 11, equivalent to the stricter of the two.

## Requirements

- Python 3.6+ (stdlib only)
- `numpy` **optional** — makes k-mer set construction ~10× faster on large pools (pure-Python fallback is automatic)

## Quick start

```bash
# generate 20 candidate primers with default criteria
# (k thresholds are auto-picked as the smallest feasible pair for your pool)
python3 oligo.py --fasta my_pool.fa

# pick a pair, e.g. P03 as forward and P11 as reverse
#   5′ flank  = P03
#   3′ flank  = P11          (its reverse complement is the actual reverse primer)

# reproducible runs (exact same output, single- or multi-threaded)
python3 oligo.py --fasta my_pool.fa --seed 42

# re-validate primers you designed earlier against an UPDATED pool
python3 oligo.py --fasta my_pool_v2.fa --check safe_primers.txt
```

Output: one primer per line (plain `safe_primers.txt`, pipe-friendly) plus a `.tsv` with GC and Wallace Tm for each primer. `--check` exits non-zero if any primer fails, so it can run in CI/snakemake.

## Build the final synthesis oligos

`build_oligos.py` takes your payload FASTA plus the two chosen primers and emits the full sequences ready to order:

```bash
python3 build_oligos.py --fasta my_pool.fa \
    --primer-f CTCGATGATGAAAACCGTCT --primer-r TCCCGACTAAGCCCATGGAT
```

produces (prefix `--out`, default `oligo_for_synthesis`):

- `out.fa` — `[primer-F] + payload + [primer-R site]`, original names kept, 80-col wrapped
- `out.csv` — `name,sequence`, vendor-order format
- `out.pcr_primers.txt` — the actual PCR primer pair (F as-is, R = reverse complement of the 3′ flank), so nobody has to figure out orientation by hand

Before writing anything it **re-verifies both primers against the exact pool you are about to submit** (3′ tail k-mer, internal k-mer, primer-dimer). If a primer has become unsafe — e.g. your pool grew after the primers were designed — it refuses to write and exits non-zero (`--force` to override, `--no-verify` to skip). It also drops sequences containing non-ACGT characters (not orderable) and warns when oligos exceed `--max-len` (default 300 nt).

Orientation conventions:

- `--primer-f` is appended at the 5′ end and doubles as the forward PCR primer.
- `--primer-r` is appended at the 3′ end as-is; the reverse PCR primer is its reverse complement. If what you have is already the reverse primer itself, pass `--primer-r-is-rc` and it will be reverse-complemented automatically.

## Parameters — `oligo.py`

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
| `--dimer-k` | 8 | pairwise 3′ complementarity check (0 = off) |
| `--gc-clamp` | 0 | require ≥ N G/C in last 5 nt (0 = off) |
| `--seed` | none | random seed, fully reproducible |
| `--max-attempts` | 2000000 | per-primer attempt cap (no infinite loops) |
| `--check FILE` | — | validate existing primers against the pool |

## Parameters — `build_oligos.py`

| Flag | Default | Meaning |
|---|---|---|
| `--fasta` | (required) | payload sequences, FASTA |
| `--primer-f` / `--primer-r` | (required) | flanking primer binding sequences |
| `--primer-r-is-rc` | off | given `--primer-r` is the reverse primer itself; RC it before appending |
| `--out` | fasta dir + `oligo_for_synthesis` | output prefix (`.fa` / `.csv` / `.pcr_primers.txt`) |
| `--no-verify` | off | skip primer-vs-pool safety verification |
| `--force` | off | write outputs even if verification fails |
| `--dimer-k` | 8 | verification: pairwise 3′ complementarity length (0 = off) |
| `--max-len` | 300 | warn when final oligos exceed this length |

## How the k thresholds are chosen

`tail-k` / `inner-k` are **not user-settable** — the tool always auto-selects the smallest feasible pair for your pool: it starts at the strictest 10 / 11, verifies feasibility with a 150,000-candidate pilot run against your actual pool, and relaxes by +1 until candidates pass. The chosen pair is printed on every run. This keeps every design at the strictest standard its pool allows, with no parameters to guess.

Reference values (total bases, both strands; from random and genome-derived test pools):

| Pool size, both strands | typical smallest feasible `--tail-k / --inner-k` |
|---|---|
| under 2,000,000 bp | 10 / 11 |
| 2,000,000 – 16,000,000 bp | 10 / 11 – 11 / 12 |
| 16,000,000 – 60,000,000 bp | 11 / 12 – 13 / 14 |
| above 60,000,000 bp | 13 / 14 or higher |

Genome-derived pools usually allow smaller k than random-sequence pools of the same size: their k-mer distribution is GC-skewed while candidate primers are GC-balanced, so the relevant k-mers are the rare ones. If no k up to the cap works, the tool reports it instead of looping forever.

## Performance

42,000 oligos (6.7 Mb pool, 13.4 Mb with reverse complements): **~7 s total** with numpy — k-mer sets built by vectorized rolling 2-bit encoding, shared with worker processes via fork copy-on-write, primers found by parallel rejection sampling.

## 中文说明

给 oligopool（MPRA 文库、CRISPR 文库、barcode 池等）设计两侧通用引物结合位点。oligo 结构为 `[引物F结合序列] + payload + [引物R结合序列(RC)]`，用一对引物即可扩增整个池。

脚本随机采样候选引物，只保留满足以下条件者：GC 40–60%、无同碱基连续≥5、无 G/C 连续≥5、3′ 端 10-mer 与全长 11-mer 均不出现在 pool 任何序列（含反向互补链）中、引物间 3′ 端无互补（防二聚体，`--dimer-k 0` 可关）。

用法：

```bash
python3 oligo.py --fasta my_pool.fa                        # k 自动选择, 生成 20 条
python3 oligo.py --fasta my_pool_v2.fa --check safe_primers.txt  # pool 更新后校验旧引物

# 选定两条后, 一条命令生成送合成的完整序列 (fa + csv + PCR 引物对)
python3 build_oligos.py --fasta my_pool.fa --primer-f <P_i序列> --primer-r <P_j序列>
```

结果中任取两条 P_i / P_j：5′ 端加 P_i，3′ 端加 P_j（其反向互补即为反向引物）。tail-k / inner-k 由程序自动选择（从最严格的 10/11 出发、预试验实测可行性后取可行的最小 k，每次运行都会打印选了什么），不接受手动指定。`build_oligos.py` 会在写出前自动校验两条引物对当前 pool 的安全性（3′ 错配 / 内部 k-mer / 二聚体），不安全则拒绝输出，避免 pool 更新后误用旧引物。同一 `--seed` 下结果可精确复现。仅依赖 Python 标准库，装了 numpy 更快。

## License

MIT — see [LICENSE](LICENSE).

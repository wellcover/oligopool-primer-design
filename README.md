# oligopool-primer-design

**English** | [简体中文](README.zh-CN.md)

Design **safe universal primer binding sites** for oligonucleotide pool (oligopool) synthesis — for MPRA, CRISPR libraries, barcode pools, and any pooled oligo order where you want to amplify the whole pool with a single primer pair.

Companion script `build_oligos.py` attaches the chosen flanks to your payload and emits synthesis-ready sequences, re-verifying primer safety against the exact pool you are about to order.

## The problem

Pooled oligo orders (Twist, IDT, etc.) are typically synthesized as:

```
[Primer F binding site] + payload (variable region) + [Primer R binding site (RC)]
```

so that the entire pool can be PCR-amplified with one primer pair. The two flanking sequences must not misprime anywhere in the pool — a single 3′-end match to some payload sequence causes off-target amplification and distorts pool representation.

`oligopool-primer-design.py` randomly samples candidate primers and keeps only those that are safe against **every sequence in your pool, on both strands**.

## Safety criteria

A candidate primer is accepted only if:

1. GC content within `[gc-min, gc-max]` (default 40–60%)
2. No homopolymer run ≥ `max-homopolymer` (default 5)
3. No G/C-only run ≥ `max-gc-run` (default 5)
4. Its 3′ terminal tail-k-mer does **not** occur anywhere in the pool (both strands) — the critical anti-mispriming filter
5. No inner-k-mer of the full primer occurs anywhere in the pool
6. No 3′-end complementarity between any two accepted primers, or of a primer with itself (`dimer-k`, default 8; set 0 to disable) — prevents primer-dimers when you pick a pair

The pool's reverse complements are always included automatically, so a primer that is safe here is safe in both orientations.

> Note: the classic "check 11-mer and 12-mer separately" scheme is redundant — any 12-mer match implies an 11-mer match. This tool applies a single internal k-mer rule, auto-selected (see below), which is equivalent to the stricter of the two.

## How the k thresholds are chosen

The tail-k / inner-k values are **fully automatic and not user-settable** — the tool always picks the smallest feasible pair for your pool: it starts at the strictest 10 / 11, verifies feasibility with a 150,000-candidate pilot run against your actual pool, and relaxes by +1 (capped at min(21, primer length)) until candidates pass. The chosen pair is printed on every run. This keeps every design at the strictest standard its pool allows, with no parameters to guess.

Reference values (total bases, both strands; from random and genome-derived test pools):

| Pool size, both strands | typical smallest feasible tail-k / inner-k |
|---|---|
| under 2,000,000 bp | 10 / 11 |
| 2,000,000 – 16,000,000 bp | 10 / 11 – 11 / 12 |
| 16,000,000 – 60,000,000 bp | 11 / 12 – 13 / 14 |
| above 60,000,000 bp | 13 / 14 or higher |

Genome-derived pools usually allow smaller k than random-sequence pools of the same size: their k-mer distribution is GC-skewed while candidate primers are GC-balanced, so the relevant k-mers are the rare ones. If no k up to the cap works, the tool reports it instead of looping forever.

## Requirements

- Python 3.6+ (stdlib only)
- `numpy` **optional** — makes k-mer set construction much faster on large pools (pure-Python fallback is automatic)

## Quick start

```bash
# generate 20 candidate primers with default criteria
# (k thresholds are auto-picked as the smallest feasible pair for your pool)
python3 oligopool-primer-design.py --fasta my_pool.fa

# pick a pair, e.g. P03 as forward and P11 as reverse
#   5′ flank  = P03
#   3′ flank  = P11          (its reverse complement is the actual reverse primer)

# reproducible runs (exact same output, single- or multi-threaded)
python3 oligopool-primer-design.py --fasta my_pool.fa --seed 42

# re-validate primers you designed earlier against an UPDATED pool
python3 oligopool-primer-design.py --fasta my_pool_v2.fa --check safe_primers.txt
```

Output: one primer per line (plain `safe_primers.txt`, pipe-friendly) plus a `.tsv` with GC and Wallace Tm for each primer. `--check` exits non-zero if any primer fails, so it can run in CI/Snakemake.

## Build the final synthesis oligos

`build_oligos.py` takes your payload FASTA plus the two chosen primers and emits the full sequences ready to order:

```bash
python3 build_oligos.py --fasta my_pool.fa \
    --primer-f CTCGATGATGAAAACCGTCT --primer-r TCCCGACTAAGCCCATGGAT
```

produces (prefix `--out`, default `oligo_for_synthesis`):

- `out.fa` — `[primer-F] + payload + [primer-R site]`, original names kept, 80-column wrapped
- `out.csv` — `name,sequence`, vendor-order format
- `out.pcr_primers.txt` — the actual PCR primer pair (F as-is, R = reverse complement of the 3′ flank), so nobody has to figure out orientation by hand

Before writing anything it **re-verifies both primers against the exact pool you are about to submit** (3′ tail k-mer, internal k-mer, primer-dimer; k auto-selected the same way). If a primer has become unsafe — e.g. your pool grew after the primers were designed — it refuses to write and exits non-zero (`--force` to override, `--no-verify` to skip). It also drops sequences containing non-ACGT characters (not orderable) and warns when final oligos exceed `--max-len` (default 300 nt).

Orientation conventions:

- `--primer-f` is appended at the 5′ end and doubles as the forward PCR primer.
- `--primer-r` is appended at the 3′ end as-is; the reverse PCR primer is its reverse complement. If what you have is already the reverse primer itself, pass `--primer-r-is-rc` and it will be reverse-complemented automatically.

## Parameters — `oligopool-primer-design.py`

| Flag | Default | Meaning |
|---|---|---|
| `--fasta` | (required) | oligopool sequences, FASTA |
| `--out` | fasta dir + `safe_primers.txt` | output file |
| `-n, --n-primers` | 20 | how many primers to generate |
| `--threads` | 20 | worker processes (1 = single-threaded) |
| `--primer-len` | 20 | primer length |
| `--gc-min / --gc-max` | 0.40 / 0.60 | GC window |
| `--max-homopolymer` | 5 | reject if any single-base run ≥ this |
| `--max-gc-run` | 5 | reject if any G/C-only run ≥ this |
| `--dimer-k` | 8 | pairwise 3′ complementarity check (0 = off) |
| `--gc-clamp` | 0 | require ≥ N G/C in the last 5 nt (0 = off) |
| `--seed` | none | random seed, fully reproducible |
| `--max-attempts` | 2,000,000 | per-primer attempt cap (no infinite loops) |
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

## Performance

42,000 160-bp oligos (6.7 Mb pool, 13.4 Mb with reverse complements): **~10 s total** with numpy on a typical lab server — k-mer sets built by vectorized rolling 2-bit encoding, shared with worker processes via fork copy-on-write, k auto-selection verified by a 150,000-candidate pilot, primers found by parallel rejection sampling.

## License

MIT — see [LICENSE](LICENSE).

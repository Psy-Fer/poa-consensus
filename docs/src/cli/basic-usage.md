# Basic Usage

## Input and output

```
poa-consensus reads.fa
poa-consensus reads.fq
zcat reads.fa.gz | poa-consensus        # stdin ('-') is the default input
```

Input is uncompressed FASTA or FASTQ, auto-detected. Reads are oriented to the seed
strand automatically. Output is a single FASTA record on stdout:

```
>consensus reads=20 seed=7 band=adaptive
CATCATCATCATCATCATCATCAT
```

- `reads`: number of input reads.
- `seed`: index of the median-length read, which the single-allele path seeds on.
- `band`: `adaptive` (the default), `unbanded`, or the fixed band width when
  `--no-adaptive-band` is given.

With `--multi`, one record is written per allele:

```
>allele_1 reads=12 total_reads=20 seed=3 band=adaptive allele=1/2
...
>allele_2 reads=8 total_reads=20 seed=3 band=adaptive allele=2/2
...
```

`reads` is the allele's own read count; use it to judge whether a minority allele has
enough depth.

With `--fastq`, each record is written as FASTQ instead, with the per-base **support
quality** as the quality string: how decisively the reads back each base over its
strongest alternative (Phred, capped at 60). Below about Q10 the reads disagree there,
for example on a repeat length; see
[Per-base support quality](../library/diagnostics.md#per-base-support-quality).

```
@consensus reads=10 seed=3 band=adaptive
CGGTCGGGCTTCTCTCTCTCCACACACACACACACACACACACACACACACACACACACACAGCGGATAATGCCTGCACTAT
+
BBBBBBBBBB..7B7BBB77..BBBBBBBBBBBBBBBBB77BBBBB7BBBBB'%B77B7BBBBBBBBBBBBBBBBBBBBBBB
```

Here the `'%` (Q6, Q4) mark a CA unit the reads split on.

## Flags

Run `poa-consensus --help` for the full list with descriptions. The main ones:

| Flag | Default | Description |
|---|---|---|
| `-m`, `--multi` | off | Multi-allele mode; one record per detected allele |
| `--fastq` | off | Write FASTQ, with per-base support quality as the quality string |
| `--left-flank SEQ`, `--right-flank SEQ` | off | Return only the region between two flank sequences, excluding partial reads automatically (both required together; use ≥ ~20 bp of unique flank) |
| `-b`, `--band-width N` | 50 | Band half-width floor. The half-width used is `max(N, b + f × read_len)` |
| `--no-adaptive-band` | off | Together with `--band-width 0`, turns banding off entirely (fully unbanded DP) |
| `--adaptive-band-b`, `--adaptive-band-f` | 10, 0.01 | Adaptive band `w = b + f × read_len` |
| `--global` | off | Global alignment (semi-global, with free terminal gaps, is the default) |
| `--match`, `--mismatch`, `--gap-open`, `--gap-extend` | 2, -4, -4, -3 | Affine scoring |
| `--consensus-mode MODE` | `best-fit` | `best-fit` (keep whichever of heaviest-path and majority-frequency the reads support better) or `majority` (force MSA-column majority) |
| `--min-reads N` | 3 | Minimum reads; error below this |
| `--min-coverage-fraction F` | 0 | Node coverage floor; 0 = strict majority `(n/2 + 1).max(2)` |
| `--min-allele-freq F` | 0.2 | Multi-allele: minimum allele frequency (raise to ~0.4 on noisy ONT) |
| `--phasing-bubble-min-span N` | 10 | Multi-allele: minimum arm span for a structural bubble |
| `--homopolymer-min-run N` | 2 | Homopolymer compression (0 = off); see below |
| `--homopolymer-gap OPEN,EXTEND` | -2,-1 | Homopolymer-aware gap costs (`off` to disable); see below |
| `--no-homopolymer` | off | Disable homopolymer handling entirely (plain engine) |
| `-s`, `--seed N` | auto | Seed read index for orientation and multi-allele. The single-allele path re-seeds on the median-length read, so this rarely changes its output |
| `-q`, `--quiet` | off | Suppress warnings and notes; errors are always printed |

The diagnostic thresholds (`--depth-warn-threshold`, `--depth-critical-threshold`,
`--depth-allele-threshold`, `--interior-support-threshold`, `--boundary-margin`,
`--truncation-ratio-threshold`, `--no-long-unbanded-warning`) are listed under
**Diagnostics** in `--help`.

## Common invocations

```
# Single-allele consensus with defaults (adaptive band, semi-global, best-fit,
# homopolymer handling)
poa-consensus reads.fa

# Multi-allele (e.g. diploid STR locus)
poa-consensus reads.fa --multi

# Only the repeat between two flanks (partial reads excluded automatically)
poa-consensus reads.fa --left-flank TCGAGTCCCTCAAGTCCTTC --right-flank CTCCTCAGCTTCCTCAGCCG

# Plain engine without homopolymer handling (about 2.4x less CPU)
poa-consensus reads.fa --no-homopolymer

# Fully unbanded (short reads, or to cross-check a banded result)
poa-consensus reads.fa --band-width 0 --no-adaptive-band

# Suppress all non-error output
poa-consensus reads.fa --quiet
```

`--band-width 0` on its own is **not** unbanded: the adaptive `b + f × read_len` band
stays on. Unbanded
memory grows with read length squared (about 4.8 GB for 20 kb reads), so keep it for
short loci.

## Homopolymer handling

A single base between two length-variable homopolymers (`A12 G A11`), or `(AAAAT)n`-type
repeats with variable A runs, is easy to miscall on noisy reads. Plain affine scoring
prefers two mismatches over a deletion plus an insertion, so the base splits across two
alignment columns. Two options, both on by default, target this:

- `--homopolymer-min-run 2` aligns runs of 2 or more identical bases as a single symbol
  and restores their lengths afterwards.
- `--homopolymer-gap -2,-1` makes indels inside a homopolymer cheaper than elsewhere.

Each produces a candidate that is spliced into the plain consensus only inside
homopolymer runs, and kept only if it fits the reads strictly better, so elsewhere the
output is the plain engine's. It applies to single-allele consensus only (ignored with
`--multi`) and costs about 2.4× the plain engine's CPU. On clean HiFi data it rarely
changes anything; use `--no-homopolymer` to skip it.

## Warnings and notes

By default the CLI prints warnings and notes to stderr. They do not affect the FASTA on
stdout:

```
poa-consensus: warning: consensus: only 7 reads — results may be unreliable; recommend ≥ 10
poa-consensus: note: consensus: 1 structural variant site(s) above frequency threshold (minority arm: 3 read(s)) — if the locus is heterozygous, re-run with --multi
```

The checks are: low depth, coverage gaps, near-zero interior support, competing
structural alleles (a note, suggesting `--multi`), and a suspected truncation (the
consensus much shorter than the median read). The CLI also warns once if reads contain
non-ACGT bases (treated as mismatches), and if automatic seed selection is inconclusive
(it then falls back to the median-length read).

Use `--quiet` to suppress all of these. Errors (empty input, insufficient depth, input
parse failures) are always printed.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Error (empty or unreadable input, insufficient depth, a flank given without the other) |
| 2 | Invalid command-line arguments |

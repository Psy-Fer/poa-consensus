# Basic Usage

## Input and output

```
poa-consensus reads.fa
poa-consensus reads.fq
poa-consensus reads.fa.gz
```

Input is FASTA or FASTQ, auto-detected. Output is a single FASTA record on stdout:

```
>consensus n_reads=20 seed=7 band_width=50
CATCATCATCATCATCATCATCAT
```

The header carries `n_reads`, `seed` (the seed read index), and `band_width` (effective
band width used for the final read).

## Flags

| Flag | Default | Description |
|---|---|---|
| `--multi` | off | Multi-allele mode; outputs one record per detected allele |
| `--band-width N` | 50 | Minimum band width; 0 = unbanded |
| `--no-adaptive-band` | off | Disable adaptive band (adaptive is on by default) |
| `--global` | off | Use global alignment (semi-global is the default) |
| `--min-reads N` | 3 | Minimum reads required; error if below this |
| `--seed N` | auto | Explicit seed index (0-based); overrides automatic selection |
| `--quiet` | off | Suppress warnings and notes; errors always printed |

## Common invocations

```
# Single-allele consensus with defaults (adaptive band, semi-global, seed auto)
poa-consensus reads.fa

# Multi-allele (e.g. diploid STR locus)
poa-consensus reads.fa --multi

# Disable adaptive band (use fixed band-width 50 only)
poa-consensus reads.fa --no-adaptive-band

# Fully unbanded (for short reads or when memory is not a concern)
poa-consensus reads.fa --band-width 0

# Suppress all non-error output
poa-consensus reads.fa --quiet

# Pipe reads from bedpull or samtools
samtools view -b bam_file region | samtools fasta | poa-consensus /dev/stdin
```

## Warnings and notes

By default, the CLI prints warnings and notes to stderr. These do not affect the stdout
FASTA output:

```
NOTE [consensus]: interior support low (min fraction 0.08)
WARN [consensus]: truncation suspected (ratio 0.43; consensus 371bp vs median read 861bp)
```

Use `--quiet` to suppress all of these. Errors (e.g. insufficient depth, empty input) are
always printed regardless of `--quiet`.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Error (empty input, insufficient depth, band too narrow after retries) |

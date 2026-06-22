# STR Analysis Workflow

This guide shows how to use `poa-consensus` as part of an end-to-end STR genotyping
pipeline using [`bedpull`](https://github.com/Psy-Fer/bedpull) to extract reads from a BAM.

## Prerequisites

- `bedpull`: extracts reads overlapping a BED region from a BAM
- `samtools`: for BAM handling
- `poa-consensus`: installed with `--features cli`

## Single-allele locus

```bash
# Extract reads at the HTT CAG repeat (hg38 coordinates)
bedpull -b reads.bam -r chr4:3074877-3074967 > htt_reads.fa

# Build consensus
poa-consensus htt_reads.fa > htt_consensus.fa

# Count CAG units in the consensus
python3 -c "
import re, sys
seq = ''.join(l.strip() for l in sys.stdin if not l.startswith('>'))
print(len(re.findall('CAG', seq)), 'CAG units')
" < htt_consensus.fa
```

## Diploid locus

```bash
bedpull -b reads.bam -r chr4:3074877-3074967 > htt_reads.fa

# Multi-allele mode outputs two FASTA records
poa-consensus htt_reads.fa --multi > htt_alleles.fa

# Count units in each allele
grep -v "^>" htt_alleles.fa | while read seq; do
    python3 -c "import re; print(len(re.findall('CAG', '$seq')), 'CAG')"
done
```

## Haplotagged BAM (HP tag)

For haplotagged BAMs, split reads by HP tag before POA:

```bash
# HP=1 reads
samtools view -b reads.bam chr4:3074877-3074967 \
    | samtools view -b -d HP:1 \
    | samtools fasta \
    | poa-consensus /dev/stdin > htt_hap1.fa

# HP=2 reads
samtools view -b reads.bam chr4:3074877-3074967 \
    | samtools view -b -d HP:2 \
    | samtools fasta \
    | poa-consensus /dev/stdin > htt_hap2.fa
```

## Long repeat expansions (RFC1 example)

For the RFC1 AAAAG repeat, banded DP can silently truncate. Use `pad=20` bp of unique
flanking context and let the automatic truncation retry handle it:

```bash
# Extract with 20 bp of unique flanking (pad=20 in bedpull or equivalent)
bedpull -b reads.bam -r chr4:39348425-39348479 --pad 20 > rfc1_reads.fa

# The CLI retries unbanded automatically if truncation is detected
poa-consensus rfc1_reads.fa
```

The truncation detector fires when the consensus is less than 60% of the median read
length, triggering an unbanded retry. The FASTA header will note the retry:

```
>consensus n_reads=15 seed=3 band_width=0 [truncation_retry]
```

## Checking result quality

```bash
# View any warnings to stderr
poa-consensus reads.fa 2>warnings.txt

# Check for specific warning types
grep "truncation" warnings.txt
grep "competing allele" warnings.txt
grep "low depth" warnings.txt
```

## Recommended flags per data type

| Data type | Recommended flags |
|---|---|
| HiFi single-allele STR | defaults (adaptive band, band-width 50, semi-global) |
| HiFi diploid STR | `--multi` |
| ONT single-allele STR | defaults |
| ONT diploid STR | `--multi` (consider raising `min_allele_freq` in library API) |
| Very long reads (>5 kb) | `--band-width 0` (unbanded; check memory) |
| Reads guaranteed spanning | `--seed 0` or `--seed N` for explicit seed control |

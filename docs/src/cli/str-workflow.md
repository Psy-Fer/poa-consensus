# STR Analysis Workflow

This guide shows how to use `poa-consensus` in an end-to-end STR genotyping pipeline,
using [`bedpull`](https://github.com/Psy-Fer/bedpull) to extract the locus from each read
in a BAM.

## Prerequisites

- `bedpull`: extracts the part of each read aligned to a BED region
- `poa-consensus`: installed with `--features cli`

`bedpull` takes regions as a BED file (`-r`/`--bed`), not a region string. By default it
writes only reads that span the whole region.

## Single-allele locus

```bash
# HTT CAG repeat (hg38). BED is 0-based, half-open.
printf 'chr4\t3074876\t3074967\n' > htt.bed

# Extract the locus plus 30 bp of reference flank on each side
bedpull -b reads.bam -r htt.bed --flanks 30 -o htt_reads.fa

# Consensus of the whole extracted region
poa-consensus htt_reads.fa > htt_consensus.fa

# Count CAG units
python3 -c "
import re, sys
seq = ''.join(l.strip() for l in sys.stdin if not l.startswith('>'))
print(len(re.findall('CAG', seq)), 'CAG units')
" < htt_consensus.fa
```

## Anchoring on flanks

Pass the unique sequence on either side of the repeat to get back just the region between
them. If a meaningful fraction of the reads do not span both flanks (partial reads), they
are excluded automatically, which is much more accurate than letting them into the
graph:

```bash
poa-consensus htt_reads.fa \
    --left-flank TCGAGTCCCTCAAGTCCTTC \
    --right-flank CTCCTCAGCTTCCTCAGCCG > htt_repeat.fa
```

Use at least ~20 bp of flank that does not occur in the repeat, and extract with enough
`--flanks` padding that the reads contain it. Both flags must be given together. If the
flanks cannot be found in the consensus, the full consensus is returned, flanks included.

## Diploid locus

```bash
# Multi-allele mode: one FASTA record per allele
poa-consensus htt_reads.fa --multi > htt_alleles.fa

# Count units in each allele
grep -v "^>" htt_alleles.fa | while read seq; do
    python3 -c "import re; print(len(re.findall('CAG', '$seq')), 'CAG')"
done
```

On noisy ONT data, raise `--min-allele-freq` (to ~0.4) so that sequencing-error bubbles
are not called as a second allele.

## Haplotagged BAM (HP tag)

If the BAM is haplotagged, let `bedpull` split the reads by HP tag and build one
consensus per haplotype:

```bash
bedpull -b reads.bam -r htt.bed --flanks 30 --hap_split -o htt.fa
# writes htt.h0.fa (untagged), htt.h1.fa and htt.h2.fa

poa-consensus htt.h1.fa > htt_hap1.fa
poa-consensus htt.h2.fa > htt_hap2.fa
```

## Homopolymer-rich repeats (RFC1 example)

The RFC1 `(AAAAG)n` repeat, and other repeats with length-variable A runs, are where the
default homopolymer handling matters most on noisy reads. RFC1 is often heterozygous for
an expansion, so combine flank anchoring with `--multi`:

```bash
printf 'chr4\t39348424\t39348479\n' > rfc1.bed
bedpull -b reads.bam -r rfc1.bed --flanks 30 -o rfc1_reads.fa

poa-consensus rfc1_reads.fa --multi \
    --left-flank AGAGCAAGACTCTGTTTCAA \
    --right-flank GCATGTTCTAAAGAGAATA
```

On HG002 ONT reads this returns two alleles of about 65 bp and 580 bp. Homopolymer
handling applies to single-allele calls only, so for the best per-allele call on noisy
reads, split by haplotype first (see above) and run without `--multi` on each file.

The right flank here skips the `AAAAA` immediately after the repeat, since repeat-like
flank sequence anchors poorly; the bases skipped are still part of the output.

For expansions longer than most reads, where few or no reads span the repeat, use
`bridged_consensus` from the library instead (see the Library Guide).

## Checking result quality

Warnings and notes go to stderr:

```bash
poa-consensus reads.fa 2> warnings.txt

grep "only .* read" warnings.txt              # low depth
grep "structural variant site" warnings.txt   # possible second allele: try --multi
grep "of median input read length" warnings.txt  # suspected truncation
```

A suspected truncation means the consensus is much shorter than the median read. Re-run
unbanded to cross-check (`--band-width 0 --no-adaptive-band`), keeping in mind that
unbanded memory grows with read length squared.

## Recommended flags per data type

| Data type | Recommended flags |
|---|---|
| HiFi single-allele STR | defaults (`--no-homopolymer` saves CPU and rarely changes the call) |
| HiFi diploid STR | `--multi` |
| ONT single-allele STR | defaults |
| ONT diploid STR | `--multi --min-allele-freq 0.4`, depth ≥ 10 per allele |
| Reads that may not span the repeat | `--left-flank SEQ --right-flank SEQ` |
| Multi-kb loci (e.g. MUC1 VNTR) | defaults; the adaptive band scales with read length |

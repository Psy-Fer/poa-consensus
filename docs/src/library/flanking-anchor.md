# Flanking-Anchor Preprocessing

## The problem

STR reads extracted from a BAM may enter the POA graph with different rotational phases.
For a GAA repeat, a read starting at the first G, one starting at the first A, and one
starting at the second A all represent the same underlying sequence but look different to the
aligner. This creates spurious branches and unreliable consensus.

The root cause is that nothing pins the read to a common anchor point before alignment.

## The fix: anchor to the flanks

The simplest path is `consensus_flanked` (and `consensus_multi_flanked`), which take the
flank sequences directly and return the consensus of the repeat region **between** them:

```rust
use poa_consensus::{consensus_flanked, PoaConfig};

let left_flank  = b"TTTCTTTCTTTC";   // unique sequence left of the repeat
let right_flank = b"GAAGAAGAAGAA";   // unique sequence right of the repeat

// Consensus of the repeat segment only (flanks excluded from the output).
let result = consensus_flanked(&reads, left_flank, right_flank, &PoaConfig::default())?;
```

These wrappers **auto-detect partial reads**: when a meaningful fraction of reads fail to
span both flanks, the consensus is built only from the reads that do span (trimmed to the
repeat), which is where anchoring wins — a raw consensus is distorted by partial reads. When
the reads already span, it builds from all of them and slices the result, so anchoring is
never worse than the plain call. It falls back to the raw all-reads consensus when too few
reads span to anchor safely. Use flanks of at least ~20 bp of unique sequence.

If you need the extracted segments themselves (e.g. to feed a different pipeline), the
lower-level `extract_flanked_region` returns the zero-copy `&[u8]` slice between the flanks:

```rust
use poa_consensus::{consensus, extract_flanked_region, PoaConfig};

let segments: Vec<&[u8]> = reads
    .iter()
    .filter_map(|read| extract_flanked_region(read, left_flank, right_flank))
    .collect();
let result = consensus(&segments, 0, &PoaConfig::default())?;
```

`extract_flanked_region` aligns the left and right flanks to the read using approximate
matching (allowing up to a few mismatches). The returned slice contains only the bases
between the two flanks. This:

1. Anchors both boundaries to the same reference points in all reads.
2. Excludes non-spanning reads (if either flank is not found, the read is skipped).
3. Eliminates rotational phase ambiguity because all reads start and end at the same
   flanking sequence.

## Design rationale

This approach is used by TRGT, LongTR, and HMMSTR. It is more robust than relying on the
POA aligner to resolve phase by itself, because phase ambiguity is a problem at the
preprocessing stage, not an alignment error the DP can correct.

The flanking sequences must be:
- Unique enough to not appear inside the repeat (a 12-mer from unique flanking context is
  usually sufficient)
- Present in all reads that should be included (reads that only partially span the repeat are
  naturally excluded because one or both flanks won't be found)

## Full pipeline example

```rust
use poa_consensus::{
    extract_flanked_region, auto_orient, select_seed, SeedSelection,
    consensus, diagnose, DiagnoseConfig, PoaConfig,
};

// Step 1: orient all reads to the forward strand
let oriented = auto_orient(&raw_reads, 0);
let oriented_slices: Vec<&[u8]> = oriented.iter().map(|r| r.as_ref()).collect();

// Step 2: extract the repeat segment from each read
let left_flank  = b"TTTCTTTCTTTC";
let right_flank = b"GAAGAAGAAGAA";
let segments: Vec<&[u8]> = oriented_slices.iter()
    .filter_map(|r| extract_flanked_region(r, left_flank, right_flank))
    .collect();

// Step 3: select seed and build consensus
let seed_idx = select_seed(&segments, &SeedSelection::Auto)?;
let result   = consensus(&segments, seed_idx, &PoaConfig::default())?;

// Step 4: check result quality
let warnings = diagnose(&result, &DiagnoseConfig::default());
```

## When not to use this

If you do not have reliable flanking sequences (e.g. the flanks are themselves repetitive,
or you are processing reads without a reference), the core POA still works without this
step. Rotational phase ambiguity on periodic repeats can degrade multi-allele phasing and
close-repeat-count resolution in edge cases, but for most clean single-allele read sets the
consensus is correct.

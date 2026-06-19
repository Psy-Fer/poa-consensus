# Flanking-Anchor Preprocessing

## The problem

STR reads extracted from a BAM may enter the POA graph with different rotational phases.
For a GAA repeat, a read starting at the first G, one starting at the first A, and one
starting at the second A all represent the same underlying sequence but look different to the
aligner. This creates spurious branches and unreliable consensus.

The root cause is that nothing pins the read to a common anchor point before alignment.

## The fix: extract the repeat segment

```rust
use poa_consensus::extract_flanked_region;

let left_flank  = b"TTTCTTTCTTTC";   // unique sequence left of the repeat
let right_flank = b"GAAGAAGAAGAA";   // unique sequence right of the repeat

for read in &reads {
    if let Some(segment) = extract_flanked_region(read, left_flank, right_flank) {
        // segment is a &[u8] slice into the original read;
        // no allocation, zero-copy
        graph.add_read(segment)?;
    }
    // reads where neither flank is found are silently skipped
}
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
step. The rotational phase bugs (CLAUDE.md bugs #1 and #2) will be present in edge cases,
but for most clean read sets the aligner resolves phase by the third or fourth read and the
consensus is correct.

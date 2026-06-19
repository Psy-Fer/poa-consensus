# Bridged Consensus

When a repeat expansion is longer than any single read, no read spans the full locus.
Reads cover the left side (up to the expansion) or the right side (from the expansion
onward) but nothing bridges both ends. A single POA call cannot produce a reliable
consensus because the graph has no shared nodes connecting the two sides.

`bridged_consensus` handles this by building two independent consensuses and joining them
with a gap of unknown size.

## Usage

```rust
use poa_consensus::{bridged_consensus, PoaConfig};

// left_reads: reads that anchor to the left flank
// right_reads: reads that anchor to the right flank
let result = bridged_consensus(
    left_reads,  left_seed_idx,
    right_reads, right_seed_idx,
    &PoaConfig::default(),
)?;

// The gap is in result.gaps
for gap in &result.gaps {
    println!("gap at {}-{}: {:?}", gap.start, gap.end, gap.kind);
}
```

## Output format

The returned `Consensus` contains:
- `sequence`: left consensus bases concatenated with right consensus bases (no bases
  represent the gap itself)
- `gaps`: a single `CoverageGap { start, end, kind: GapKind::Unknown }` where
  `start == end` (zero-length in the sequence, because the gap size is unknown)
- `n_reads`: `left.n_reads + right.n_reads`

Typical rendering:

```
left_seq...(gap:unknown)...right_seq
```

## Detecting the need for bridged consensus

`SeedSelection::Auto` returns `Err(PoaError::NoSpanningReads { left_depth, right_depth })`
when it detects that reads form two non-overlapping groups:

```rust
use poa_consensus::{select_seed, SeedSelection, PoaError};

match select_seed(&reads, &SeedSelection::Auto) {
    Ok(seed_idx) => {
        // normal path
    }
    Err(PoaError::NoSpanningReads { left_depth, right_depth }) => {
        // Split reads into left and right groups and call bridged_consensus
        // left_depth and right_depth are the sizes of the two groups
    }
    Err(e) => return Err(e),
}
```

## Limitations

The gap size is completely unknown. The output indicates that the expansion is too large for
any single read to span, but gives no size estimate. The caller should record the gap as an
unresolved region rather than treating the joined sequence as contiguous.

For long STR expansions (hundreds to thousands of repeat units), anchored assembly is
probably more appropriate than consensus POA. `bridged_consensus` is a best-effort result
for cases where the caller prefers a partial answer over no answer.

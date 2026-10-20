# Seed Selection

The consensus engine seeds itself. `consensus()` and `consensus_multi()` build the graph on
the **median-length read**, whichever `seed_idx` you pass. A length-outlier seed would bias
a periodic repeat toward that outlier's length, and the median read avoids that. The
`seed_idx` argument is only checked for bounds (`Err(PoaError::SeedOutOfBounds)`).

`select_seed` is still useful for two things:

- **Detecting loci that no read spans.** `SeedSelection::Auto` returns
  `Err(PoaError::NoSpanningReads)` when reads form separate left-only and right-only groups,
  which is the signal to use `bridged_consensus`.
- **Choosing an orientation reference.** `auto_orient(&reads, seed_idx)` orients every read
  to the strand of `reads[seed_idx]`, so the reference should be a representative, fully
  spanning read (see [Read Orientation](orientation.md)).

## Automatic selection (`SeedSelection::Auto`)

```rust
use poa_consensus::{select_seed, SeedSelection};

let seed_idx = select_seed(&reads, &SeedSelection::Auto)?;
```

The heuristic uses **terminal k-mer anchors** and performs no alignment:

1. Sample the first and last 50 bp of each read with 15-mers.
2. A k-mer is a valid anchor only if it is end-specific: shared by at least
   `max(2, floor(0.3 × n))` reads at one end and by fewer than that at the other. This
   discards interior repeat k-mers that appear at both ends.
3. Reads with an anchor at both ends are spanning candidates.

**Decision order:**
1. Spanning candidates exist → the shortest spanning read.
2. Otherwise, if there are both left-only and right-only reads →
   `Err(PoaError::NoSpanningReads { left_depth, right_depth })`.
3. Otherwise (no cluster structure, e.g. reads too short or flanks too repetitive) → the
   longest read.

`NoSpanningReads` is a signal to use `bridged_consensus` instead (see
[Bridged Consensus](bridged-consensus.md)).

## Explicit and shortest

```rust
// Validate a caller-chosen index (e.g. a medoid); Err(SeedOutOfBounds) if idx >= reads.len().
let seed_idx = select_seed(&reads, &SeedSelection::Explicit(3))?;

// Shortest read, no k-mer analysis.
let seed_idx = select_seed(&reads, &SeedSelection::Shortest)?;
```

`Shortest` is a fast choice when all reads are known to be spanning (e.g. after
flanking-anchor extraction). None of these choices changes the backbone `consensus()` builds
on.

## Advice for repetitive loci

For STR loci with many short repeat units, interior k-mers are non-specific by definition.
The end-specific filter in `Auto` mode is essential here: without it, interior repeat k-mers
would rank all reads as "spanning" regardless of whether they actually reach both boundaries.

If `Auto` returns `NoSpanningReads` for a locus where reads do overlap, the flanking
sequence is probably itself repetitive and the end-specific filter is too aggressive. In
that case use `Shortest` or `Explicit`, or anchor with known flanks via `consensus_flanked`
(see [Flanking-Anchor Preprocessing](flanking-anchor.md)).

# Seed Selection

The seed read initialises the graph as a linear chain. All other reads are aligned into
this initial structure, so the seed choice affects alignment quality for every subsequent
read.

## Automatic selection (`SeedSelection::Auto`)

```rust
use poa_consensus::{select_seed, SeedSelection};

let seed_idx = select_seed(&reads, &SeedSelection::Auto)?;
```

The automatic heuristic uses **terminal k-mer anchors**:

1. Sample the first and last 50 bp of each read with 15-mers.
2. A k-mer is a valid anchor only if it is end-specific: common at one end of the read set
   and rare at the other. This suppresses interior repeat k-mers that appear at all
   positions.
3. Among reads with anchors at both ends (fully spanning), the shortest is chosen as the
   seed. A shorter seed means fewer initial nodes, which reduces the initial band width
   needed to align longer reads.

**Fallback hierarchy:**
- Reads with anchors at both ends → shortest spanning read
- No spanning reads → longest read
- Reads split into two non-overlapping groups → `Err(PoaError::NoSpanningReads)`

`NoSpanningReads` is a signal to use `bridged_consensus` instead (see
[Bridged Consensus](bridged-consensus.md)).

## Explicit seed

```rust
// Use read index 3 as the seed.
let result = consensus(&reads, 3, &PoaConfig::default())?;
```

Useful when you know which read is the most representative (e.g. the medoid from a
pairwise distance matrix). Pass the index directly and skip `select_seed`.

## Shortest seed

```rust
let seed_idx = select_seed(&reads, &SeedSelection::Shortest)?;
```

Always picks the shortest read regardless of k-mer spanning. Useful as a fast heuristic
when all reads are known to be spanning (e.g. after flanking-anchor extraction).

## Advice for repetitive loci

For STR loci with many short repeat units, interior k-mers are non-specific by definition.
The end-specific filter in `Auto` mode is essential here: without it, interior repeat k-mers
would rank all reads as "spanning" regardless of whether they actually reach both boundaries.

If `Auto` consistently returns `NoSpanningReads` for a locus where reads do overlap, it
likely means the flanking sequence is itself repetitive and the end-specific filter is too
aggressive. In that case, use `Shortest` or pass an explicit index.

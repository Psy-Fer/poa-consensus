# Algorithm Overview

The algorithm runs in three stages for each read set:

1. **Seed initialisation** -- the seed read becomes a linear chain of nodes. Every subsequent
   read is aligned into this growing graph.
2. **Iterative alignment** -- for each read: (a) topological sort, (b) banded DP alignment
   against the graph, (c) update the graph from the alignment traceback.
3. **Consensus extraction** -- follow the heaviest (most-supported) path through the final
   graph, trim low-coverage boundaries, and return the sequence with coverage and statistics.

## Data structures

### Node

Each node carries:

- `base: u8` -- the nucleotide at this position
- `coverage: u32` -- number of reads that matched (not deleted) this node
- `out_edges: Vec<(usize, i32)>` -- `(target_node_index, edge_weight)` pairs

Edge weight records how many reads traversed from this node to the target. Weights are
updated on every alignment: a match or insert op increments the outgoing edge weight.

### PoaGraph

The graph owns all nodes in a flat `Vec<Node>`. Node indices are stable across the life of
the graph; topological order is recomputed before each alignment. The graph also caches a
**spine** -- the current heaviest-path consensus -- which is used by the banded aligner to
centre the DP band.

## Alignment scoring

| Operation | Score |
|---|---|
| Match | +1 |
| Mismatch | -1 |
| Gap open | -2 |
| Gap extend | -1 (per additional base) |

Affine gap penalties are used from the start. The two-component recurrence (M/I/D matrices)
handles gap-open separately from gap-extend, giving meaningfully better accuracy on
homopolymer runs than linear penalties would.

## Critical design decisions

These choices are load-bearing -- changing them without understanding the tests is likely to
break things.

**Delete ops do not increment coverage.** A read that skips a node via a delete traverses
it but does not count as evidence for it. This prevents phase-shifted reads from inflating
boundary-node coverage and blocking the boundary trim from removing nodes that should not be
in the consensus.

**Edge weight normalisation `(weight - 1)`.** The heaviest-path DP accumulates
`edge_weight - 1` rather than `edge_weight`. Without this, a single outlier read with extra
trailing bases always extends the consensus because its edge still scores positively. With
normalisation, only edges supported by two or more reads score positively; single-read
branches are penalised.

**Seed selection is the caller's responsibility.** The library takes an explicit seed index.
The recommended heuristic is the shortest fully-spanning read (see [Seed Selection](../library/seed-selection.md)).
A poor seed degrades alignment quality for every subsequent read because the seed becomes the
initial graph backbone.

**Semi-global alignment for STR reads.** Global alignment has a subtle failure mode at
homopolymer-flanked loci. Use `AlignmentMode::SemiGlobal` (the default) for extracted STR
reads. See [Single-Allele Consensus](../library/single-allele.md) for details.

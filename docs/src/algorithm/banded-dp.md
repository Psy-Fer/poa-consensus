# Banded DP Alignment

## Why banding

The full affine-gap DP over a graph with V nodes and a read of length L requires O(V × L)
cells. For a 150 bp read against a 160-node graph that is ~24,000 cells: trivial. For a
20 kb read against a 20,000-node graph it is ~400 million cells; at 12 bytes per cell
(three i32 scores, M/I/D) that is at least ~4.8 GB. Banding limits each row's DP to a
corridor of columns around a predicted diagonal, reducing memory and work to O(V × band).

## Band width

Band width is controlled by `PoaConfig`:

```rust
pub struct PoaConfig {
    pub band_width: usize,      // minimum half-width
    pub adaptive_band: bool,    // banding switch (see below)
    pub adaptive_band_b: usize, // base term (default 10)
    pub adaptive_band_f: f32,   // fraction of read length (default 0.01)
    // ...
}
```

Banding is **on** when `band_width > 0 || adaptive_band`. To run unbanded you must set
**both** `band_width = 0` and `adaptive_band = false`; `band_width = 0` alone still bands.

When banding is on, the half-width for a read of length `L` is:

```
half = max(band_width, adaptive_band_b + ceil(adaptive_band_f * L), 4)   // adaptive_band
half = max(band_width, 4)                                                // fixed band
```

`adaptive_band` both adds the length term and turns banding on when `band_width` is 0.

**Default:** `band_width = 50`, `adaptive_band = true`. The 50 bp floor prevents silent
repeat-unit loss in long repetitive alleles where the length term alone (16 for a 600 bp
read) would be too narrow. The length term takes over above ~4 kb (e.g. 210 at 20 kb).

## Static-diagonal-union band (anti-fold)

Each row (one graph node, indexed by topological rank) scores only the query columns in a
window `[lo, hi]`. That window is the **union** of two centres, each padded by `half`:

- an **adaptive centre**: one past the best-scoring column of the node's predecessors
  (`max(best_j[pred]) + 1`), which follows wherever the alignment is actually going; and
- a **static graph-geometry diagonal**: `L - remain[t]`, where `remain[t]` is the node's
  distance to the sink along the heaviest out-edge path. This is fixed by graph geometry and
  does **not** move with the adaptive centre.

Because the static diagonal is always included, the sink node's window always contains the
terminal column `L`, so the end cell is **always reachable within the band**. This is the
key property (abPOA's `GET_AD_DP_*` idea): a too-narrow band can never fold the read into a
bogus all-gap alignment and collapse the consensus. At worst it yields a *suboptimal*
alignment, never an empty or truncated-to-nothing one. There is no separate edge detector
and no multi-pass retry.

## Partial reads are unbanded

The static diagonal `L - remain[t]` assumes the read spans the graph. A short partial read
violates that assumption, so banding is disabled for it: a read is banded only when
`L >= 0.8 × graph_span` (the longest node path through the graph). Partial reads are short,
so the unbanded cost is small.

## How the DP is filled

Rows are filled in topological order. For each node:

1. **Accumulate over predecessors.** Each predecessor's banded row is scanned once
   (predecessor outer, column inner), keeping three running maxima per column: the best
   M/I/D at the previous column (feeds a match), the best M/I at the same column (opens a
   delete) and the best D at the same column (extends a delete). These loops run over
   contiguous slices without branches.
2. **M and D pass.** M and D have no dependency along the row, so they are computed in one
   branch-free pass from the accumulators (the match or mismatch score is added once).
3. **I pass.** Inserts depend on the previous column of the same row, so this pass is
   sequential; it also records the row's best column, which centres the next rows' bands.

No backpointers are stored. Traceback re-derives each step's predecessor by replaying the
fill's comparisons in the same order with the same tie-breaking, which costs O(in-degree)
per step instead of extra bytes per cell. Storage is truly banded: each row keeps only its
`[lo, hi]` window in flat M/I/D arrays (12 bytes per cell) with per-row offsets, and the
arrays are kept and reused across reads rather than reallocated.

Under `SemiGlobal` (the default) the alignment may start at any node and end at any node
in the last read column, so graph regions the read does not reach cost nothing; the read
itself is always fully aligned. Under `Global` it must start at a source and end at a sink.

## Memory and speed

DP storage is about `nodes × (2 × half + 8) × 12` bytes, more where the adaptive and static
centres diverge and the window widens. For a 600 bp locus with ~700 graph nodes that is
roughly 0.9 MB at `half = 50` and 1.7 MB at `half = 100`.

Measured peak RSS (whole process, banded defaults). "Plain engine" is homopolymer handling
off; "default" includes it (see
[Consensus Extraction](consensus-extraction.md#homopolymer-handling)):

| Input | Plain engine | Default |
|---|---|---|
| ~5 kb MUC1 VNTR haplotype | ~55 MB | ~65 MB |
| 5 kb random sequence, ONT-like reads | ~174 MB | ~174 MB |
| 20 kb random sequence | ~265 MB | ~272 MB |
| 20 kb, unbanded | at least ~4.8 GB (DP matrices alone) | |

Speed, median per `consensus()` call:

| Input | Plain engine | Default |
|---|---|---|
| HG002 ONT STR loci | ~5 ms | ~10 ms |
| ~5 kb MUC1 haplotype | ~0.2 s | ~0.55 s |
| 20 kb random sequence | ~0.75 s CPU | ~3.3 s CPU |

Homopolymer handling builds extra candidate consensuses, which is where the default's extra
time goes. Alignment dominates runtime. Building with `-C target-cpu=native` (AVX2) adds only about
1.2x, so explicit SIMD is a low-value lever on this structure.

## Known limitation: wrong diagonal on repetitive sequence

In pure-repetitive sequence (e.g. CAG×100) the DP landscape has many near-equally-scoring
diagonals. The band cannot collapse the alignment (the static diagonal keeps the end cell
in-band), but among several equal-cost alignments the DP still selects one, which may differ
from the truth by a repeat unit or two. This is inherent to aligning periodic sequence and is
not specific to banding.

Gross truncation is surfaced by the **truncation detector** in `diagnose()`: if the consensus
is markedly shorter than the median input read length, the truncation warning fires and the
caller can rebuild unbanded (`band_width = 0` and `adaptive_band = false`) and compare. See
[Diagnostics and Warnings](../library/diagnostics.md). For periodic-repeat phase ambiguity
that no band width resolves, flanking-anchor preprocessing
([Flanking-Anchor Preprocessing](../library/flanking-anchor.md)) is the intended fix.

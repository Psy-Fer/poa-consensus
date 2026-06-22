# Banded DP Alignment

## Why banding

The full affine-gap DP over a graph with V nodes and a read of length L requires O(V × L)
cells. For a 150 bp read against a 160-node graph that is ~24,000 cells: trivial. For a
20 kb read against a 21,000-node graph it is ~400 million cells, requiring ~9.6 GB for three
i32 matrices. Banding limits the DP to a corridor of width 2w + 1 around a predicted
diagonal, reducing memory and work to O(w × L).

## Band width

Band width is controlled by `PoaConfig`:

```rust
pub struct PoaConfig {
    pub band_width: usize,     // fixed minimum; 0 = unbanded
    pub adaptive_band: bool,   // use abPOA formula: w = b + f*L
    pub adaptive_band_b: usize, // base (default 10)
    pub adaptive_band_f: f32,  // fraction of read length (default 0.01)
    // ...
}
```

The **adaptive formula** computes `w = max(band_width, b + f * max(read_len, graph_nodes))`.
Using `graph_nodes` rather than just `read_len` corrects for accumulated insert nodes from
prior reads that shift topological ranks away from the read-position diagonal.

**Default:** `adaptive_band = true`, `band_width = 50`. The 50 bp floor prevents silent
repeat-unit loss in long repetitive alleles where the raw adaptive formula (w ≈ 15 for 500 bp
reads) would be too narrow for even a single wrong diagonal.

**Memory table (3 × i32 matrices):**

| Read length | Band width | Memory |
|---|---|---|
| 600 bp | 50 | ~1.4 MB |
| 600 bp | 100 | ~2.9 MB |
| 20 kb | 210 (adaptive) | ~200 MB |
| 20 kb | unbanded | ~9.6 GB |

## Tracking band

The band is not centred on the fixed diagonal (read_pos = graph_rank). Instead it tracks
`best_j` (the read position of the highest-scoring cell in the previous row). This
self-corrects for:

- Phase shifts (reads with a different number of leading repeat units than the seed)
- SV events (a deletion in the read shifts the alignment diagonal by the deletion size)
- Graph expansion from prior reads (insert nodes push topological ranks upward)

The tracking band means that a graph built from reads with CAG×20 and CAG×22 does not
require a band wide enough to hold both simultaneously; the band re-centres at the phase
shift boundary.

## Diagonal skip

For spine nodes (single predecessor, single successor, not at a bubble boundary) the
full three-matrix DP is often unnecessary. If the predecessor's best read position `bj`
matches the node's base `node_base == query[bj]`, the correct cell is a pure match:

```
M[t][bj+1] = M[pred][bj] + 1
```

This O(1) path is taken instead of computing the full row. It fires on roughly 99% of
spine nodes once the graph has converged (after ~3 reads on a clean set). The saving is
significant: a 150 bp spine with 95% skip rate reduces per-read work from ~150 × 2w
cells to ~8 × 2w cells plus 142 O(1) operations.

## Bubble-range j-window

At bubble entry nodes (nodes with two or more outgoing edges above the allele-frequency
threshold) the band is widened by the bubble span. Each arm is processed in its own
narrow window anchored to the bubble entry's j position, then the band contracts back to
spine width at the bubble exit. This concentrates DP work at divergence points without
widening the global band.

## Approaching-edge detector and retry

If `best_j` gets within `GAP_MARGIN` cells of either band edge during a non-trivial row,
the alignment is approaching the band boundary. Rather than completing a potentially wrong
alignment, the aligner returns `Err(BandTooNarrow { required })` early.

`add_read` wraps the alignment in a 3-pass retry loop:

1. **Pass 1** with the configured band width.
2. **Pass 2** with `w = required` (the width the approaching-edge detector estimated).
3. **Pass 3** unbanded (always correct; used only when passes 1 and 2 both approach an edge).

Pass 1 exits early on `BandTooNarrow` so wasted work is minimal. Pass 3 is rare in
practice.

## Stale spine

The spine (heaviest-path consensus cached in the graph) is used to centre the band and to
anchor minimizer-based alignment. Recomputing it after every read is wasteful; the spine
rarely changes once the graph has converged. The library recomputes only when:

- `n_reads - spine_updated_at >= spine_interval`
- After each recompute, if the spine changed by ≤ 3 bases, `spine_interval` doubles
  (capped at 32). Otherwise it resets to 1.

For 100 reads on a stable sequence this produces ~9 recomputes instead of 99. Correctness
is unconditional: the final `consensus()` call always runs a fresh heaviest-path from the
complete graph.

## Known limitation: silent wrong alignment on repetitive sequence

In pure-repetitive sequence (e.g. CAG×100) the DP landscape has many equally-scoring
diagonals. The tracking band can converge to the wrong diagonal without approaching the
band edge, producing a silently truncated consensus. The approaching-edge detector cannot
catch this because the band never gets close to the edge.

The workaround is the **truncation detector** in `diagnose()` and the automatic retry in
`consensus_adaptive`: if the consensus is less than 60% of the median input read length and
banded alignment was used, the alignment is retried unbanded. See
[Adaptive Two-Pass Mode](adaptive.md).

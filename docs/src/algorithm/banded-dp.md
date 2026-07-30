# Banded DP Alignment

## Why banding

The full affine-gap DP over a graph with V nodes and a read of length L requires O(V × L)
cells. For a 150 bp read against a 160-node graph that is ~24,000 cells: trivial. For a
20 kb read against a 21,000-node graph it is ~400 million cells, requiring ~9.6 GB for three
i32 matrices. Banding limits each row's DP to a corridor of columns around a predicted
diagonal, reducing memory and work to O(V × band).

## Band width

Band width is controlled by `PoaConfig`:

```rust
pub struct PoaConfig {
    pub band_width: usize,      // fixed minimum half-width; 0 = unbanded
    pub adaptive_band: bool,    // add the abPOA-style length term
    pub adaptive_band_b: usize, // base (default 10)
    pub adaptive_band_f: f32,   // fraction of read length (default 0.01)
    // ...
}
```

The band half-width for a read of length `L` is:

```
half = max(band_width, adaptive_band_b + ceil(adaptive_band_f * L), 4)
```

**Default:** `adaptive_band = true`, `band_width = 50`. The 50 bp floor prevents silent
repeat-unit loss in long repetitive alleles where the raw adaptive term alone
(≈ 15 for a 500 bp read) would be too narrow.

**Memory table (3 × i32 matrices):**

| Read length | Band half-width | Memory |
|---|---|---|
| 600 bp | 50 | ~1.4 MB |
| 600 bp | 100 | ~2.9 MB |
| 20 kb | ~210 (adaptive) | ~200 MB |
| 20 kb | unbanded | ~9.6 GB |

## Static-diagonal-union band (anti-fold)

Each row (one graph node, indexed by topological rank) scores only the query columns in a
window `[lo, hi]`. That window is the **union** of two centres, each padded by `half`:

- an **adaptive centre** — one past the best-scoring column of the node's predecessors
  (`max(best_j[pred]) + 1`), which follows wherever the alignment is actually going; and
- a **static graph-geometry diagonal** — `L - remain[t]`, where `remain[t]` is the node's
  distance to the sink along the heaviest out-edge path. This is fixed by graph geometry and
  does **not** move with the adaptive centre.

Because the static diagonal is always included, the sink node's window always contains the
terminal column `L`, so the global end cell is **always reachable within the band**. This is
the key property (abPOA's `GET_AD_DP_*` idea): a too-narrow band can never fold the read into
a bogus all-gap alignment and collapse the consensus — at worst it yields a *suboptimal*
alignment, never an empty or truncated-to-nothing one. There is no separate approaching-edge
detector and no multi-pass retry; the static diagonal makes them unnecessary.

Storage is genuinely banded: each row keeps only its `[lo, hi]` window in flat arrays with
per-row offsets, so memory is O(V × band), not O(V × L).

## Partial reads are unbanded

The static diagonal `L - remain[t]` assumes the read spans the graph. A short partial read
(one that does not reach both boundaries) violates that assumption, so banding is disabled
for it: a read is banded only when `L >= 0.8 × graph_span`. Partial reads are short, so the
unbanded cost is small.

## Known limitation: wrong diagonal on repetitive sequence

In pure-repetitive sequence (e.g. CAG×100) the DP landscape has many near-equally-scoring
diagonals. The band cannot collapse the alignment (the static diagonal keeps the end cell
in-band), but among several equal-cost alignments the DP still selects one, which may differ
from the truth by a repeat unit or two. This is inherent to aligning periodic sequence and is
not specific to banding.

Gross truncation is surfaced by the **truncation detector** in `diagnose()`: if the consensus
is markedly shorter than the median input read length, the truncation warning fires and the
caller can rebuild unbanded (`band_width = 0`) and compare. See
[Diagnostics and Warnings](../library/diagnostics.md). For periodic-repeat phase ambiguity
that no band width resolves, flanking-anchor preprocessing
([Flanking-Anchor Preprocessing](../library/flanking-anchor.md)) is the intended fix.

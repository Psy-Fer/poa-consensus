# Visualisation

The `plot` feature provides SVG helpers for inspecting graph structure and alignment
behaviour. All functions return a `String` containing a self-contained SVG.

## Enabling the feature

```toml
[dependencies]
poa-consensus = { version = "0.1", features = ["plot"] }
```

## Graph network

The network plot shows the POA graph topology: spine nodes laid out horizontally, arm
(insert) nodes above or below their attachment points, labels inside nodes.

```rust
use poa_consensus::plot::graph_network_svg;

// Basic graph (spine highlighted in blue)
let svg = graph_network_svg(&graph, None);
std::fs::write("/tmp/graph.svg", &svg)?;

// Overlay a specific read's alignment path (highlighted in orange)
let probe = b"CATCATCAT";
let svg   = graph_network_svg(&graph, Some(probe));
std::fs::write("/tmp/graph_overlay.svg", &svg)?;
```

Edges between spine nodes and arm nodes are drawn as curves so they don't obscure the
spine. The overlay uses a distinct colour to show which nodes a specific read traverses,
useful for debugging alignment of an outlier read.

## Coverage plot

```rust
use poa_consensus::plot::coverage_svg;

let svg = coverage_svg(&result);
std::fs::write("/tmp/coverage.svg", &svg)?;
```

Per-position read depth along the consensus sequence as a line chart. Coverage drops at
gaps or partial-read boundaries are clearly visible.

## Graph statistics

```rust
use poa_consensus::plot::graph_stats_svg;

let svg = graph_stats_svg(&result.graph_stats);
std::fs::write("/tmp/stats.svg", &svg)?;
```

Two-panel bar chart: count statistics on the left (bubble count, node count, edge count),
rate statistics on the right (single-support fraction, edge weight Gini).

## Edge weight histogram

```rust
use poa_consensus::plot::edge_weight_histogram_svg;

let svg = edge_weight_histogram_svg(&graph);
```

Distribution of edge weights across the full graph. A clean single-allele graph has a spike
at the read depth (all spine edges) and a smaller spike at 1 (noise edges). A bimodal
distribution at two different depths suggests two alleles.

## Node coverage histogram

```rust
use poa_consensus::plot::node_coverage_histogram_svg;

let svg = node_coverage_histogram_svg(&graph);
```

Distribution of per-node coverage. Similar interpretation to the edge weight histogram but
counts Match ops rather than traversals.

## Alignment density heatmap

```rust
use poa_consensus::plot::alignment_density_svg;

let svg = alignment_density_svg(&graph);
```

2D histogram of (graph_rank, read_position) alignment density across all reads. The density
should be concentrated along the main diagonal. Off-diagonal density indicates reads that
align to the wrong part of the graph (band too narrow, phase shift, or genuine SV).

## Band with reads overlay

```rust
use poa_consensus::plot::band_with_reads_svg;

let svg = band_with_reads_svg(&graph);
```

Shows the DP band corridor for each read, with per-read alignment paths drawn on top.
Out-of-band segments are shown in red with a connection line to the nearest in-band segment.
Useful for diagnosing `BandTooNarrow` errors or band-width tuning.

## Running the example

```
cargo run --example network_plot --features plot
```

Writes five SVG files to `/tmp/`:
- `poa_network.svg` -- two-allele SNV graph
- `poa_network_with_read.svg` -- same graph with allele-B read overlaid
- `poa_network_noisy.svg` -- single-allele graph with noise branches
- `poa_network_deletion.svg` -- graph showing a deletion in minority reads
- `poa_network_deletion_overlay.svg` -- deletion graph with the deleting read overlaid

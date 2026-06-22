//! Visualisation helpers for POA graphs and consensus output (requires feature `plot`).
//!
//! All functions return an SVG string.  Write it to a file or embed it in HTML.
//!
//! # Example
//!
//! ```rust,no_run
//! use poa_consensus::{PoaGraph, PoaConfig};
//! use poa_consensus::plot::{coverage_svg, graph_network_svg, graph_stats_svg, edge_weight_histogram_svg};
//!
//! let reads: &[&[u8]] = &[b"CATCATCAT", b"CATCATCAT", b"CATCGTCAT"];
//! let mut graph = PoaGraph::new(reads[0], PoaConfig::default()).unwrap();
//! for r in &reads[1..] { graph.add_read(r).unwrap(); }
//!
//! let consensus = graph.consensus().unwrap();
//! std::fs::write("coverage.svg", coverage_svg(&consensus)).unwrap();
//! std::fs::write("stats.svg", graph_stats_svg(&graph.stats())).unwrap();
//! std::fs::write("weights.svg", edge_weight_histogram_svg(&graph.edge_weights())).unwrap();
//! ```

use kuva::backend::svg::SvgBackend;
use kuva::plot::histogram2d::ColorMap as Histogram2DColorMap;
use kuva::prelude::*;

use crate::graph::AlignOp;
use crate::{Consensus, GraphStats, PoaGraph};

// ── Coverage plot ─────────────────────────────────────────────────────────────

/// SVG line plot of per-position coverage along the consensus sequence.
///
/// The x-axis is the consensus position (0-based) and the y-axis is the number
/// of reads with a Match operation at that node.  Low-coverage positions at the
/// boundaries (after trim) are typically absent; a deep trough mid-sequence
/// indicates a bubble or alignment artefact worth investigating.
pub fn coverage_svg(consensus: &Consensus) -> String {
    let data: Vec<(f64, f64)> = consensus
        .coverage
        .iter()
        .enumerate()
        .map(|(i, &c)| (i as f64, c as f64))
        .collect();

    let mean = if data.is_empty() {
        0.0
    } else {
        data.iter().map(|(_, y)| y).sum::<f64>() / data.len() as f64
    };

    let line = LinePlot::new()
        .with_data(data)
        .with_color("steelblue")
        .with_stroke_width(1.5)
        .with_fill()
        .with_fill_opacity(0.15)
        .with_legend("coverage");

    let mean_line = LinePlot::new()
        .with_data(vec![
            (0.0, mean),
            (consensus.coverage.len().saturating_sub(1) as f64, mean),
        ])
        .with_color("#e05c5c")
        .with_stroke_width(1.0)
        .with_dashed()
        .with_legend("mean");

    let plots: Vec<Plot> = vec![line.into(), mean_line.into()];
    let layout = Layout::auto_from_plots(&plots)
        .with_title("Consensus coverage")
        .with_x_label("Position (bp)")
        .with_y_label("Coverage (reads)");

    render_to_svg(plots, layout)
}

// ── Graph statistics bar chart ────────────────────────────────────────────────

/// SVG bar chart summarising the key `GraphStats` fields.
///
/// Two charts are returned in a single SVG: one for raw counts (nodes, edges,
/// bubbles) and one for normalised fractions (×100 so they sit on a [0,100]
/// scale).  Both are rendered as separate `BarPlot` structs.
///
/// Fractions multiplied by 100 for readability:
/// - `single_support_%`: fraction of nodes with only one read
/// - `gini×100`: edge-weight Gini coefficient
pub fn graph_stats_svg(stats: &GraphStats) -> String {
    let counts = BarPlot::new()
        .with_bars(vec![
            ("nodes", stats.node_count as f64),
            ("edges", stats.edge_count as f64),
            ("bubbles", stats.bubble_count as f64),
            ("max_bubble_depth", stats.max_bubble_depth as f64),
        ])
        .with_color("steelblue");

    let fractions = BarPlot::new()
        .with_bars(vec![
            ("cov_mean", stats.coverage_mean),
            ("single_support_%", stats.single_support_fraction * 100.0),
            ("gini×100", stats.edge_weight_gini * 100.0),
            ("entropy_mean", stats.mean_column_entropy * 100.0),
        ])
        .with_color("seagreen");

    let counts_plots: Vec<Plot> = vec![counts.into()];
    let fractions_plots: Vec<Plot> = vec![fractions.into()];

    let layout_counts = Layout::auto_from_plots(&counts_plots)
        .with_title("Counts")
        .with_y_label("Count");
    let layout_fractions = Layout::auto_from_plots(&fractions_plots)
        .with_title("Rates (×100)")
        .with_y_label("Value");

    let scene = Figure::new(1, 2)
        .with_plots(vec![counts_plots, fractions_plots])
        .with_layouts(vec![layout_counts, layout_fractions])
        .with_title("Graph statistics")
        .render();

    SvgBackend.render_scene(&scene)
}

// ── Edge weight histogram ─────────────────────────────────────────────────────

/// SVG histogram of the edge weight distribution in the graph.
///
/// Each bar represents the count of edges with a weight in that bin.  A
/// heavily right-skewed distribution (most edges at weight 1, a few at high
/// weight) is typical of noisy or heterozygous read sets.  A more uniform
/// distribution suggests a clean single-allele set.
///
/// `weights` is `PoaGraph::edge_weights()`.
pub fn edge_weight_histogram_svg(weights: &[i32]) -> String {
    if weights.is_empty() {
        return empty_svg("Edge weights (empty graph)");
    }

    let data: Vec<f64> = weights.iter().map(|&w| w as f64).collect();
    let min = data.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    let n_bins = ((max - min + 1.0) as usize).clamp(1, 50);

    let hist = Histogram::new()
        .with_data(data)
        .with_bins(n_bins)
        .with_range((min - 0.5, max + 0.5))
        .with_color("steelblue")
        .with_legend("edge weight");

    let plots: Vec<Plot> = vec![hist.into()];
    let layout = Layout::auto_from_plots(&plots)
        .with_title("Edge weight distribution")
        .with_x_label("Weight (reads)")
        .with_y_label("Count");

    render_to_svg(plots, layout)
}

// ── Node coverage histogram ───────────────────────────────────────────────────

/// SVG histogram of per-node coverage across the full graph (in topological
/// order).
///
/// Complements [`coverage_svg`] which plots only the consensus-path nodes.
/// This histogram shows the full graph, including nodes on minority bubble arms
/// that were not selected by the heaviest-path consensus.
///
/// `coverages` is `PoaGraph::node_coverages()`.
pub fn node_coverage_histogram_svg(coverages: &[u32]) -> String {
    if coverages.is_empty() {
        return empty_svg("Node coverage (empty graph)");
    }

    let data: Vec<f64> = coverages.iter().map(|&c| c as f64).collect();
    let max = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let n_bins = (max as usize + 1).clamp(1, 50);

    let hist = Histogram::new()
        .with_data(data)
        .with_bins(n_bins)
        .with_range((0.5, max + 0.5))
        .with_color("mediumpurple")
        .with_legend("node coverage");

    let plots: Vec<Plot> = vec![hist.into()];
    let layout = Layout::auto_from_plots(&plots)
        .with_title("Node coverage distribution")
        .with_x_label("Coverage (reads)")
        .with_y_label("Node count");

    render_to_svg(plots, layout)
}

// ── Alignment density (2D histogram) ─────────────────────────────────────────

/// SVG 2D histogram of where all reads align in (graph_rank, read_position)
/// space.
///
/// Each alignment operation is projected to a point:
/// - **Match** at graph node with topological rank `r` while consuming read
///   position `p` → point `(r, p)`.
/// - **Delete** (read skips a node) → point `(r, p)` with `p` unchanged.
/// - **Insert** (read base with no graph node) → omitted; no graph rank to
///   plot against.
///
/// Density reveals structure:
/// - **Clean single-allele set** — points cluster on a tight diagonal.
/// - **Bubble / heterozygous site** — the diagonal splits into two parallel
///   corridors at the bubble position.
/// - **Insertions in some reads** — vertical smear (read position advances
///   without graph rank advancing).
/// - **Deletions in some reads** — horizontal smear (graph rank advances
///   without read position advancing).
///
/// Reads that fail alignment (e.g. `BandTooNarrow`) are silently skipped.
pub fn alignment_density_svg(graph: &PoaGraph, reads: &[&[u8]]) -> String {
    let n_nodes = graph.node_count();
    if n_nodes == 0 || reads.is_empty() {
        return empty_svg("Alignment density (empty)");
    }

    let max_read_len = reads.iter().map(|r| r.len()).max().unwrap_or(0);

    let mut points: Vec<(f64, f64)> = Vec::new();

    for read in reads {
        let Ok((ops, _, rank_of)) = graph.align_read_ops(read) else {
            continue;
        };
        let mut read_pos: usize = 0;

        for op in &ops {
            match op {
                AlignOp::Match(node_idx) => {
                    points.push((rank_of[*node_idx] as f64, read_pos as f64));
                    read_pos += 1;
                }
                AlignOp::Delete(node_idx) => {
                    points.push((rank_of[*node_idx] as f64, read_pos as f64));
                }
                AlignOp::Insert(_) => {
                    read_pos += 1;
                }
            }
        }
    }

    if points.is_empty() {
        return empty_svg("Alignment density (no alignments)");
    }

    // Bin size: aim for ~1 cell per node on x, scaled proportionally on y.
    let bins_x = n_nodes.clamp(10, 200);
    let bins_y = max_read_len.clamp(10, 200);

    let hist = Histogram2D::new()
        .with_data(
            points,
            (0.0, n_nodes as f64),
            (0.0, max_read_len as f64),
            bins_x,
            bins_y,
        )
        .with_color_map(Histogram2DColorMap::Viridis);

    let plots: Vec<Plot> = vec![hist.into()];
    let layout = Layout::auto_from_plots(&plots)
        .with_title("Alignment density")
        .with_x_label("Graph node (topological rank)")
        .with_y_label("Read position");

    render_to_svg(plots, layout)
}

// ── Band corridor ─────────────────────────────────────────────────────────────

/// SVG showing the alignment band corridor and traceback path for a single
/// read.
///
/// The band (grey shaded region) spans `[rank - w, rank + w]` in read
/// coordinates at each graph rank `rank`, where `w` is the effective band
/// width for this read.  The traceback path (coloured line) shows where the
/// alignment actually went within that corridor.
///
/// Useful for diagnosing:
/// - Whether the band is unnecessarily wide (lots of empty space between path
///   and band edges).
/// - Whether the band is too narrow (path hugging one edge — a sign that
///   `BandTooNarrow` is imminent on harder reads).
///
/// Returns an empty SVG with an explanatory title when the graph is empty or
/// alignment fails.
pub fn band_svg(graph: &PoaGraph, read: &[u8]) -> String {
    let n_nodes = graph.node_count();
    if n_nodes == 0 {
        return empty_svg("Band (empty graph)");
    }

    let (ops, w, rank_of) = match graph.align_read_ops(read) {
        Ok(v) => v,
        Err(e) => return empty_svg(&format!("Band (alignment error: {e})")),
    };

    // Build the traceback path as (graph_rank, read_pos) pairs.
    let mut path: Vec<(f64, f64)> = Vec::new();
    let mut read_pos: usize = 0;

    for op in &ops {
        match op {
            AlignOp::Match(node_idx) => {
                path.push((rank_of[*node_idx] as f64, read_pos as f64));
                read_pos += 1;
            }
            AlignOp::Delete(node_idx) => {
                path.push((rank_of[*node_idx] as f64, read_pos as f64));
            }
            AlignOp::Insert(_) => {
                read_pos += 1;
            }
        }
    }

    // Band corridor: at each graph rank x, valid read rows are [x-w, x+w].
    let unbanded = w == usize::MAX;
    let x_vals: Vec<f64> = (0..n_nodes).map(|i| i as f64).collect();
    let lower: Vec<f64> = x_vals
        .iter()
        .map(|&x| {
            if unbanded {
                0.0
            } else {
                (x - w as f64).max(0.0)
            }
        })
        .collect();
    let upper: Vec<f64> = x_vals
        .iter()
        .map(|&x| {
            if unbanded {
                read.len() as f64
            } else {
                (x + w as f64).min(read.len() as f64)
            }
        })
        .collect();

    let band = BandPlot::new(x_vals, lower, upper)
        .with_color("#aaaaaa")
        .with_opacity(0.35)
        .with_legend("band corridor");

    let path_line = LinePlot::new()
        .with_data(path)
        .with_color("steelblue")
        .with_stroke_width(1.5)
        .with_legend("alignment path");

    let title = if unbanded {
        "Band corridor (unbanded)".to_string()
    } else {
        format!("Band corridor (w = {w})")
    };

    let plots: Vec<Plot> = vec![band.into(), path_line.into()];
    let layout = Layout::auto_from_plots(&plots)
        .with_title(title)
        .with_x_label("Graph node (topological rank)")
        .with_y_label("Read position");

    render_to_svg(plots, layout)
}

// ── Band corridor with all reads overlaid ─────────────────────────────────────

/// SVG showing the band corridor with every read's alignment path overlaid.
///
/// The band corridor is computed from `reads[seed_idx]` (the seed read that
/// was used to initialise the graph).  Every read's traceback path is drawn
/// as a separate coloured line using the Tol muted palette.
///
/// This is the primary diagnostic for band-related problems:
/// - A path that **hugs the top or bottom edge** of the corridor is at risk of
///   `BandTooNarrow` on harder reads.
/// - A path that **touches the edge** and then snaps back is a sign that the
///   band clipped the optimal route — the consensus may be subtly wrong for
///   that read.
/// - **Outlier reads** whose paths diverge sharply from the majority diagonal
///   are visible immediately as isolated colour streaks far from the centre.
///
/// Reads that fail alignment (e.g. already `BandTooNarrow`) are drawn in a
/// dashed red line so they are visible rather than silently absent.
pub fn band_with_reads_svg(graph: &PoaGraph, reads: &[&[u8]], seed_idx: usize) -> String {
    let n_nodes = graph.node_count();
    if n_nodes == 0 || reads.is_empty() {
        return empty_svg("Band + reads (empty graph)");
    }

    let seed = reads[seed_idx.min(reads.len() - 1)];
    let w = match graph.align_read_ops(seed) {
        Ok((_, w, _)) => w,
        Err(e) => return empty_svg(&format!("Band + reads (seed alignment error: {e})")),
    };

    // Band corridor from the seed.
    let unbanded = w == usize::MAX;
    let max_read_len = reads.iter().map(|r| r.len()).max().unwrap_or(0);
    let x_vals: Vec<f64> = (0..n_nodes).map(|i| i as f64).collect();
    let lower: Vec<f64> = x_vals
        .iter()
        .map(|&x| {
            if unbanded {
                0.0
            } else {
                (x - w as f64).max(0.0)
            }
        })
        .collect();
    let upper: Vec<f64> = x_vals
        .iter()
        .map(|&x| {
            if unbanded {
                max_read_len as f64
            } else {
                (x + w as f64).min(max_read_len as f64)
            }
        })
        .collect();

    let band = BandPlot::new(x_vals, lower, upper)
        .with_color("#bbbbbb")
        .with_opacity(0.30)
        .with_legend("band corridor");

    let mut plots: Vec<Plot> = vec![band.into()];

    // One LinePlot per read, coloured from the Tol muted palette.
    let palette = Palette::tol_muted();
    let mut colors = palette.iter();

    for (i, read) in reads.iter().enumerate() {
        let color = colors.next().unwrap_or("#888888");
        match graph.align_read_ops(read) {
            Ok((ops, _, rank_of)) => {
                let mut path: Vec<(f64, f64)> = Vec::new();
                let mut read_pos: usize = 0;
                for op in &ops {
                    match op {
                        AlignOp::Match(node_idx) => {
                            path.push((rank_of[*node_idx] as f64, read_pos as f64));
                            read_pos += 1;
                        }
                        AlignOp::Delete(node_idx) => {
                            path.push((rank_of[*node_idx] as f64, read_pos as f64));
                        }
                        AlignOp::Insert(_) => {
                            read_pos += 1;
                        }
                    }
                }
                // Highlight seed in a thicker line, others thinner.
                let sw = if i == seed_idx { 2.5 } else { 1.0 };
                let line = LinePlot::new()
                    .with_data(path)
                    .with_color(color)
                    .with_stroke_width(sw)
                    .with_legend(format!("read {i}"));
                plots.push(line.into());
            }
            Err(_) => {
                // Alignment failed (BandTooNarrow). Retry unbanded to get the
                // true path, then split it into in-band and out-of-band
                // segments. Insert ops are plotted at the last known graph rank
                // so the path is continuous (no gaps). Segments share their
                // boundary point so in-band and out-of-band connect visually.
                if let Ok((ops, rank_of)) = graph.align_read_ops_unbanded(read) {
                    // Build a flat list of (graph_rank, read_pos, in_band).
                    // Insert ops advance read_pos but keep graph_rank fixed,
                    // producing vertical movement in the plot.
                    let mut tagged: Vec<(f64, f64, bool)> = Vec::new();
                    let mut read_pos: usize = 0;
                    let mut last_gr: usize = 0;

                    for op in &ops {
                        let (gr, advance) = match op {
                            AlignOp::Match(n) => (rank_of[*n], true),
                            AlignOp::Delete(n) => (rank_of[*n], false),
                            AlignOp::Insert(_) => (last_gr, true),
                        };
                        last_gr = gr;
                        let inside = w == usize::MAX || read_pos.abs_diff(gr) <= w;
                        tagged.push((gr as f64, read_pos as f64, inside));
                        if advance {
                            read_pos += 1;
                        }
                    }

                    // Split into contiguous same-band segments. Each segment
                    // boundary point is shared with its neighbour so segments
                    // connect without a gap.
                    let mut in_band: Vec<(f64, f64)> = Vec::new();
                    let mut out_band: Vec<(f64, f64)> = Vec::new();
                    let mut prev_inside: Option<bool> = None;

                    for &(gx, ry, inside) in &tagged {
                        match prev_inside {
                            Some(prev) if prev != inside => {
                                // Transition: append the boundary point to the
                                // segment that's ending so they share it.
                                if inside {
                                    out_band.push((gx, ry)); // close out-band
                                } else {
                                    in_band.push((gx, ry)); // close in-band
                                }
                            }
                            _ => {}
                        }
                        if inside {
                            in_band.push((gx, ry));
                        } else {
                            out_band.push((gx, ry));
                        }
                        prev_inside = Some(inside);
                    }

                    if !in_band.is_empty() {
                        plots.push(
                            LinePlot::new()
                                .with_data(in_band)
                                .with_color(color)
                                .with_stroke_width(1.0)
                                .with_legend(format!("read {i} (in band)"))
                                .into(),
                        );
                    }
                    if !out_band.is_empty() {
                        plots.push(
                            LinePlot::new()
                                .with_data(out_band)
                                .with_color("#cc0000")
                                .with_stroke_width(2.0)
                                .with_dashed()
                                .with_legend(format!("read {i} (out of band)"))
                                .into(),
                        );
                    }
                }
            }
        }
    }

    let title = if unbanded {
        "Band + reads (unbanded)".to_string()
    } else {
        format!("Band + reads (w = {w})")
    };

    let layout = Layout::auto_from_plots(&plots)
        .with_title(title)
        .with_x_label("Graph node (topological rank)")
        .with_y_label("Read position");

    render_to_svg(plots, layout)
}

// ── POA graph network diagram ─────────────────────────────────────────────────

/// SVG network diagram of the POA graph.
///
/// The layout matches the standard POA paper convention: the consensus spine
/// runs horizontally through the centre, and bubble arms branch above or below.
/// All node positions are pinned so the layout algorithm is not used.
///
/// - Spine nodes are blue, sized by coverage.
/// - Off-spine (bubble arm) nodes are grey.
/// - When `read` is `Some`, that sequence is aligned against the graph:
///   Match nodes turn red, Delete nodes turn orange.
/// - Labels show `"{rank}:{base}"` (e.g. `"3:G"`) for graphs with ≤ 80 nodes.
pub fn graph_network_svg(graph: &PoaGraph, read: Option<&[u8]>) -> String {
    graph_network_svg_inner(graph, read, false)
}

/// Same as [`graph_network_svg`] but labels each edge with its read count
/// (edge weight). Useful for illustrating the heaviest-path DP and boundary
/// trim in documentation or debugging.
pub fn graph_network_svg_labeled(graph: &PoaGraph, read: Option<&[u8]>) -> String {
    graph_network_svg_inner(graph, read, true)
}

fn graph_network_svg_inner(
    graph: &PoaGraph,
    read: Option<&[u8]>,
    show_edge_weights: bool,
) -> String {
    use std::collections::HashMap;

    let topology = graph.graph_topology();
    if topology.nodes.is_empty() {
        return empty_svg("POA graph (empty)");
    }

    // ── Read path ─────────────────────────────────────────────────────────
    let mut match_ranks: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut delete_ranks: std::collections::HashSet<usize> = std::collections::HashSet::new();
    if let Some(r) = read {
        if let Ok((ops, _, rank_of)) = graph.align_read_ops(r) {
            for op in &ops {
                match op {
                    AlignOp::Match(ni) => {
                        match_ranks.insert(rank_of[*ni]);
                    }
                    AlignOp::Delete(ni) => {
                        delete_ranks.insert(rank_of[*ni]);
                    }
                    AlignOp::Insert(_) => {}
                }
            }
        }
    }

    // ── Spine index ───────────────────────────────────────────────────────
    let n = topology.nodes.len();
    let spine_len = topology.spine_ranks.len();
    let spine_set: std::collections::HashSet<usize> =
        topology.spine_ranks.iter().copied().collect();
    let spine_idx: HashMap<usize, usize> = topology
        .spine_ranks
        .iter()
        .enumerate()
        .map(|(i, &rank)| (rank, i))
        .collect();

    // ── Pinned layout ─────────────────────────────────────────────────────
    let pad = 0.05_f64;
    let spine_y_val = 0.5_f64;

    let spine_x = |si: usize| -> f64 {
        if spine_len <= 1 {
            0.5
        } else {
            pad + (si as f64 / (spine_len - 1) as f64) * (1.0 - 2.0 * pad)
        }
    };

    let mut positions = vec![(0.5_f64, spine_y_val); n];

    for &rank in &topology.spine_ranks {
        positions[rank] = (spine_x(spine_idx[&rank]), spine_y_val);
    }

    // Group arm nodes into chains: directly connected arm nodes (multi-node
    // insertion arms) must land on the same side of the spine. We identify
    // chains via arm-only edges, then assign one y slot per chain.
    let arm_arm_set: std::collections::HashSet<(usize, usize)> = topology
        .edges
        .iter()
        .filter(|e| !spine_set.contains(&e.from_rank) && !spine_set.contains(&e.to_rank))
        .map(|e| (e.from_rank, e.to_rank))
        .collect();

    // chain_of[rank] = root rank of the arm chain. Nodes are visited in
    // topological order so predecessors are already inserted.
    let mut chain_of: HashMap<usize, usize> = HashMap::new();
    for node in &topology.nodes {
        let rank = node.topo_rank;
        if spine_set.contains(&rank) {
            continue;
        }
        let root = chain_of
            .iter()
            .filter(|&(&pred, _)| arm_arm_set.contains(&(pred, rank)))
            .map(|(_, &c)| c)
            .next()
            .unwrap_or(rank);
        chain_of.insert(rank, root);
    }

    // chains: root → [ranks in topological order].
    let mut chains: HashMap<usize, Vec<usize>> = HashMap::new();
    for node in &topology.nodes {
        let rank = node.topo_rank;
        if spine_set.contains(&rank) {
            continue;
        }
        chains.entry(chain_of[&rank]).or_default().push(rank);
    }

    // Assign a y slot per chain. Process in root-rank order for a deterministic
    // layout. Chains at the same (entry, exit) bracket alternate above/below.
    let mut chain_roots: Vec<usize> = chains.keys().copied().collect();
    chain_roots.sort_unstable();

    let mut arm_slots: HashMap<(Option<usize>, Option<usize>), usize> = HashMap::new();
    let mut chain_y: HashMap<usize, f64> = HashMap::new();

    for &root in &chain_roots {
        let chain_ranks = &chains[&root];
        let first = chain_ranks[0];
        let last = *chain_ranks.last().unwrap();
        let entry = (0..first).rev().find(|&r| spine_set.contains(&r));
        let exit = (last + 1..n).find(|&r| spine_set.contains(&r));
        let slot = arm_slots.entry((entry, exit)).or_insert(0);
        let y = if *slot % 2 == 0 { 0.2 } else { 0.8 };
        *slot += 1;
        chain_y.insert(root, y);
    }

    // Position each arm node. Center the chain on the bracket midpoint and
    // use the spine step as the horizontal increment so nodes in a tight
    // bracket (e.g. two adjacent spine nodes) still have room to breathe.
    let spine_step = if spine_len <= 1 {
        0.0
    } else {
        (1.0 - 2.0 * pad) / (spine_len - 1) as f64
    };

    for &root in &chain_roots {
        let chain_ranks = &chains[&root];
        let k = chain_ranks.len();
        let y = chain_y[&root];
        let first = chain_ranks[0];
        let last = *chain_ranks.last().unwrap();
        let entry = (0..first).rev().find(|&r| spine_set.contains(&r));
        let exit = (last + 1..n).find(|&r| spine_set.contains(&r));
        let x_entry = entry
            .and_then(|r| spine_idx.get(&r).copied())
            .map(spine_x)
            .unwrap_or(pad);
        let x_exit = exit
            .and_then(|r| spine_idx.get(&r).copied())
            .map(spine_x)
            .unwrap_or(1.0 - pad);
        let mid = (x_entry + x_exit) / 2.0;
        for (i, &rank) in chain_ranks.iter().enumerate() {
            let x = mid + (i as f64 - (k - 1) as f64 / 2.0) * spine_step;
            positions[rank] = (x, y);
        }
    }

    // ── Node sizing ───────────────────────────────────────────────────────
    let max_cov = topology
        .nodes
        .iter()
        .map(|nd| nd.coverage)
        .max()
        .unwrap_or(1)
        .max(1) as f64;
    let base_r = 22.0_f64;

    // ── Build NetworkPlot ─────────────────────────────────────────────────
    let show_labels = n <= 80;

    let mut net = NetworkPlot::new()
        .with_directed()
        .with_legend("type")
        .with_layout(NetworkLayout::ForceDirected);

    if show_labels {
        net = net.with_labels_inside();
    }

    for node in &topology.nodes {
        let label = fmt_label(node);
        let (x, y) = positions[node.topo_rank];
        net = net.with_node(label.clone()).with_node_position(label, x, y);
    }

    // Arm nodes that have been weight-annotated (avoids double-labelling when
    // both the entry and exit edge of an arm are iterated).
    let mut annotated_arm_nodes: std::collections::HashSet<usize> =
        std::collections::HashSet::new();
    // (norm_x, above_spine, weight_str) — one entry per arm node, used for
    // SVG post-processing after kuva renders the network.
    let mut arm_labels: Vec<(f64, bool, String)> = Vec::new();

    for edge in &topology.edges {
        let src_rank = edge.from_rank;
        let tgt_rank = edge.to_rank;
        let src = fmt_label(&topology.nodes[src_rank]);
        let tgt = fmt_label(&topology.nodes[tgt_rank]);
        let src_on_spine = spine_set.contains(&src_rank);
        let tgt_on_spine = spine_set.contains(&tgt_rank);

        if src_on_spine != tgt_on_spine {
            let arm_rank = if src_on_spine { tgt_rank } else { src_rank };
            let arm_y = positions[arm_rank].1;
            let curve = if arm_y < 0.5 { -0.25 } else { 0.25 };
            net = net.with_edge_curved(src, tgt, edge.weight as f64, curve);

            if show_edge_weights && annotated_arm_nodes.insert(arm_rank) {
                let (ax, _) = positions[arm_rank];
                arm_labels.push((ax, arm_y < 0.5, edge.weight.to_string()));
            }
        } else if show_edge_weights {
            net = net.with_edge_label(src, tgt, edge.weight as f64, edge.weight.to_string());
        } else {
            net = net.with_edge(src, tgt, edge.weight as f64);
        }
    }

    for node in &topology.nodes {
        let label = fmt_label(node);
        let size = base_r + (node.coverage as f64 / max_cov) * 6.0;
        let (color, group) = if match_ranks.contains(&node.topo_rank) {
            ("#d73027", "read match")
        } else if delete_ranks.contains(&node.topo_rank) {
            ("#fdae61", "read delete")
        } else if spine_set.contains(&node.topo_rank) {
            ("#2166ac", "spine")
        } else {
            ("#aaaaaa", "other")
        };
        net = net
            .with_node_color(&label, color)
            .with_node_size(&label, size)
            .with_node_group(&label, group);
    }

    let plots: Vec<Plot> = vec![net.into()];
    let title = if read.is_some() {
        format!("POA graph ({n} nodes, {} spine) [read overlay]", spine_len)
    } else {
        format!("POA graph ({n} nodes, {} spine)", spine_len)
    };
    let canvas_w = (spine_len as f64) * 120.0 + 200.0;
    let layout = Layout::auto_from_plots(&plots)
        .with_title(title)
        .with_width(canvas_w)
        .with_height(500.0);
    let svg = render_to_svg(plots, layout);
    if show_edge_weights && !arm_labels.is_empty() {
        inject_arm_labels(svg, &arm_labels)
    } else {
        svg
    }
}

fn fmt_label(node: &crate::GraphNodeInfo) -> String {
    format!("{}:{}", node.topo_rank, node.base as char)
}

// ── SVG post-processing for arm edge labels ───────────────────────────────────

/// Extract a single f64 attribute value from an SVG element fragment.
/// Looks for `attr_name` (e.g. `cx="`) and parses up to the closing `"`.
fn svg_attr_f64(fragment: &str, attr: &str) -> Option<f64> {
    let start = fragment.find(attr)? + attr.len();
    let end = fragment[start..].find('"')?;
    fragment[start..start + end].parse().ok()
}

/// Inject weight labels onto curved arm edges by parsing the rendered SVG for
/// arm node circle positions and overlaying `<text>` elements directly.
///
/// `TextAnnotation` uses the chart data-axis coordinate system which does not
/// align with the network renderer's separate pixel mapping, so post-processing
/// is the only reliable way to place text at the correct pixel position.
///
/// `arm_labels`: `(norm_x, above_spine, weight_string)` per arm node.
fn inject_arm_labels(mut svg: String, arm_labels: &[(f64, bool, String)]) -> String {
    // Parse arm node circles: fill="#aaaaaa", radius > 15 px (excludes legend dots).
    let mut circles: Vec<(f64, f64, f64)> = Vec::new(); // (cx, cy, r)
    for chunk in svg.split("<circle ") {
        if !chunk.contains("fill=\"#aaaaaa\"") {
            continue;
        }
        if let (Some(cx), Some(cy), Some(r)) = (
            svg_attr_f64(chunk, "cx=\""),
            svg_attr_f64(chunk, "cy=\""),
            svg_attr_f64(chunk, "r=\""),
        ) {
            if r > 15.0 {
                circles.push((cx, cy, r));
            }
        }
    }
    if circles.is_empty() {
        return svg;
    }

    // Sort circles by cx and arm_labels by norm_x — both increase monotonically
    // in the same direction, so positional pairing is correct.
    circles.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let mut labels: Vec<_> = arm_labels.iter().collect();
    labels.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let mut texts = String::new();
    for ((_norm_x, above, weight), &(cx, cy, r)) in labels.iter().zip(circles.iter()) {
        let ty = if *above { cy - r - 4.0 } else { cy + r + 12.0 };
        texts.push_str(&format!(
            "<text x=\"{cx:.1}\" y=\"{ty:.1}\" text-anchor=\"middle\" \
             font-size=\"11\" font-family=\"sans-serif\" fill=\"#444444\">{weight}</text>\n"
        ));
    }

    if let Some(pos) = svg.rfind("</svg>") {
        svg.insert_str(pos, &texts);
    }
    svg
}

// ── Internal helpers ──────────────────────────────────────────────────────────

fn empty_svg(title: &str) -> String {
    let plots: Vec<Plot> = vec![];
    let layout = Layout::auto_from_plots(&plots).with_title(title);
    render_to_svg(plots, layout)
}

//! Visualisation helpers for POA graphs and consensus output (requires feature `plot`).
//!
//! All functions return an SVG string.  Write it to a file or embed it in HTML.
//!
//! # Example
//!
//! ```rust,no_run
//! use poa_consensus::{PoaGraph, PoaConfig};
//! use poa_consensus::plot::{coverage_svg, graph_stats_svg, edge_weight_histogram_svg};
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

// ── Internal helpers ──────────────────────────────────────────────────────────

fn empty_svg(title: &str) -> String {
    let plots: Vec<Plot> = vec![];
    let layout = Layout::auto_from_plots(&plots).with_title(title);
    render_to_svg(plots, layout)
}

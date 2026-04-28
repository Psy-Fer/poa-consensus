//! Smoke tests for the `plot` feature.  Each function must return valid SVG.
//!
//! Run with:
//!   cargo test --features plot --test plot

use poa_consensus::plot::{
    coverage_svg, edge_weight_histogram_svg, graph_stats_svg, node_coverage_histogram_svg,
};
use poa_consensus::{PoaConfig, PoaGraph};

fn simple_graph() -> PoaGraph {
    let reads: &[&[u8]] = &[
        b"CATCATCATCAT",
        b"CATCATCATCAT",
        b"CATCATCATCAT",
        b"CATCATCGTCAT",
        b"CATCATCATCAT",
    ];
    let mut g = PoaGraph::new(reads[0], PoaConfig::default()).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    g
}

#[test]
fn coverage_svg_is_valid_svg() {
    let g = simple_graph();
    let c = g.consensus().unwrap();
    let svg = coverage_svg(&c);
    assert!(
        svg.contains("<svg"),
        "not SVG: {}",
        &svg[..svg.len().min(100)]
    );
    assert!(svg.contains("</svg>"), "unclosed SVG");
}

#[test]
fn graph_stats_svg_is_valid_svg() {
    let g = simple_graph();
    let stats = g.stats();
    let svg = graph_stats_svg(&stats);
    // graph_stats_svg returns two concatenated SVGs.
    assert!(svg.contains("<svg"), "not SVG");
}

#[test]
fn edge_weight_histogram_svg_is_valid_svg() {
    let g = simple_graph();
    let weights = g.edge_weights();
    assert!(!weights.is_empty(), "no edge weights");
    let svg = edge_weight_histogram_svg(&weights);
    assert!(svg.contains("<svg"), "not SVG");
}

#[test]
fn node_coverage_histogram_svg_is_valid_svg() {
    let g = simple_graph();
    let covs = g.node_coverages();
    assert!(!covs.is_empty(), "no node coverages");
    let svg = node_coverage_histogram_svg(&covs);
    assert!(svg.contains("<svg"), "not SVG");
}

#[test]
fn edge_weight_histogram_empty_graph_does_not_panic() {
    let svg = edge_weight_histogram_svg(&[]);
    assert!(svg.contains("<svg"));
}

#[test]
fn node_coverage_histogram_empty_does_not_panic() {
    let svg = node_coverage_histogram_svg(&[]);
    assert!(svg.contains("<svg"));
}

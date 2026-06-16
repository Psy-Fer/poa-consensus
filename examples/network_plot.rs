//! Quick demo: render the POA graph from a fixture as a network SVG.
//!
//! Run with:
//!   cargo run --example network_plot --features plot

use poa_consensus::plot::graph_network_svg;
use poa_consensus::{PoaConfig, PoaGraph};
use std::fs;

fn main() {
    // ── Two-allele set: SNP bubble at position 11 (A vs G) ────────────────
    let allele_a: &[&[u8]] = &[
        b"GCTAGCTAGCTACTAGCTAGCT",
        b"GCTAGCTAGCTACTAGCTAGCT",
        b"GCTAGCTAGCTACTAGCTAGCT",
        b"GCTAGCTAGCTACTAGCTAGCT",
        b"GCTAGCTAGCTACTAGCTAGCT",
        b"GCTAGCTAGCTACTAGCTAGCT",
        b"GCTAGCTTGCTACTAGCTAGCT", // one noisy read
        b"GCTAGCTAGCTACTAGCTAGCT",
        b"GCTAGCTAGCTACTAGCTAGCT",
        b"GCTAGCTAGCTACTAGCTAGCT",
    ];
    let allele_b: &[&[u8]] = &[
        b"GCTAGCTAGCTGCTAGCTAGCT",
        b"GCTAGCTAGCTGCTAGCTAGCT",
        b"GCTAGCTAGCTGCTAGCTAGCT",
        b"GCTAGCTAGCTGCTAGCTAGCT",
        b"GCTAGCTAGCTGCTAGCTAGCT",
        b"GCTAGCTAGCTGCTAGCTAGCT",
        b"GCTAGCTAGCTGCTAGCTAGCT",
        b"GCTAGCTAGCTGCTAGCTAGCT",
        b"GCTAGCTAGCTGCTAGCTAGCT",
        b"GCTAGCTAGCTGCTAGCTAGCT",
    ];

    let mut reads: Vec<&[u8]> = allele_a.to_vec();
    reads.extend_from_slice(allele_b);

    let config = PoaConfig::default();
    let mut graph = PoaGraph::new(reads[0], config).unwrap();
    for r in &reads[1..] {
        graph.add_read(r).unwrap();
    }

    // Basic graph (spine highlighted in blue)
    let svg = graph_network_svg(&graph, None);
    fs::write("/tmp/poa_network.svg", &svg).unwrap();
    println!("poa_network.svg written ({} bytes)", svg.len());

    // Graph with a specific read's path overlaid
    let probe = b"GCTAGCTAGCTGCTAGCTAGCT"; // allele-B read
    let svg2 = graph_network_svg(&graph, Some(probe));
    fs::write("/tmp/poa_network_with_read.svg", &svg2).unwrap();
    println!("poa_network_with_read.svg written ({} bytes)", svg2.len());

    // ── Noisy single-allele set ────────────────────────────────────────────
    let noisy: &[&[u8]] = &[
        b"GCTAGCTAGCTAGCTA",
        b"GCTAGCTAGCTAGCTA",
        b"GCTAACTAGCTAGCTA",
        b"GCTAGCTAGCTAGCTA",
        b"GCTAGCTAGCTAGCTA",
        b"GCTAGCTAGCTAGCTA",
        b"GCTAGCTAGCTAGCTA",
        b"GCTAGCTAGCTATCTA",
        b"GCTAGCTAGCTAGCTA",
        b"GCTAGCTAGCTAGCTA",
    ];

    let mut g2 = PoaGraph::new(noisy[0], PoaConfig::default()).unwrap();
    for r in &noisy[1..] {
        g2.add_read(r).unwrap();
    }

    let svg3 = graph_network_svg(&g2, None);
    fs::write("/tmp/poa_network_noisy.svg", &svg3).unwrap();
    println!("poa_network_noisy.svg written ({} bytes)", svg3.len());
}

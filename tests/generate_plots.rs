//! Tests that generate SVG plots into `test_plots/`.
//!
//! Run with:
//!   cargo test --features plot --test generate_plots -- --nocapture
//!
//! Inspect the output in `test_plots/` to visually verify plot quality.

#![cfg(feature = "plot")]

use std::fs;
use std::path::PathBuf;

use poa_consensus::plot::{
    alignment_density_svg, band_svg, band_with_reads_svg, coverage_svg, edge_weight_histogram_svg,
    graph_stats_svg, node_coverage_histogram_svg,
};
use poa_consensus::{AlignmentMode, PoaConfig, PoaGraph, consensus, consensus_multi};

fn out_dir() -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test_plots");
    fs::create_dir_all(&dir).expect("create test_plots/");
    dir
}

fn save(name: &str, svg: &str) {
    let path = out_dir().join(name);
    fs::write(&path, svg).unwrap_or_else(|e| panic!("write {}: {e}", path.display()));
    eprintln!("wrote {}", path.display());
}

fn default_config() -> PoaConfig {
    PoaConfig::default()
}

// ── simple fixture ────────────────────────────────────────────────────────────

fn simple_reads() -> Vec<Vec<u8>> {
    vec![
        b"ACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        b"ACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        b"ACGTACGTACTTACGTACGTACGTACGT".to_vec(),
        b"ACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        b"ACGTACGTACGTACGTACGTACGGACGT".to_vec(),
        b"ACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        b"ACGTAGGTACGTACGTACGTACGTACGT".to_vec(),
        b"ACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        b"ACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        b"ACGTACGTACGTACGTACGTACGTACGT".to_vec(),
    ]
}

#[test]
fn plot_simple_coverage() {
    let reads_owned = simple_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let c = consensus(&reads, 0, &default_config()).unwrap();
    save("simple_coverage.svg", &coverage_svg(&c));
    assert!(!c.coverage.is_empty());
}

#[test]
fn plot_simple_graph_stats() {
    let reads_owned = simple_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let mut g = PoaGraph::new(reads[0], default_config()).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    save("simple_graph_stats.svg", &graph_stats_svg(&g.stats()));
}

#[test]
fn plot_simple_edge_weights() {
    let reads_owned = simple_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let mut g = PoaGraph::new(reads[0], default_config()).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    let weights = g.edge_weights();
    save(
        "simple_edge_weights.svg",
        &edge_weight_histogram_svg(&weights),
    );
}

#[test]
fn plot_simple_node_coverage() {
    let reads_owned = simple_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let mut g = PoaGraph::new(reads[0], default_config()).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    let covs = g.node_coverages();
    save(
        "simple_node_coverage.svg",
        &node_coverage_histogram_svg(&covs),
    );
}

// ── noisy fixture ─────────────────────────────────────────────────────────────

fn noisy_reads() -> Vec<Vec<u8>> {
    vec![
        b"GCTAGCTAGCTAGCTA".to_vec(),
        b"GCTAGCTAGCTAGCTA".to_vec(),
        b"GCTAACTAGCTAGCTA".to_vec(),
        b"GCTAGCTAGCTAGCTA".to_vec(),
        b"GCTAGCTAGCTAGCTA".to_vec(),
        b"GCTAGCTAGCTAGCTA".to_vec(),
        b"GCTAGCTAGCTAGCTA".to_vec(),
        b"GCTAGCTAGCTATCTA".to_vec(),
        b"GCTAGCTAGCTAGCTA".to_vec(),
        b"GCTAGCTAGCTAGCTA".to_vec(),
    ]
}

#[test]
fn plot_noisy_coverage() {
    let reads_owned = noisy_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let c = consensus(&reads, 0, &default_config()).unwrap();
    save("noisy_coverage.svg", &coverage_svg(&c));
}

#[test]
fn plot_noisy_edge_weights() {
    let reads_owned = noisy_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let mut g = PoaGraph::new(reads[0], default_config()).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    save(
        "noisy_edge_weights.svg",
        &edge_weight_histogram_svg(&g.edge_weights()),
    );
}

// ── two-allele fixture ────────────────────────────────────────────────────────

fn two_allele_reads() -> Vec<Vec<u8>> {
    // 10 reads per allele; alleles differ at position 11 (A vs G)
    let mut v: Vec<Vec<u8>> = Vec::new();
    for _ in 0..9 {
        v.push(b"GCTAGCTAGCTACTAGCTAGCT".to_vec());
    }
    v.push(b"GCTAGCTTGCTACTAGCTAGCT".to_vec()); // SNP noise in allele 1
    for _ in 0..10 {
        v.push(b"GCTAGCTAGCTGCTAGCTAGCT".to_vec());
    }
    v
}

#[test]
fn plot_two_allele_per_allele_coverage() {
    let reads_owned = two_allele_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let alleles = consensus_multi(&reads, 0, &default_config()).unwrap();
    assert_eq!(alleles.len(), 2, "expected 2 alleles");
    for (i, allele) in alleles.iter().enumerate() {
        save(
            &format!("two_allele_allele{}_coverage.svg", i + 1),
            &coverage_svg(allele),
        );
    }
}

#[test]
fn plot_two_allele_graph_stats() {
    let reads_owned = two_allele_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let mut g = PoaGraph::new(reads[0], default_config()).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    save("two_allele_graph_stats.svg", &graph_stats_svg(&g.stats()));
    save(
        "two_allele_node_coverage.svg",
        &node_coverage_histogram_svg(&g.node_coverages()),
    );
    save(
        "two_allele_edge_weights.svg",
        &edge_weight_histogram_svg(&g.edge_weights()),
    );
}

// ── banded vs unbanded — longer reads ────────────────────────────────────────

fn long_reads() -> Vec<Vec<u8>> {
    // 200 bp reads: ground truth + 5% SNP noise
    let gt: Vec<u8> = b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\
                         ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\
                         ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\
                         ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        .to_vec();
    let mut reads = vec![gt.clone(); 8];
    // introduce one SNP per noisy read
    reads[2][30] = b'G';
    reads[5][80] = b'C';
    reads[7][150] = b'A';
    reads
}

#[test]
fn plot_long_reads_coverage() {
    let reads_owned = long_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let c = consensus(&reads, 0, &default_config()).unwrap();
    save("long_reads_coverage.svg", &coverage_svg(&c));
}

#[test]
fn plot_long_reads_banded_coverage() {
    let reads_owned = long_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let config = PoaConfig {
        band_width: 50,
        ..PoaConfig::default()
    };
    let c = consensus(&reads, 0, &config).unwrap();
    save("long_reads_banded_coverage.svg", &coverage_svg(&c));
}

// ── semi-global — partial reads ───────────────────────────────────────────────

fn partial_reads() -> Vec<Vec<u8>> {
    vec![
        b"ATCGATCGATCGATCG".to_vec(), // full
        b"ATCGATCGATCGATCG".to_vec(),
        b"ATCGATCGATCGATCG".to_vec(),
        b"ATCGATCGATCGATCG".to_vec(),
        b"ATCGATCGATCGATCG".to_vec(),
        b"ATCGATCGATCG".to_vec(), // partial left
        b"ATCGATCGATCG".to_vec(),
        b"TCGATCGATCGATCG".to_vec(), // partial right
        b"CGATCGATCGATCG".to_vec(),
        b"GATCGATCGATCGATCG".to_vec(),
    ]
}

#[test]
fn plot_partial_reads_semi_global_coverage() {
    let reads_owned = partial_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let config = PoaConfig {
        alignment_mode: AlignmentMode::SemiGlobal,
        ..PoaConfig::default()
    };
    let c = consensus(&reads, 0, &config).unwrap();
    save("partial_semi_global_coverage.svg", &coverage_svg(&c));
    // Coverage should drop at the boundaries where partial reads don't reach.
    eprintln!(
        "partial semi-global consensus: {}",
        String::from_utf8_lossy(&c.sequence)
    );
}

// ── alignment density (2D histogram) ─────────────────────────────────────────

/// Clean single-allele input — expect a tight diagonal.
#[test]
fn plot_alignment_density_simple() {
    let reads_owned = simple_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let mut g = PoaGraph::new(reads[0], default_config()).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    save(
        "alignment_density_simple.svg",
        &alignment_density_svg(&g, &reads),
    );
}

/// Two-allele input — expect the diagonal to split at the bubble.
#[test]
fn plot_alignment_density_two_allele() {
    let reads_owned = two_allele_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let mut g = PoaGraph::new(reads[0], default_config()).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    save(
        "alignment_density_two_allele.svg",
        &alignment_density_svg(&g, &reads),
    );
}

/// Noisy reads — diagonal with scattered off-diagonal points from SNPs.
#[test]
fn plot_alignment_density_noisy() {
    let reads_owned = noisy_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let mut g = PoaGraph::new(reads[0], default_config()).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    save(
        "alignment_density_noisy.svg",
        &alignment_density_svg(&g, &reads),
    );
}

/// Partial reads with length variation — shows diagonal with boundary dropout.
#[test]
fn plot_alignment_density_partial() {
    let reads_owned = partial_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let config = PoaConfig {
        alignment_mode: AlignmentMode::SemiGlobal,
        ..PoaConfig::default()
    };
    let mut g = PoaGraph::new(reads[0], config).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    save(
        "alignment_density_partial.svg",
        &alignment_density_svg(&g, &reads),
    );
}

// ── band corridor ─────────────────────────────────────────────────────────────

/// Unbanded: corridor fills the full plot area; path runs diagonally.
#[test]
fn plot_band_unbanded() {
    let reads_owned = simple_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let mut g = PoaGraph::new(reads[0], default_config()).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    save("band_unbanded.svg", &band_svg(&g, reads[0]));
}

/// Narrow band: corridor is visibly constrained around the diagonal.
#[test]
fn plot_band_narrow() {
    let reads_owned = simple_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let config = PoaConfig {
        band_width: 5,
        ..PoaConfig::default()
    };
    let mut g = PoaGraph::new(reads[0], config).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    save("band_narrow.svg", &band_svg(&g, reads[0]));
}

/// Adaptive band: corridor width grows with read length (abPOA formula).
#[test]
fn plot_band_adaptive_long_reads() {
    let reads_owned = long_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let config = PoaConfig {
        adaptive_band: true,
        adaptive_band_b: 10,
        adaptive_band_f: 0.05,
        ..PoaConfig::default()
    };
    let mut g = PoaGraph::new(reads[0], config).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    // Show band for the noisy read (index 2) so path drifts slightly.
    save("band_adaptive_long.svg", &band_svg(&g, reads[2]));
}

// ── band + reads overlaid ─────────────────────────────────────────────────────

/// Clean reads, generous band — all paths cluster tightly on the diagonal.
#[test]
fn plot_band_with_reads_clean() {
    let reads_owned = simple_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let config = PoaConfig {
        band_width: 8,
        ..PoaConfig::default()
    };
    let mut g = PoaGraph::new(reads[0], config).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    save(
        "band_with_reads_clean.svg",
        &band_with_reads_svg(&g, &reads, 0),
    );
}

/// Noisy reads, tight band — SNP reads drift but stay inside; seed bold.
#[test]
fn plot_band_with_reads_noisy_tight() {
    let reads_owned = noisy_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let config = PoaConfig {
        band_width: 4,
        ..PoaConfig::default()
    };
    let mut g = PoaGraph::new(reads[0], config).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    save(
        "band_with_reads_noisy_tight.svg",
        &band_with_reads_svg(&g, &reads, 0),
    );
}

/// Two-allele reads, moderate band — two path clusters diverge at the bubble.
#[test]
fn plot_band_with_reads_two_allele() {
    let reads_owned = two_allele_reads();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let config = PoaConfig {
        band_width: 6,
        ..PoaConfig::default()
    };
    let mut g = PoaGraph::new(reads[0], config).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    save(
        "band_with_reads_two_allele.svg",
        &band_with_reads_svg(&g, &reads, 0),
    );
}

/// Shows three groups of reads against a narrow band:
///   - majority (6 reads, 24 bp): paths sit on the diagonal
///   - near-edge (1 read, 27 bp, 3-bp insertion): path drifts to the band edge
///     but succeeds — visible as a line touching the corridor boundary
///   - failing (1 read, 30 bp, 6-bp insertion): exceeds the band → red dashed line
///
/// The test intentionally builds the graph *without* the failing read so the
/// graph topology stays clean; `band_with_reads_svg` then tries to align all
/// reads including the failing one and shows the failure visually.
#[test]
fn plot_band_with_reads_outlier_deletion() {
    //          |-------- 24 bp majority --------|
    //          |---------- 24 bp majority ----------|
    let majority: &[u8] = b"GCTAGCTAGCTAGCTAGCTAGCTA";
    //          |--- 13 bp ---|-- 3 insert --|-- 11 bp ---|  = 27 bp
    let near_edge: &[u8] = b"GCTAGCTAGCTAGAAACTAGCTAGCTA";
    //          |--- 13 bp ---|---- 6 insert ----|-- 11 bp ---|  = 30 bp
    let failing: &[u8] = b"GCTAGCTAGCTAGAAAAAACTAGCTAGCTA";

    let reads_owned: Vec<Vec<u8>> = std::iter::repeat(majority.to_vec())
        .take(6)
        .chain(std::iter::once(near_edge.to_vec()))
        .chain(std::iter::once(failing.to_vec()))
        .collect();
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();

    // Band of 5: 3-bp insertion drifts 3 cells (fits); 6-bp drifts 6 cells (fails).
    let config = PoaConfig {
        band_width: 5,
        ..PoaConfig::default()
    };
    // Build graph from majority only — keeps it at 24 nodes so the band
    // constraint is tight.  Neither outlier is added: they are shown only
    // via align_read_ops in band_with_reads_svg, so the graph topology stays
    // clean and the band boundaries are unambiguous.
    let mut g = PoaGraph::new(reads[0], config).unwrap();
    for r in &reads[1..6] {
        g.add_read(r).unwrap();
    }

    save(
        "band_with_reads_outlier.svg",
        &band_with_reads_svg(&g, &reads, 0),
    );
}

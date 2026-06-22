//! Demo: render POA graph fixtures as network SVGs.
//!
//! Run with:
//!   cargo run --example network_plot --features plot
//!
//! Writes to /tmp/ for quick inspection, and to docs/src/diagrams/ for the
//! documentation website.

use poa_consensus::plot::{graph_network_svg, graph_network_svg_labeled};
use poa_consensus::{PoaConfig, PoaGraph};
use std::fs;
use std::path::Path;

fn write(path: &str, svg: &str) {
    if let Some(parent) = Path::new(path).parent() {
        fs::create_dir_all(parent).ok();
    }
    fs::write(path, svg).unwrap();
    println!("{} ({} bytes)", path, svg.len());
}

fn main() {
    // ── Linear consensus (all reads identical, no variants) ───────────────
    // Uses 5 reads so min_cov = ceil(5 * 0.5) = 3, all nodes have coverage 5.
    // Shows the simplest possible graph: a single linear spine, no arms.
    let linear: &[&[u8]] = &[
        b"CATCATCAT",
        b"CATCATCAT",
        b"CATCATCAT",
        b"CATCATCAT",
        b"CATCATCAT",
    ];
    let mut g_linear = PoaGraph::new(linear[0], PoaConfig::default()).unwrap();
    for r in &linear[1..] {
        g_linear.add_read(r).unwrap();
    }
    let svg = graph_network_svg(&g_linear, None);
    write("/tmp/poa_linear.svg", &svg);
    write("docs/src/diagrams/poa_linear.svg", &svg);

    // ── Heaviest-path edge weights: 9 majority + 1 minority ──────────────
    // One read has a mismatch at position 3 (A→G), creating a single-node arm.
    // Spine edges carry weight 10; the arm's entry/exit edges carry weight 1.
    // Edge labels make the heaviest path obvious: follow weight 10, not weight 1.
    let hp_reads: &[&[u8]] = &[
        b"CATCAT", b"CATCAT", b"CATCAT", b"CATCAT", b"CATCAT", b"CATCAT", b"CATCAT", b"CATCAT",
        b"CATCAT", b"CATGAT", // mismatch at pos 3: A→G
    ];
    let mut g_hp = PoaGraph::new(hp_reads[0], PoaConfig::default()).unwrap();
    for r in &hp_reads[1..] {
        g_hp.add_read(r).unwrap();
    }
    let svg_hp = graph_network_svg_labeled(&g_hp, None);
    write("/tmp/poa_heaviest_path.svg", &svg_hp);
    write("docs/src/diagrams/poa_heaviest_path.svg", &svg_hp);

    // ── Boundary trim: full-span vs partial reads ─────────────────────────
    // 3 reads span the full GGGCATCATGGG (12 bp).
    // 7 reads cover only the interior CATCAT (6 bp).
    // With semi-global alignment (default), partial reads take free terminal
    // gaps. Boundary nodes (the GGG flanks, coverage 3) fall below
    // min_cov = ceil(10 * 0.5) = 5 and are trimmed; interior nodes survive.
    let full_span: &[&[u8]] = &[b"GGGCATCATGGG", b"GGGCATCATGGG", b"GGGCATCATGGG"];
    let partial: &[&[u8]] = &[
        b"CATCAT", b"CATCAT", b"CATCAT", b"CATCAT", b"CATCAT", b"CATCAT", b"CATCAT",
    ];
    let mut all_reads: Vec<&[u8]> = full_span.to_vec();
    all_reads.extend_from_slice(partial);
    let mut g_trim = PoaGraph::new(all_reads[0], PoaConfig::default()).unwrap();
    for r in &all_reads[1..] {
        g_trim.add_read(r).unwrap();
    }
    let svg_trim = graph_network_svg_labeled(&g_trim, None);
    write("/tmp/poa_boundary_trim.svg", &svg_trim);
    write("docs/src/diagrams/poa_boundary_trim.svg", &svg_trim);

    // ── CAG repeat: normal vs pathogenic allele lengths ───────────────────
    // The heaviest-path algorithm always selects the path with the highest
    // cumulative edge weight. In a mixed two-allele read set, the longer allele
    // accumulates more edge weight and wins the spine calculation. To clearly
    // show both allele lengths, we render them as separate single-allele graphs.
    //
    //   Normal:     (CAG)×4 = 12 bp   (e.g. HTT healthy range)
    //   Pathogenic: (CAG)×7 = 21 bp   (pathogenic; real disease: >36 units)
    fn build_single_allele(reads: &[&[u8]]) -> PoaGraph {
        let mut g = PoaGraph::new(reads[0], PoaConfig::default()).unwrap();
        for r in &reads[1..] {
            g.add_read(r).unwrap();
        }
        g
    }

    let cag4: &[&[u8]] = &[
        b"CAGCAGCAGCAG",
        b"CAGCAGCAGCAG",
        b"CAGCAGCAGCAG",
        b"CAGCAGCAGCAG",
        b"CAGCAGCAGCAG",
        b"CAGCAGCAGCAG",
        b"CAGCAGCAGCAG",
        b"CAGCAGCAGCAG",
        b"CAGCAGCAGCAG",
        b"CAGCAGCAGCAG",
    ];
    let cag7: &[&[u8]] = &[
        b"CAGCAGCAGCAGCAGCAGCAG",
        b"CAGCAGCAGCAGCAGCAGCAG",
        b"CAGCAGCAGCAGCAGCAGCAG",
        b"CAGCAGCAGCAGCAGCAGCAG",
        b"CAGCAGCAGCAGCAGCAGCAG",
        b"CAGCAGCAGCAGCAGCAGCAG",
        b"CAGCAGCAGCAGCAGCAGCAG",
        b"CAGCAGCAGCAGCAGCAGCAG",
        b"CAGCAGCAGCAGCAGCAGCAG",
        b"CAGCAGCAGCAGCAGCAGCAG",
    ];
    let svg_cag_normal = graph_network_svg(&build_single_allele(cag4), None);
    write("/tmp/poa_cag_normal.svg", &svg_cag_normal);
    write("docs/src/diagrams/poa_cag_normal.svg", &svg_cag_normal);

    let svg_cag_expanded = graph_network_svg(&build_single_allele(cag7), None);
    write("/tmp/poa_cag_expanded.svg", &svg_cag_expanded);
    write("docs/src/diagrams/poa_cag_expanded.svg", &svg_cag_expanded);

    // Legacy combined name (kept so docs references don't break if referenced elsewhere)
    write("docs/src/diagrams/poa_cag_repeat.svg", &svg_cag_normal);

    // ── GAA repeat: normal vs FRDA-expanded allele ────────────────────────
    // Models the FXN intron-1 GAA repeat (Friedreich's ataxia).
    //   Normal:     (GAA)×4 = 12 bp
    //   Pathogenic: (GAA)×8 = 24 bp   (disease range: >66 units; scaled for clarity)
    let gaa4: &[&[u8]] = &[
        b"GAAGAAGAAGAA",
        b"GAAGAAGAAGAA",
        b"GAAGAAGAAGAA",
        b"GAAGAAGAAGAA",
        b"GAAGAAGAAGAA",
        b"GAAGAAGAAGAA",
        b"GAAGAAGAAGAA",
        b"GAAGAAGAAGAA",
        b"GAAGAAGAAGAA",
        b"GAAGAAGAAGAA",
    ];
    let gaa8: &[&[u8]] = &[
        b"GAAGAAGAAGAAGAAGAAGAAGAA",
        b"GAAGAAGAAGAAGAAGAAGAAGAA",
        b"GAAGAAGAAGAAGAAGAAGAAGAA",
        b"GAAGAAGAAGAAGAAGAAGAAGAA",
        b"GAAGAAGAAGAAGAAGAAGAAGAA",
        b"GAAGAAGAAGAAGAAGAAGAAGAA",
        b"GAAGAAGAAGAAGAAGAAGAAGAA",
        b"GAAGAAGAAGAAGAAGAAGAAGAA",
        b"GAAGAAGAAGAAGAAGAAGAAGAA",
        b"GAAGAAGAAGAAGAAGAAGAAGAA",
    ];
    let svg_gaa_normal = graph_network_svg(&build_single_allele(gaa4), None);
    write("/tmp/poa_gaa_normal.svg", &svg_gaa_normal);
    write("docs/src/diagrams/poa_gaa_normal.svg", &svg_gaa_normal);

    let svg_gaa_expanded = graph_network_svg(&build_single_allele(gaa8), None);
    write("/tmp/poa_gaa_expanded.svg", &svg_gaa_expanded);
    write("docs/src/diagrams/poa_gaa_expanded.svg", &svg_gaa_expanded);

    write("docs/src/diagrams/poa_gaa_repeat.svg", &svg_gaa_normal);

    // ── Two-allele set: SNP bubble ─────────────────────────────────────────
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

    let mut graph = PoaGraph::new(reads[0], PoaConfig::default()).unwrap();
    for r in &reads[1..] {
        graph.add_read(r).unwrap();
    }

    let svg_net = graph_network_svg(&graph, None);
    write("/tmp/poa_network.svg", &svg_net);
    write("docs/src/diagrams/poa_network.svg", &svg_net);

    let probe = b"GCTAGCTAGCTGCTAGCTAGCT";
    let svg_overlay = graph_network_svg(&graph, Some(probe));
    write("/tmp/poa_network_with_read.svg", &svg_overlay);
    write("docs/src/diagrams/poa_network_with_read.svg", &svg_overlay);

    // ── Large bubble: 2-node insertion arm ────────────────────────────────
    // 15 reference reads + 5 variant reads with a 2-base GG insertion at
    // position 4. With only 2 arm nodes, the arm path (3 edges × weight 5 =
    // score 12) scores below the direct spine edge (weight 15 = score 14), so
    // the reference stays on the spine and the insertion arm stays grey.
    // min_cov = ceil(20 * 0.5) = 10; arm coverage 5 < 10 (below threshold).
    let ref_reads: &[&[u8]] = &[
        b"ACGTACGTACGT",
        b"ACGTACGTACGT",
        b"ACGTACGTACGT",
        b"ACGTACGTACGT",
        b"ACGTACGTACGT",
        b"ACGTACGTACGT",
        b"ACGTACGTACGT",
        b"ACGTACGTACGT",
        b"ACGTACGTACGT",
        b"ACGTACGTACGT",
        b"ACGTACGTACGT",
        b"ACGTACGTACGT",
        b"ACGTACGTACGT",
        b"ACGTACGTACGT",
        b"ACGTACGTACGT",
    ];
    // 2-base GG insertion between T(3) and A(4): ACGT | GG | ACGT ACGT
    let ins_reads: &[&[u8]] = &[
        b"ACGTGGACGTACGT",
        b"ACGTGGACGTACGT",
        b"ACGTGGACGTACGT",
        b"ACGTGGACGTACGT",
        b"ACGTGGACGTACGT",
    ];
    let mut all_bubble: Vec<&[u8]> = ref_reads.to_vec();
    all_bubble.extend_from_slice(ins_reads);
    let mut g_bubble = PoaGraph::new(all_bubble[0], PoaConfig::default()).unwrap();
    for r in &all_bubble[1..] {
        g_bubble.add_read(r).unwrap();
    }
    let svg_bubble = graph_network_svg(&g_bubble, None);
    write("/tmp/poa_large_bubble.svg", &svg_bubble);
    write("docs/src/diagrams/poa_large_bubble.svg", &svg_bubble);

    // ── Noisy single-allele: multiple insertion errors at different positions
    // 6 clean reads + 4 error reads each with a single-base (or two-base)
    // insertion at a distinct position. This produces 4 grey arm nodes
    // scattered across the spine rather than 2 clustered errors. Each arm
    // has coverage 1, well below min_cov = ceil(10 * 0.5) = 5.
    //
    //   Error 1: T inserted after C(1)   -> 1-node arm at position 1-2
    //   Error 2: G inserted after G(4)   -> 1-node arm at position 4-5
    //   Error 3: A inserted after A(7)   -> 1-node arm at position 7-8
    //   Error 4: CC inserted after A(11) -> 2-node arm at position 11-12
    let noisy: &[&[u8]] = &[
        b"GCTAGCTAGCTAGCTA", // clean ×6
        b"GCTAGCTAGCTAGCTA",
        b"GCTAGCTAGCTAGCTA",
        b"GCTAGCTAGCTAGCTA",
        b"GCTAGCTAGCTAGCTA",
        b"GCTAGCTAGCTAGCTA",
        b"GCTTAGCTAGCTAGCTA",  // T inserted after C(1)
        b"GCTAGGCTAGCTAGCTA",  // G inserted after G(4)
        b"GCTAGCTAAGCTAGCTA",  // A inserted after A(7)
        b"GCTAGCTAGCTACCGCTA", // CC inserted after A(11) -> 2-node arm
    ];
    let mut g_noisy = PoaGraph::new(noisy[0], PoaConfig::default()).unwrap();
    for r in &noisy[1..] {
        g_noisy.add_read(r).unwrap();
    }
    let svg_noisy = graph_network_svg(&g_noisy, None);
    write("/tmp/poa_network_noisy.svg", &svg_noisy);
    write("docs/src/diagrams/poa_network_noisy.svg", &svg_noisy);

    // ── RFC1-style mixed allele: AAAAG normal + AAAGT expanded ────────────
    // Models an unphased mixed read set from the RFC1 CANVAS locus.
    //   Normal (75%):     (AAAAG)×4 = 20 bp  -- majority allele
    //   Pathogenic (25%): (AAAGT)×6 = 30 bp  -- minority expanded allele
    //
    // The two motifs are 2 substitutions per 5-mer apart; the aligner prefers
    // 8 mismatches (score -8) over Insert+Delete pairs (score -6 each x8=-48).
    // After aligning the shared 20-bp region, the pathogenic reads still have
    // 10 bases (2 extra AAAGT units) to place: these become 10 Insert arm nodes
    // extending the spine. With semi-global mode (default), the 15 normal reads
    // take free terminal gaps for the 10 extension nodes, leaving them with
    // coverage 5 -- below min_cov = ceil(20 * 0.5) = 10.
    // The heaviest path extends through the arm (accumulated edge score beats
    // early termination), so all 30 nodes are on the spine; the coverage drop
    // is visible in the smaller node sizes at positions 20-29.
    let rfc1_normal: &[&[u8]] = &[
        b"AAAAGAAAAGAAAAGAAAAG",
        b"AAAAGAAAAGAAAAGAAAAG",
        b"AAAAGAAAAGAAAAGAAAAG",
        b"AAAAGAAAAGAAAAGAAAAG",
        b"AAAAGAAAAGAAAAGAAAAG",
        b"AAAAGAAAAGAAAAGAAAAG",
        b"AAAAGAAAAGAAAAGAAAAG",
        b"AAAAGAAAAGAAAAGAAAAG",
        b"AAAAGAAAAGAAAAGAAAAG",
        b"AAAAGAAAAGAAAAGAAAAG",
        b"AAAAGAAAAGAAAAGAAAAG",
        b"AAAAGAAAAGAAAAGAAAAG",
        b"AAAAGAAAAGAAAAGAAAAG",
        b"AAAAGAAAAGAAAAGAAAAG",
        b"AAAAGAAAAGAAAAGAAAAG",
    ];
    let rfc1_patho: &[&[u8]] = &[
        b"AAAGTAAAGTAAAGTAAAGTAAAGTAAAGT",
        b"AAAGTAAAGTAAAGTAAAGTAAAGTAAAGT",
        b"AAAGTAAAGTAAAGTAAAGTAAAGTAAAGT",
        b"AAAGTAAAGTAAAGTAAAGTAAAGTAAAGT",
        b"AAAGTAAAGTAAAGTAAAGTAAAGTAAAGT",
    ];
    let mut rfc1_reads: Vec<&[u8]> = rfc1_normal.to_vec();
    rfc1_reads.extend_from_slice(rfc1_patho);
    let mut g_rfc1 = PoaGraph::new(rfc1_reads[0], PoaConfig::default()).unwrap();
    for r in &rfc1_reads[1..] {
        g_rfc1.add_read(r).unwrap();
    }
    let svg_rfc1 = graph_network_svg(&g_rfc1, None);
    write("/tmp/poa_rfc1_mixed.svg", &svg_rfc1);
    write("docs/src/diagrams/poa_rfc1_mixed.svg", &svg_rfc1);

    // ── Deletion fixture ───────────────────────────────────────────────────
    let with_del: &[&[u8]] = &[
        b"GCTAGCTAGCTAGCTAGCT",
        b"GCTAGCTAGCTAGCTAGCT",
        b"GCTAGCTAGCTAGCTAGCT",
        b"GCTAGCTAGCTAGCTAGCT",
        b"GCTAGCTAGCTAGCTAGCT",
        b"GCTAGCTAGCTAGCTAGCT",
        b"GCTAGCTAGCTAGCTAGCT",
        b"GCTAGCTAGCTAGCTAGCT",
        b"GCTAGCTAGCTAGCTAGCT",
        b"GCTAGCTAGCTAGCTAGCT",
        b"GCTAGCTAGCTAGCT",
        b"GCTAGCTAGCTAGCT",
        b"GCTAGCTAGCTAGCT",
    ];
    let mut g_del = PoaGraph::new(with_del[0], PoaConfig::default()).unwrap();
    for r in &with_del[1..] {
        g_del.add_read(r).unwrap();
    }
    let svg_del = graph_network_svg(&g_del, None);
    write("/tmp/poa_network_deletion.svg", &svg_del);
    write("docs/src/diagrams/poa_network_deletion.svg", &svg_del);

    let svg_del_ov = graph_network_svg(&g_del, Some(b"GCTAGCTAGCTAGCT"));
    write("/tmp/poa_network_deletion_overlay.svg", &svg_del_ov);
    write(
        "docs/src/diagrams/poa_network_deletion_overlay.svg",
        &svg_del_ov,
    );
}

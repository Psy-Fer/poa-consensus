//! Accuracy tests for POA consensus under simulated error models.
//!
//! Tests measure how well the library recovers a known ground-truth sequence
//! across error rates, depths, band widths, and repeat structures.
//!
//! Run with:
//!   cargo test --test accuracy -- --nocapture

use poa_consensus::{PoaConfig, consensus, consensus_multi};

// ── Simulator ─────────────────────────────────────────────────────────────────

struct Sim {
    sub_rate: f64,
    ins_rate: f64,
    del_rate: f64,
    seed: u64,
}

impl Sim {
    fn new(sub_rate: f64, ins_rate: f64, del_rate: f64, seed: u64) -> Self {
        Self {
            sub_rate,
            ins_rate,
            del_rate,
            seed,
        }
    }

    fn simulate(&self, template: &[u8], n_reads: usize) -> Vec<Vec<u8>> {
        let mut state = self.seed;
        let mut reads = Vec::with_capacity(n_reads);
        for _ in 0..n_reads {
            let mut read = Vec::with_capacity(template.len());
            for &base in template {
                if rand_f64(&mut state) < self.ins_rate {
                    read.push(random_base(&mut state));
                }
                if rand_f64(&mut state) < self.del_rate {
                    continue;
                }
                if rand_f64(&mut state) < self.sub_rate {
                    read.push(random_base_not(base, &mut state));
                } else {
                    read.push(base);
                }
            }
            if !read.is_empty() {
                reads.push(read);
            }
        }
        reads
    }
}

fn xorshift(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

fn rand_f64(state: &mut u64) -> f64 {
    xorshift(state) as f64 / u64::MAX as f64
}

fn random_base(state: &mut u64) -> u8 {
    b"ACGT"[(xorshift(state) % 4) as usize]
}

fn random_base_not(exclude: u8, state: &mut u64) -> u8 {
    let opts: [u8; 3] = match exclude {
        b'A' => [b'C', b'G', b'T'],
        b'C' => [b'A', b'G', b'T'],
        b'G' => [b'A', b'C', b'T'],
        _ => [b'A', b'C', b'G'],
    };
    opts[(xorshift(state) % 3) as usize]
}

// ── Metrics ───────────────────────────────────────────────────────────────────

fn levenshtein(a: &[u8], b: &[u8]) -> usize {
    let (m, n) = (a.len(), b.len());
    let mut prev: Vec<usize> = (0..=n).collect();
    let mut curr = vec![0usize; n + 1];
    for i in 1..=m {
        curr[0] = i;
        for j in 1..=n {
            curr[j] = if a[i - 1] == b[j - 1] {
                prev[j - 1]
            } else {
                1 + prev[j - 1].min(prev[j]).min(curr[j - 1])
            };
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[n]
}

fn identity(truth: &[u8], got: &[u8]) -> f64 {
    let denom = truth.len().max(got.len());
    if denom == 0 {
        return 1.0;
    }
    1.0 - levenshtein(truth, got) as f64 / denom as f64
}

fn refs(reads: &[Vec<u8>]) -> Vec<&[u8]> {
    reads.iter().map(|r| r.as_slice()).collect()
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn repeat(unit: &[u8], n: usize) -> Vec<u8> {
    unit.iter().cycle().take(unit.len() * n).copied().collect()
}

// ── Accuracy tests (all use band_width=0, i.e. unbanded) ─────────────────────
//
// Banded DP accuracy degrades as reads are added because insert operations
// create new graph nodes, shifting topological ranks away from the read-position
// diagonal.  Unbanded DP isolates the POA algorithm from band-width effects.
// The band width tradeoff is covered separately in `band_width_accuracy_tradeoff`.

#[test]
fn zero_error_exact_consensus() {
    let template = b"GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA";
    let reads = Sim::new(0.0, 0.0, 0.0, 1).simulate(template, 10);
    let result = consensus(&refs(&reads), 0, &PoaConfig::default()).unwrap();
    assert_eq!(
        result.sequence, template,
        "zero-error reads must recover the exact template"
    );
}

#[test]
fn hifi_error_high_identity() {
    // 0.1% substitution, 0.05% ins/del -- typical PacBio HiFi profile
    let template = repeat(b"GCTAGCT", 10); // 70 bp
    let reads = Sim::new(0.001, 0.0005, 0.0005, 2).simulate(&template, 20);
    let result = consensus(&refs(&reads), 0, &PoaConfig::default()).unwrap();
    let id = identity(&template, &result.sequence);
    assert!(
        id >= 0.99,
        "HiFi error model: expected ≥99% identity, got {:.1}%",
        id * 100.0
    );
}

#[test]
fn ont_error_sufficient_depth() {
    // 5% substitution, 2% ins/del -- typical ONT profile
    let template = repeat(b"GCTAGCT", 10); // 70 bp
    let reads = Sim::new(0.05, 0.02, 0.02, 3).simulate(&template, 30);
    let result = consensus(&refs(&reads), 0, &PoaConfig::default()).unwrap();
    let id = identity(&template, &result.sequence);
    assert!(
        id >= 0.90,
        "ONT error model (30x): expected ≥90% identity, got {:.1}%",
        id * 100.0
    );
}

#[test]
fn str_majority_repeat_count_wins() {
    // 10 reads with 6 CAT repeats (majority), 4 reads with 7 repeats (minority)
    let majority = repeat(b"CAT", 6); // 18 bp
    let minority = repeat(b"CAT", 7); // 21 bp
    let mut reads: Vec<Vec<u8>> = vec![majority.clone(); 10];
    reads.extend(vec![minority; 4]);
    let result = consensus(&refs(&reads), 0, &PoaConfig::default()).unwrap();
    assert_eq!(
        result.sequence.len(),
        majority.len(),
        "majority repeat count should win: expected {}bp, got {}bp",
        majority.len(),
        result.sequence.len()
    );
}

#[test]
fn str_with_errors_majority_length_wins() {
    // Majority: 8 GAA repeats; minority: 10 repeats. Both with ONT errors.
    let template_a = repeat(b"GAA", 8); // 24 bp
    let template_b = repeat(b"GAA", 10); // 30 bp
    let mut reads = Sim::new(0.05, 0.02, 0.02, 4).simulate(&template_a, 15);
    reads.extend(Sim::new(0.05, 0.02, 0.02, 5).simulate(&template_b, 5));
    let result = consensus(&refs(&reads), 0, &PoaConfig::default()).unwrap();
    let id_a = identity(&template_a, &result.sequence);
    let id_b = identity(&template_b, &result.sequence);
    assert!(
        id_a > id_b,
        "majority allele (8×GAA) should be closer to consensus than minority (10×GAA): \
         id_a={:.1}% id_b={:.1}%",
        id_a * 100.0,
        id_b * 100.0
    );
}

#[test]
fn two_allele_snp_recovery() {
    // Two 22bp alleles differing at position 12 (A vs G).
    // SNP-based divergence creates a proper bubble fork in the graph;
    // length-based divergence creates a tail extension, not a detectable bubble.
    let allele_a = b"ACGTACGTACGTACGTACGTAC"; // pos 12 = A
    let allele_b = b"ACGTACGTACGTGCGTACGTAC"; // pos 12 = G
    let sim_a = Sim::new(0.01, 0.005, 0.005, 6);
    let sim_b = Sim::new(0.01, 0.005, 0.005, 7);
    let mut reads = sim_a.simulate(allele_a, 10);
    reads.extend(sim_b.simulate(allele_b, 10));
    let alleles = consensus_multi(&refs(&reads), 0, &PoaConfig::default()).unwrap();
    assert_eq!(
        alleles.len(),
        2,
        "expected 2 alleles, got {}",
        alleles.len()
    );
    let best_a = alleles
        .iter()
        .map(|a| identity(allele_a, &a.sequence))
        .fold(0.0_f64, f64::max);
    let best_b = alleles
        .iter()
        .map(|a| identity(allele_b, &a.sequence))
        .fold(0.0_f64, f64::max);
    assert!(
        best_a >= 0.90,
        "allele A recovery: expected ≥90%, got {:.1}%",
        best_a * 100.0
    );
    assert!(
        best_b >= 0.90,
        "allele B recovery: expected ≥90%, got {:.1}%",
        best_b * 100.0
    );
}

#[test]
fn frda_gaa_simulated_ont() {
    // FRDA locus: flanking sequence + GAA repeat + flanking sequence, ONT error profile.
    let template: Vec<u8> = [
        b"ACGTACGTACGT".as_slice(),
        &repeat(b"GAA", 20),
        b"TGCATGCATGCA",
    ]
    .concat();
    let reads = Sim::new(0.05, 0.02, 0.02, 8).simulate(&template, 20);
    let result = consensus(&refs(&reads), 0, &PoaConfig::default()).unwrap();
    let id = identity(&template, &result.sequence);
    assert!(
        id >= 0.90,
        "FRDA GAA (ONT, 20x): expected ≥90% identity, got {:.1}%",
        id * 100.0
    );
}

// ── Graph-size-aware adaptive band ───────────────────────────────────────────

#[test]
fn adaptive_band_graph_size_aware_closes_accuracy_gap() {
    // With a fixed narrow band, accuracy degrades as reads are added because
    // Insert ops from each read expand the graph, shifting topological ranks
    // away from the read-position diagonal.  The adaptive band formula
    // w = b + f×max(read_len, graph_nodes) compensates by widening the band
    // as the graph grows, matching unbanded accuracy at much lower per-read cost.
    let template = repeat(b"GCTAGCT", 10); // 70 bp
    let reads = Sim::new(0.05, 0.03, 0.03, 11).simulate(&template, 30);
    let refs = refs(&reads);

    let narrow = PoaConfig {
        band_width: 20,
        ..PoaConfig::default()
    };
    let adaptive = PoaConfig {
        adaptive_band: true,
        adaptive_band_b: 10,
        adaptive_band_f: 0.5, // w = 10 + 0.5×max(read_len, graph_nodes)
        ..PoaConfig::default()
    };
    let unbanded = PoaConfig::default();

    let id_narrow = consensus(&refs, 0, &narrow)
        .map(|r| identity(&template, &r.sequence))
        .unwrap_or(0.0);
    let id_adaptive = consensus(&refs, 0, &adaptive)
        .map(|r| identity(&template, &r.sequence))
        .unwrap_or(0.0);
    let id_unbanded = consensus(&refs, 0, &unbanded)
        .map(|r| identity(&template, &r.sequence))
        .unwrap_or(0.0);

    println!("\nGraph-size-aware adaptive band (30 reads, 5% sub, 3% ins/del, 70bp)");
    println!("  narrow (w=20):   {:.1}%", id_narrow * 100.0);
    println!("  adaptive:        {:.1}%", id_adaptive * 100.0);
    println!("  unbanded:        {:.1}%", id_unbanded * 100.0);

    assert!(
        id_adaptive >= id_narrow,
        "adaptive band should match or beat narrow band: adaptive={:.1}% narrow={:.1}%",
        id_adaptive * 100.0,
        id_narrow * 100.0
    );
    assert!(
        (id_adaptive - id_unbanded).abs() <= 0.05,
        "adaptive band should be within 5% of unbanded: adaptive={:.1}% unbanded={:.1}%",
        id_adaptive * 100.0,
        id_unbanded * 100.0
    );
}

// ── Informational (print tables, assert at key thresholds) ────────────────────

#[test]
fn depth_accuracy_curve() {
    // Unbanded so depth genuinely improves accuracy without band-expansion effects.
    let template = repeat(b"GCTAGCT", 10); // 70 bp, ONT error model
    let sim = Sim::new(0.05, 0.02, 0.02, 9);
    println!("\nDepth accuracy curve  (5% sub, 2% ins/del, unbanded, 70bp template)");
    println!("{:>6}  {:>10}  {:>8}", "depth", "identity%", "len");
    for &depth in &[3usize, 5, 10, 20, 50] {
        let reads = sim.simulate(&template, depth);
        match consensus(&refs(&reads), 0, &PoaConfig::default()) {
            Ok(r) => {
                let id = identity(&template, &r.sequence);
                println!(
                    "{:>6}  {:>9.1}%  {:>8}",
                    depth,
                    id * 100.0,
                    r.sequence.len()
                );
                if depth >= 10 {
                    assert!(
                        id >= 0.90,
                        "at depth {depth}: expected ≥90% identity, got {:.1}%",
                        id * 100.0
                    );
                }
            }
            Err(e) => println!("{:>6}  {:>10}  {:>8}", depth, format!("Err({e})"), "-"),
        }
    }
}

#[test]
fn band_width_accuracy_tradeoff() {
    // Shows how accuracy degrades as band width narrows.
    // Banded DP is stressed more than per-read drift alone suggests because each
    // read's insert operations add new graph nodes, shifting topological ranks.
    // With 20 reads at 3% ins/del on a 70bp template, the graph accumulates
    // ~42 extra insert nodes, requiring an effective band of ~40+ for full accuracy.
    let template = repeat(b"GCTAGCT", 10); // 70 bp
    let reads = Sim::new(0.05, 0.03, 0.03, 10).simulate(&template, 20);
    println!("\nBand width tradeoff  (5% sub, 3% ins/del, 20 reads, 70bp template)");
    println!("{:>12}  {:>10}  {:>8}", "band_width", "identity%", "len");
    for &bw in &[0usize, 5, 10, 20, 50, 100] {
        let config = PoaConfig {
            band_width: bw,
            warn_on_long_unbanded: false,
            ..PoaConfig::default()
        };
        match consensus(&refs(&reads), 0, &config) {
            Ok(r) => {
                let id = identity(&template, &r.sequence);
                let label = if bw == 0 {
                    "unbanded".to_string()
                } else {
                    bw.to_string()
                };
                println!(
                    "{:>12}  {:>9.1}%  {:>8}",
                    label,
                    id * 100.0,
                    r.sequence.len()
                );
            }
            Err(e) => {
                let label = if bw == 0 {
                    "unbanded".to_string()
                } else {
                    bw.to_string()
                };
                println!("{:>12}  {:>10}  {:>8}", label, format!("Err({e})"), "-");
            }
        }
    }
}

// ── Plots (require --features plot) ──────────────────────────────────────────

#[cfg(feature = "plot")]
mod accuracy_plots {
    use super::*;
    use kuva::prelude::*;
    use std::fs;
    use std::path::PathBuf;

    fn save(name: &str, svg: String) {
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test_plots");
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        fs::write(&path, svg).unwrap();
        println!("wrote {}", path.display());
    }

    /// Identity% vs depth for HiFi and ONT error models, averaged over N seeds.
    ///
    /// Single-seed runs show non-monotonic accuracy because min_cov = n/2+1 scales
    /// with depth, and specific unlucky seeds produce reads that interact with the
    /// boundary trim at particular depths.  Averaging over multiple seeds reveals
    /// the true monotonically-increasing expected trend.
    #[test]
    fn plot_depth_accuracy_curve() {
        const N_SEEDS: u64 = 20;
        let template = repeat(b"GCTAGCT", 10); // 70 bp
        let depths = [3usize, 5, 10, 20, 30, 50];

        let mean_identity = |sub: f64, ins: f64, del: f64, depth: usize| -> f64 {
            let sum: f64 = (0..N_SEEDS)
                .filter_map(|seed| {
                    let reads =
                        Sim::new(sub, ins, del, (seed + 1) * 100).simulate(&template, depth);
                    consensus(&refs(&reads), 0, &PoaConfig::default())
                        .ok()
                        .map(|r| identity(&template, &r.sequence) * 100.0)
                })
                .sum();
            sum / N_SEEDS as f64
        };

        let hifi_pts: Vec<(f64, f64)> = depths
            .iter()
            .map(|&d| (d as f64, mean_identity(0.001, 0.0005, 0.0005, d)))
            .collect();

        let ont_pts: Vec<(f64, f64)> = depths
            .iter()
            .map(|&d| (d as f64, mean_identity(0.05, 0.02, 0.02, d)))
            .collect();

        let hifi_line: Plot = LinePlot::new()
            .with_data(hifi_pts)
            .with_color("steelblue")
            .with_legend(&format!("HiFi (0.1% sub, mean of {N_SEEDS} seeds)"))
            .into();
        let ont_line: Plot = LinePlot::new()
            .with_data(ont_pts)
            .with_color("darkorange")
            .with_legend(&format!("ONT (5% sub, 2% indel, mean of {N_SEEDS} seeds)"))
            .into();

        let plots = vec![hifi_line, ont_line];
        let layout = Layout::auto_from_plots(&plots)
            .with_title("Consensus accuracy vs depth (70 bp template, unbanded)")
            .with_x_label("Depth (reads)")
            .with_y_label("Identity (%)");

        save("accuracy_depth_curve.svg", render_to_svg(plots, layout));
    }

    /// Identity% vs band width, showing graph-expansion degradation, averaged over N seeds.
    /// Each series is a different read depth; dashed lines mark the unbanded reference.
    #[test]
    fn plot_band_width_accuracy_curve() {
        const N_SEEDS: u64 = 20;
        let template = repeat(b"GCTAGCT", 10); // 70 bp, ONT 5% sub 3% indel
        let band_widths = [5usize, 10, 20, 50, 100];
        let depths = [(10usize, "steelblue"), (20, "darkorange"), (30, "seagreen")];

        let mean_id = |depth: usize, bw: usize| -> Option<f64> {
            let config = PoaConfig {
                band_width: bw,
                warn_on_long_unbanded: false,
                ..PoaConfig::default()
            };
            let vals: Vec<f64> = (0..N_SEEDS)
                .filter_map(|seed| {
                    let reads =
                        Sim::new(0.05, 0.03, 0.03, (seed + 1) * 100).simulate(&template, depth);
                    consensus(&refs(&reads), 0, &config)
                        .ok()
                        .map(|r| identity(&template, &r.sequence) * 100.0)
                })
                .collect();
            if vals.is_empty() {
                None
            } else {
                Some(vals.iter().sum::<f64>() / vals.len() as f64)
            }
        };

        let mut plots: Vec<Plot> = Vec::new();

        for (depth, color) in depths {
            let banded_pts: Vec<(f64, f64)> = band_widths
                .iter()
                .filter_map(|&bw| mean_id(depth, bw).map(|id| (bw as f64, id)))
                .collect();

            let unbanded_id = mean_id(depth, 0).unwrap_or(0.0);

            let label = format!("{depth}x reads (mean of {N_SEEDS})");
            plots.push(
                LinePlot::new()
                    .with_data(banded_pts)
                    .with_color(color)
                    .with_legend(&label)
                    .into(),
            );
            let unbanded_line: Plot = LinePlot::new()
                .with_data(vec![(5.0, unbanded_id), (100.0, unbanded_id)])
                .with_color(color)
                .with_stroke_width(1.0)
                .with_dashed()
                .into();
            plots.push(unbanded_line);
        }

        let layout = Layout::auto_from_plots(&plots)
            .with_title("Consensus accuracy vs band width (70 bp, ONT 5% sub 3% indel)\n(dashed = unbanded reference)")
            .with_x_label("Band width (cells)")
            .with_y_label("Identity (%)");

        save(
            "accuracy_band_width_curve.svg",
            render_to_svg(plots, layout),
        );
    }

    /// Identity% vs substitution rate at fixed depth (20 reads, unbanded), averaged over N seeds.
    #[test]
    fn plot_error_rate_sweep() {
        const N_SEEDS: u64 = 20;
        let template = repeat(b"GCTAGCT", 10); // 70 bp
        let sub_rates = [0.0f64, 0.01, 0.02, 0.05, 0.08, 0.10];

        let mean_id = |sub: f64, ins: f64, del: f64| -> f64 {
            let sum: f64 = (0..N_SEEDS)
                .filter_map(|seed| {
                    let reads = Sim::new(sub, ins, del, (seed + 1) * 100).simulate(&template, 20);
                    consensus(&refs(&reads), 0, &PoaConfig::default())
                        .ok()
                        .map(|r| identity(&template, &r.sequence) * 100.0)
                })
                .sum();
            sum / N_SEEDS as f64
        };

        let low_indel_pts: Vec<(f64, f64)> = sub_rates
            .iter()
            .map(|&s| (s * 100.0, mean_id(s, 0.005, 0.005)))
            .collect();

        let high_indel_pts: Vec<(f64, f64)> = sub_rates
            .iter()
            .map(|&s| (s * 100.0, mean_id(s, s / 2.0, s / 2.0)))
            .collect();

        let low_line: Plot = LinePlot::new()
            .with_data(low_indel_pts)
            .with_color("steelblue")
            .with_legend(&format!("low indel (0.5% ins/del, mean of {N_SEEDS})"))
            .into();
        let high_line: Plot = LinePlot::new()
            .with_data(high_indel_pts)
            .with_color("darkorange")
            .with_legend(&format!("high indel (ins/del = sub/2, mean of {N_SEEDS})"))
            .into();

        let plots = vec![low_line, high_line];
        let layout = Layout::auto_from_plots(&plots)
            .with_title("Consensus accuracy vs error rate (70 bp, 20 reads, unbanded)")
            .with_x_label("Substitution rate (%)")
            .with_y_label("Identity (%)");

        save(
            "accuracy_error_rate_sweep.svg",
            render_to_svg(plots, layout),
        );
    }
}


// ── SCA4_ZFHX3 Hap2 diagnostic ────────────────────────────────────────────────

/// Reproduce the SCA4_ZFHX3 Hap2 consensus failure.
///
/// 5 of 8 reads are identical (63 bp). The correct consensus is
/// GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC.
/// Bladerunner produced CCGCCGCGCCGCCGCCGCCGCCCCGCCGCCGCCGCACTGCCACCGCCGCCG (51 bp)
/// which is wrong -- phase-shifted and truncated.
///
/// This test sweeps every read as seed to identify which seed causes failure
/// and which gives the correct result.
#[test]
fn diag_sca4_zfhx3_hap2_seed_sweep() {
    let reads: Vec<&[u8]> = vec![
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // read 8  (63 bp)
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // read 9  (63 bp)
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // read 10 (63 bp)
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC",     // read 11 (60 bp)
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCC",        // read 12 (57 bp)
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // read 13 (63 bp)
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // read 14 (63 bp)
        b"GCCGCCGCGCCGCCGCCGCCGCCCCGCGCCGCCGCCGCACTGCCACCGCCGCCGCCGCC",     // read 15 (59 bp, outlier)
    ];
    let expected = b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC";

    let cfg = poa_consensus::PoaConfig {
        adaptive_band: true,
        band_width: 50,
        alignment_mode: poa_consensus::AlignmentMode::SemiGlobal,
        ..poa_consensus::PoaConfig::default()
    };

    println!("\n=== SCA4_ZFHX3 Hap2 seed sweep ===");
    println!("Expected ({} bp): {}", expected.len(), std::str::from_utf8(expected).unwrap());
    println!();

    for seed_idx in 0..reads.len() {
        let result = poa_consensus::consensus(&reads, seed_idx, &cfg)
            .expect("consensus failed");
        let cons = &result.sequence;
        let dist = levenshtein(cons, expected);
        let matched = cons == expected;
        println!(
            "seed={} ({}bp): got {}bp dist={} {}  -- {}",
            seed_idx,
            reads[seed_idx].len(),
            cons.len(),
            dist,
            if matched { "CORRECT" } else { "WRONG" },
            std::str::from_utf8(cons).unwrap(),
        );
    }
    println!();

    // At least one seed (any of reads 0-4 which are the 63 bp majority) must give correct result.
    let any_correct = (0..reads.len()).any(|seed_idx| {
        poa_consensus::consensus(&reads, seed_idx, &cfg)
            .map(|r| r.sequence == expected.to_vec())
            .unwrap_or(false)
    });
    assert!(any_correct, "no seed produced the correct consensus");
}


/// Try to reproduce the bladerunner wrong output for SCA4_ZFHX3 Hap2
/// by sweeping configs (alignment mode, band settings).
#[test]
fn diag_sca4_zfhx3_hap2_config_sweep() {
    let reads: Vec<&[u8]> = vec![
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC",
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC",
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC",
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC",
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCC",
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC",
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC",
        b"GCCGCCGCGCCGCCGCCGCCGCCCCGCGCCGCCGCCGCACTGCCACCGCCGCCGCCGCC",
    ];
    let expected = b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC";
    let bladerunner_got = b"CCGCCGCGCCGCCGCCGCCGCCCCGCCGCCGCCGCACTGCCACCGCCGCCG";

    use poa_consensus::AlignmentMode;

    let configs: &[(&str, poa_consensus::PoaConfig)] = &[
        ("adaptive+semi-global+bw50 (recommended)", poa_consensus::PoaConfig {
            adaptive_band: true, band_width: 50,
            alignment_mode: AlignmentMode::SemiGlobal,
            ..poa_consensus::PoaConfig::default()
        }),
        ("adaptive+global+bw50", poa_consensus::PoaConfig {
            adaptive_band: true, band_width: 50,
            alignment_mode: AlignmentMode::Global,
            ..poa_consensus::PoaConfig::default()
        }),
        ("fixed-bw50+semi-global", poa_consensus::PoaConfig {
            adaptive_band: false, band_width: 50,
            alignment_mode: AlignmentMode::SemiGlobal,
            ..poa_consensus::PoaConfig::default()
        }),
        ("fixed-bw50+global", poa_consensus::PoaConfig {
            adaptive_band: false, band_width: 50,
            alignment_mode: AlignmentMode::Global,
            ..poa_consensus::PoaConfig::default()
        }),
        ("adaptive+semi-global+bw10 (default adaptive_band_b)", poa_consensus::PoaConfig {
            adaptive_band: true, band_width: 10,
            alignment_mode: AlignmentMode::SemiGlobal,
            ..poa_consensus::PoaConfig::default()
        }),
        ("adaptive+global+bw10", poa_consensus::PoaConfig {
            adaptive_band: true, band_width: 10,
            alignment_mode: AlignmentMode::Global,
            ..poa_consensus::PoaConfig::default()
        }),
    ];

    println!("\n=== SCA4_ZFHX3 Hap2 config sweep (seed=7, the outlier) ===");
    println!("Expected ({} bp):       {}", expected.len(), std::str::from_utf8(expected).unwrap());
    println!("Bladerunner got ({} bp): {}", bladerunner_got.len(), std::str::from_utf8(bladerunner_got).unwrap());
    println!();

    for (name, cfg) in configs {
        // Try both the outlier (seed=7) and a 63bp read (seed=0) as seed
        for seed_idx in [7usize, 0] {
            let result = poa_consensus::consensus(&reads, seed_idx, cfg)
                .expect("consensus failed");
            let cons = &result.sequence;
            let dist_expected = levenshtein(cons, expected);
            let dist_bladerunner = levenshtein(cons, bladerunner_got);
            println!(
                "cfg={:50} seed={}: got {}bp  dist_expected={} dist_bladerunner={}  -- {}",
                name, seed_idx, cons.len(), dist_expected, dist_bladerunner,
                std::str::from_utf8(cons).unwrap(),
            );
        }
    }
}


/// Test POA with padded reads (BED ± 20 bp flanks) as bladerunner sends them.
///
/// Left flank (hg38):  CACGCCAGGCAGTGGTACGA
/// Right flank (hg38): GTGGGGACGTGAAGCACCAT
///
/// The Hap2 consensus from bladerunner was 51 bp in CCG phase. This test
/// checks whether POA on the padded reads also produces a wrong result.
#[test]
fn diag_sca4_zfhx3_hap2_padded() {
    let lf = b"CACGCCAGGCAGTGGTACGA";
    let rf = b"GTGGGGACGTGAAGCACCAT";

    // unpadded Hap2 reads (from bladerunner_str_output.tsv)
    let repeats: &[&[u8]] = &[
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // read 8  (63)
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // read 9  (63)
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // read 10 (63)
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC",    // read 11 (60)
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCC",       // read 12 (57)
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // read 13 (63)
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // read 14 (63)
        b"GCCGCCGCGCCGCCGCCGCCGCCCCGCGCCGCCGCCGCACTGCCACCGCCGCCGCCGCC",    // read 15 (59)
    ];

    // build padded sequences
    let padded: Vec<Vec<u8>> = repeats.iter().map(|r| {
        let mut v = Vec::with_capacity(lf.len() + r.len() + rf.len());
        v.extend_from_slice(lf);
        v.extend_from_slice(r);
        v.extend_from_slice(rf);
        v
    }).collect();

    let expected_repeat = b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC";
    let bladerunner_got_repeat = b"CCGCCGCGCCGCCGCCGCCGCCCCGCCGCCGCCGCACTGCCACCGCCGCCG";

    let cfg = poa_consensus::PoaConfig::default();

    println!("\n=== SCA4_ZFHX3 Hap2 with padded reads (BED ± 20 bp) ===");
    println!("Padded lengths: {:?}", padded.iter().map(|r| r.len()).collect::<Vec<_>>());

    for seed_idx in 0..padded.len() {
        let slices: Vec<&[u8]> = padded.iter().map(|r| r.as_slice()).collect();
        let result = poa_consensus::consensus(&slices, seed_idx, &cfg).expect("consensus failed");
        let cons = &result.sequence;

        // strip flanks from result for comparison (if present)
        let repeat_part: &[u8] = if cons.len() > lf.len() + rf.len()
            && cons.starts_with(lf) && cons.ends_with(rf)
        {
            &cons[lf.len()..cons.len() - rf.len()]
        } else {
            cons.as_slice()
        };

        let dist_exp = levenshtein(repeat_part, expected_repeat);
        let dist_br  = levenshtein(repeat_part, bladerunner_got_repeat);
        println!(
            "seed={} ({}bp total): cons={}bp  repeat_part={}bp  dist_expected={}  dist_br={}  -- {}",
            seed_idx, padded[seed_idx].len(), cons.len(), repeat_part.len(),
            dist_exp, dist_br,
            std::str::from_utf8(cons).unwrap(),
        );
    }
}


/// Step-by-step addition for padded reads with seed=0 vs seed=3 to find which
/// read corrupts the graph.
#[test]
fn diag_sca4_zfhx3_hap2_padded_stepwise() {
    let lf = b"CACGCCAGGCAGTGGTACGA";
    let rf = b"GTGGGGACGTGAAGCACCAT";
    let repeats: &[&[u8]] = &[
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // 0 r8  63
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // 1 r9  63
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // 2 r10 63
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC",    // 3 r11 60
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCC",       // 4 r12 57
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // 5 r13 63
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // 6 r14 63
        b"GCCGCCGCGCCGCCGCCGCCGCCCCGCGCCGCCGCCGCACTGCCACCGCCGCCGCCGCC",    // 7 r15 59 (outlier)
    ];
    let padded: Vec<Vec<u8>> = repeats.iter().map(|r| {
        let mut v = Vec::with_capacity(lf.len() + r.len() + rf.len());
        v.extend_from_slice(lf); v.extend_from_slice(r); v.extend_from_slice(rf);
        v
    }).collect();
    let expected = b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC";
    let cfg = poa_consensus::PoaConfig::default();

    println!("\n=== Stepwise padded POA for seed=0 and seed=3 ===");

    for seed_idx in [0usize, 3, 4, 5] {
        println!("\n--- seed={} ({}bp padded) ---", seed_idx, padded[seed_idx].len());
        let mut graph = poa_consensus::PoaGraph::new(&padded[seed_idx], cfg.clone()).unwrap();
        for (i, r) in padded.iter().enumerate() {
            if i == seed_idx { continue; }
            graph.add_read(r).unwrap();
            let c = graph.consensus().unwrap();
            let seq = &c.sequence;
            // strip flanks if present
            let rep: &[u8] = if seq.len() > lf.len() + rf.len() && seq.starts_with(lf) && seq.ends_with(rf) {
                &seq[lf.len()..seq.len()-rf.len()]
            } else { seq.as_slice() };
            let dist = levenshtein(rep, expected);
            println!(
                "  after adding r{} ({}bp): cons={}bp repeat={}bp dist={} -- {}",
                i, padded[i].len(), seq.len(), rep.len(), dist,
                std::str::from_utf8(seq).unwrap()
            );
        }
    }
}


/// For seed=0, show the alignment ops that each read produces,
/// to identify which read creates wrong INSERT nodes.
#[test]
fn diag_sca4_zfhx3_hap2_alignment_ops() {
    use poa_consensus::AlignOp;
    let lf = b"CACGCCAGGCAGTGGTACGA";
    let rf = b"GTGGGGACGTGAAGCACCAT";
    let repeats: &[&[u8]] = &[
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // 0 r8  63
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // 1 r9  63
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // 2 r10 63
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC",    // 3 r11 60
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCC",       // 4 r12 57
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // 5 r13 63
        b"GCCGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCACTGCCACCGCCGCCGCCGCC", // 6 r14 63
        b"GCCGCCGCGCCGCCGCCGCCGCCCCGCGCCGCCGCCGCACTGCCACCGCCGCCGCCGCC",    // 7 r15 59 (outlier)
    ];
    let padded: Vec<Vec<u8>> = repeats.iter().map(|r| {
        let mut v = Vec::with_capacity(lf.len() + r.len() + rf.len());
        v.extend_from_slice(lf); v.extend_from_slice(r); v.extend_from_slice(rf);
        v
    }).collect();
    let cfg = poa_consensus::PoaConfig::default();

    println!("\n=== Alignment ops for seed=0, one read at a time ===");
    let seed_idx = 0;
    let mut graph = poa_consensus::PoaGraph::new(&padded[seed_idx], cfg.clone()).unwrap();

    for (i, r) in padded.iter().enumerate() {
        if i == seed_idx { continue; }

        let (ops, _, _) = graph.align_read_ops(r).unwrap();

        // Summarize ops: count M/I/D and flag inserts
        let matches = ops.iter().filter(|o| matches!(o, AlignOp::Match(_))).count();
        let inserts: Vec<u8> = ops.iter().filter_map(|o| if let AlignOp::Insert(b) = o { Some(*b) } else { None }).collect();
        let deletes = ops.iter().filter(|o| matches!(o, AlignOp::Delete(_))).count();

        let insert_str = if inserts.is_empty() {
            String::from("none")
        } else {
            String::from_utf8_lossy(&inserts).to_string()
        };

        println!(
            "  r{} ({}bp): M={} I={} D={}  inserts={}  graph_nodes={}",
            i, padded[i].len(), matches, inserts.len(), deletes, insert_str,
            graph.node_count()
        );

        let (add_ops, anchor_count, anchor_list) = graph.add_read_debug(r).unwrap();
        let add_matches = add_ops.iter().filter(|o| matches!(o, AlignOp::Match(_))).count();
        let add_inserts: Vec<u8> = add_ops.iter().filter_map(|o| if let AlignOp::Insert(b) = o { Some(b) } else { None }).cloned().collect();
        let add_deletes = add_ops.iter().filter(|o| matches!(o, AlignOp::Delete(_))).count();
        let add_insert_str = if add_inserts.is_empty() { String::from("none") } else { String::from_utf8_lossy(&add_inserts).to_string() };
        println!(
            "    add_read: M={} I={} D={}  inserts={}  anchors={}  graph_nodes={}",
            add_matches, add_inserts.len(), add_deletes, add_insert_str, anchor_count,
            graph.node_count()
        );
        if i == 3 {  // r3 is the first problem case
            println!("    anchor_list (read_pos, topo_rank): {:?}", anchor_list);
        }
        let c = graph.consensus().unwrap();
        println!("    -> cons={}bp", c.sequence.len());
    }
}

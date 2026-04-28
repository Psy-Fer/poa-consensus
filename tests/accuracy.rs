//! Accuracy tests for POA consensus under simulated error models.
//!
//! Tests measure how well the library recovers a known ground-truth sequence
//! across error rates, depths, band widths, and repeat structures.
//!
//! Run with:
//!   cargo test --test accuracy -- --nocapture

use poa_consensus::{consensus, consensus_multi, PoaConfig};

// ── Simulator ─────────────────────────────────────────────────────────────────

struct Sim {
    sub_rate: f64,
    ins_rate: f64,
    del_rate: f64,
    seed: u64,
}

impl Sim {
    fn new(sub_rate: f64, ins_rate: f64, del_rate: f64, seed: u64) -> Self {
        Self { sub_rate, ins_rate, del_rate, seed }
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
        _    => [b'A', b'C', b'G'],
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
    let template_a = repeat(b"GAA", 8);  // 24 bp
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
    assert_eq!(alleles.len(), 2, "expected 2 alleles, got {}", alleles.len());
    let best_a = alleles.iter().map(|a| identity(allele_a, &a.sequence)).fold(0.0_f64, f64::max);
    let best_b = alleles.iter().map(|a| identity(allele_b, &a.sequence)).fold(0.0_f64, f64::max);
    assert!(best_a >= 0.90, "allele A recovery: expected ≥90%, got {:.1}%", best_a * 100.0);
    assert!(best_b >= 0.90, "allele B recovery: expected ≥90%, got {:.1}%", best_b * 100.0);
}

#[test]
fn frda_gaa_simulated_ont() {
    // FRDA locus: flanking sequence + GAA repeat + flanking sequence, ONT error profile.
    let template: Vec<u8> =
        [b"ACGTACGTACGT".as_slice(), &repeat(b"GAA", 20), b"TGCATGCATGCA"].concat();
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

    let narrow = PoaConfig { band_width: 20, ..PoaConfig::default() };
    let adaptive = PoaConfig {
        adaptive_band: true,
        adaptive_band_b: 10,
        adaptive_band_f: 0.5, // w = 10 + 0.5×max(read_len, graph_nodes)
        ..PoaConfig::default()
    };
    let unbanded = PoaConfig::default();

    let id_narrow   = consensus(&refs, 0, &narrow).map(|r| identity(&template, &r.sequence)).unwrap_or(0.0);
    let id_adaptive = consensus(&refs, 0, &adaptive).map(|r| identity(&template, &r.sequence)).unwrap_or(0.0);
    let id_unbanded = consensus(&refs, 0, &unbanded).map(|r| identity(&template, &r.sequence)).unwrap_or(0.0);

    println!("\nGraph-size-aware adaptive band (30 reads, 5% sub, 3% ins/del, 70bp)");
    println!("  narrow (w=20):   {:.1}%", id_narrow   * 100.0);
    println!("  adaptive:        {:.1}%", id_adaptive * 100.0);
    println!("  unbanded:        {:.1}%", id_unbanded * 100.0);

    assert!(
        id_adaptive >= id_narrow,
        "adaptive band should match or beat narrow band: adaptive={:.1}% narrow={:.1}%",
        id_adaptive * 100.0, id_narrow * 100.0
    );
    assert!(
        (id_adaptive - id_unbanded).abs() <= 0.05,
        "adaptive band should be within 5% of unbanded: adaptive={:.1}% unbanded={:.1}%",
        id_adaptive * 100.0, id_unbanded * 100.0
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
                println!("{:>6}  {:>9.1}%  {:>8}", depth, id * 100.0, r.sequence.len());
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
                let label = if bw == 0 { "unbanded".to_string() } else { bw.to_string() };
                println!("{:>12}  {:>9.1}%  {:>8}", label, id * 100.0, r.sequence.len());
            }
            Err(e) => {
                let label = if bw == 0 { "unbanded".to_string() } else { bw.to_string() };
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
                    let reads = Sim::new(sub, ins, del, (seed + 1) * 100).simulate(&template, depth);
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
                    let reads = Sim::new(0.05, 0.03, 0.03, (seed + 1) * 100).simulate(&template, depth);
                    consensus(&refs(&reads), 0, &config)
                        .ok()
                        .map(|r| identity(&template, &r.sequence) * 100.0)
                })
                .collect();
            if vals.is_empty() { None } else { Some(vals.iter().sum::<f64>() / vals.len() as f64) }
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
                LinePlot::new().with_data(banded_pts).with_color(color).with_legend(&label).into(),
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

        save("accuracy_band_width_curve.svg", render_to_svg(plots, layout));
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

        save("accuracy_error_rate_sweep.svg", render_to_svg(plots, layout));
    }
}

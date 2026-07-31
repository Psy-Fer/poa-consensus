//! Comprehensive flanked-vs-raw characterization sweep for the flank-anchoring
//! API (`consensus_flanked` / `consensus_multi_flanked`). Reports where anchoring
//! WINS, is NEUTRAL, or REGRESSES vs the raw engine across motif × units × depth
//! × error × partial-fraction × flank-size, aggregated by partial fraction (the
//! axis where the win is expected). Ignored by default (report, not assertion).
//!
//! Run: `cargo test --release --test flank_anchor_sweep -- --ignored --nocapture`

mod synth;

use poa_consensus::{
    PoaConfig, consensus, consensus_flanked, consensus_multi, consensus_multi_flanked,
};
use synth::{AlleleSpec, Spec};

const MOTIFS: &[(&str, &[u8])] = &[
    ("CAG", b"CAG"),
    ("GAA", b"GAA"),
    ("AAAAG", b"AAAAG"),
    ("CATCAT", b"CATCAT"),
];
const SEEDS: &[u64] = &[1, 2, 3];

fn cfg() -> PoaConfig {
    PoaConfig {
        min_reads: 3,
        min_allele_freq: 0.2,
        phasing_bubble_min_span: 10,
        ..PoaConfig::default()
    }
}

fn unique_flank(size: usize, seed: u64) -> Vec<u8> {
    synth::random_seq(size, seed)
}

fn spec(
    left: &[u8],
    right: &[u8],
    alleles: Vec<AlleleSpec>,
    err: f64,
    partial: f64,
    unit_len: usize,
    seed: u64,
) -> Spec {
    Spec {
        left: left.to_vec(),
        right: right.to_vec(),
        alleles,
        sub: err,
        ins: err,
        del: err,
        partial_frac: partial,
        unit_jitter: false,
        unit_len,
        seed,
    }
}

/// Bucketed win/neutral/loss counts keyed by partial fraction.
#[derive(Default, Clone)]
struct Tally {
    win: u32,
    neutral: u32,
    loss: u32,
    raw_sum: f64,
    anch_sum: f64,
    n: u32,
}

impl Tally {
    fn record(&mut self, raw: f64, anch: f64) {
        self.raw_sum += raw;
        self.anch_sum += anch;
        self.n += 1;
        // +1 tolerance so sub-unit noise isn't counted as a win/loss.
        if anch + 1.0 < raw {
            self.win += 1;
        } else if raw + 1.0 < anch {
            self.loss += 1;
        } else {
            self.neutral += 1;
        }
    }
    fn line(&self, label: &str) {
        println!(
            "  {label:<14} win {:>3}  neutral {:>3}  loss {:>3}   mean raw {:>5.1} -> anch {:>5.1}",
            self.win,
            self.neutral,
            self.loss,
            if self.n > 0 {
                self.raw_sum / self.n as f64
            } else {
                0.0
            },
            if self.n > 0 {
                self.anch_sum / self.n as f64
            } else {
                0.0
            },
        );
    }
}

#[test]
#[ignore = "characterization sweep; run with --ignored --nocapture"]
fn single_allele_flanked_vs_raw() {
    let units_set = [20usize, 40];
    let depths = [5usize, 10, 20, 40];
    let errs = [0.005f64, 0.01, 0.03];
    let partials = [0.0f64, 0.2, 0.4];
    let flanks = [20usize, 34];

    let mut by_partial: std::collections::BTreeMap<u32, Tally> = Default::default();
    for (_mn, m) in MOTIFS {
        for &units in &units_set {
            let mid = synth::repeat(m, units);
            for &depth in &depths {
                for &err in &errs {
                    for &partial in &partials {
                        for &flank in &flanks {
                            let left = unique_flank(flank, 7001);
                            let right = unique_flank(flank, 8002);
                            let truth_full: Vec<u8> = left
                                .iter()
                                .chain(mid.iter())
                                .chain(right.iter())
                                .copied()
                                .collect();
                            let (mut raw, mut anch) = (0.0, 0.0);
                            for &s in SEEDS {
                                let reads = synth::generate(&spec(
                                    &left,
                                    &right,
                                    vec![AlleleSpec {
                                        mid: mid.clone(),
                                        n_reads: depth,
                                    }],
                                    err,
                                    partial,
                                    m.len(),
                                    s,
                                ));
                                let seqs = synth::Read::seqs(&reads);
                                let r = consensus(&seqs, 0, &cfg())
                                    .map(|c| synth::edit_distance(&c.sequence, &truth_full))
                                    .unwrap_or(truth_full.len());
                                let a = consensus_flanked(&seqs, &left, &right, &cfg())
                                    .map(|c| synth::edit_distance(&c.sequence, &mid))
                                    .unwrap_or(mid.len());
                                raw += r as f64;
                                anch += a as f64;
                            }
                            let n = SEEDS.len() as f64;
                            by_partial
                                .entry((partial * 100.0) as u32)
                                .or_default()
                                .record(raw / n, anch / n);
                        }
                    }
                }
            }
        }
    }
    println!("\n===== SINGLE-ALLELE: flanked vs raw (edit-to-repeat-truth; lower=better) =====");
    for (p, t) in &by_partial {
        t.line(&format!("partial={}%", p));
    }
    println!("(win = anchoring beat raw by >1 edit; loss = raw beat anchoring by >1)\n");
}

#[test]
#[ignore = "characterization sweep; run with --ignored --nocapture"]
fn multi_length_flanked_vs_raw() {
    let gaps: &[(usize, usize)] = &[(20, 25), (20, 40), (40, 45)];
    let depths = [10usize, 20];
    let errs = [0.01f64, 0.03];
    let partials = [0.0f64, 0.3];
    let flanks = [20usize, 34];

    // "correct2" = exactly 2 alleles at ~{ua,ub} units (±1). Score 0 if correct
    // (so lower=better, comparable to Tally), 1 if not.
    fn score(units: &[i64], ua: i64, ub: i64) -> f64 {
        if units.len() == 2 {
            let (lo, hi) = (units[0].min(units[1]), units[0].max(units[1]));
            let (ta, tb) = (ua.min(ub), ua.max(ub));
            if (lo - ta).abs() <= 1 && (hi - tb).abs() <= 1 {
                return 0.0;
            }
        }
        1.0
    }

    let mut by_partial: std::collections::BTreeMap<u32, Tally> = Default::default();
    for (_mn, m) in &MOTIFS[..3] {
        for &(ua, ub) in gaps {
            for &depth in &depths {
                for &err in &errs {
                    for &partial in &partials {
                        for &flank in &flanks {
                            let left = unique_flank(flank, 7001);
                            let right = unique_flank(flank, 8002);
                            let (mut raw, mut anch) = (0.0, 0.0);
                            for &s in SEEDS {
                                let reads = synth::generate(&spec(
                                    &left,
                                    &right,
                                    vec![
                                        AlleleSpec {
                                            mid: synth::repeat(m, ua),
                                            n_reads: depth,
                                        },
                                        AlleleSpec {
                                            mid: synth::repeat(m, ub),
                                            n_reads: depth,
                                        },
                                    ],
                                    err,
                                    partial,
                                    m.len(),
                                    s,
                                ));
                                let seqs = synth::Read::seqs(&reads);
                                let ul = m.len();
                                let r = consensus_multi(&seqs, 0, &cfg())
                                    .map(|v| {
                                        let u: Vec<i64> = v
                                            .iter()
                                            .map(|c| {
                                                ((c.sequence.len().saturating_sub(2 * flank))
                                                    as f64
                                                    / ul as f64)
                                                    .round()
                                                    as i64
                                            })
                                            .collect();
                                        score(&u, ua as i64, ub as i64)
                                    })
                                    .unwrap_or(1.0);
                                let a = consensus_multi_flanked(&seqs, &left, &right, &cfg())
                                    .map(|v| {
                                        let u: Vec<i64> = v
                                            .iter()
                                            .map(|c| {
                                                (c.sequence.len() as f64 / ul as f64).round() as i64
                                            })
                                            .collect();
                                        score(&u, ua as i64, ub as i64)
                                    })
                                    .unwrap_or(1.0);
                                raw += r;
                                anch += a;
                            }
                            let n = SEEDS.len() as f64;
                            // scale to 0..10 so the >1 tolerance in Tally is meaningful per-config
                            by_partial
                                .entry((partial * 100.0) as u32)
                                .or_default()
                                .record(raw / n * 10.0, anch / n * 10.0);
                        }
                    }
                }
            }
        }
    }
    println!("\n===== MULTI-ALLELE LENGTH: flanked vs raw (miscall rate ×10; lower=better) =====");
    for (p, t) in &by_partial {
        t.line(&format!("partial={}%", p));
    }
    println!("(win = anchoring called 2 correct alleles where raw didn't)\n");
}

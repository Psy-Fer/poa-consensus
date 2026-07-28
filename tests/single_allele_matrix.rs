//! Single-allele consensus robustness matrix — the release correctness check for
//! the SHIPPING default path (`consensus()` → poa2). Sweeps a clean single allele
//! across motif × length × depth × error × partial-fraction × seed and scores the
//! consensus by edit-distance and unit-count vs ground truth. Same approach that
//! exposed the multi-allele fragility; here it certifies the single-allele path.
//!
//! Report: cargo test --release --features cli --test single_allele_matrix -- --ignored --nocapture
//! The plain `#[test]`s assert must-hold accuracy bounds.

mod synth;

use poa_consensus::{PoaConfig, consensus};
use synth::{AlleleSpec, Read, Spec};

fn cfg() -> PoaConfig {
    PoaConfig {
        min_reads: 3,
        ..PoaConfig::default()
    }
}

/// One clean single allele: flanks + `mid`, per-base error, optional partials.
/// `err` is the TOTAL per-base error, split ONT-style across substitution /
/// insertion / deletion channels (the generator applies the three rates
/// cumulatively, so they must sum to `err`, not each equal it).
fn spec(mid: Vec<u8>, unit_len: usize, depth: usize, err: f64, partial: f64, seed: u64) -> Spec {
    Spec {
        left: synth::LEFT.to_vec(),
        right: synth::RIGHT.to_vec(),
        alleles: vec![AlleleSpec {
            mid,
            n_reads: depth,
        }],
        sub: err * 0.34,
        ins: err * 0.33,
        del: err * 0.33,
        partial_frac: partial,
        unit_jitter: false, // a single TRUE allele has no unit jitter; error only
        unit_len,
        seed,
    }
}

/// True allele sequence (no error): left + mid + right.
fn truth(s: &Spec) -> Vec<u8> {
    let mut t = s.left.clone();
    t.extend_from_slice(&s.alleles[0].mid);
    t.extend_from_slice(&s.right);
    t
}

/// Returns (edit distance to truth, consensus len, truth len). edit=usize::MAX on error.
fn eval(s: &Spec) -> (usize, usize, usize) {
    let reads = synth::generate(s);
    let seqs = Read::seqs(&reads);
    let t = truth(s);
    match consensus(&seqs, 0, &cfg()) {
        Ok(c) => (
            synth::edit_distance(&c.sequence, &t),
            c.sequence.len(),
            t.len(),
        ),
        Err(_) => (usize::MAX, 0, t.len()),
    }
}

const MOTIFS: &[(&str, &[u8])] = &[
    ("CAG", b"CAG"),
    ("GAA", b"GAA"),
    ("AAAAG", b"AAAAG"),
    ("CATCAT", b"CATCAT"),
];

// ── Report sweep ──────────────────────────────────────────────────────────────

struct Row {
    label: String,
    edits: Vec<usize>, // per seed
    truth_len: usize,
}

fn report(title: &str, rows: &[Row]) {
    let n: usize = rows.iter().map(|r| r.edits.len()).sum();
    let perfect = rows
        .iter()
        .flat_map(|r| &r.edits)
        .filter(|&&e| e == 0)
        .count();
    let mean: f64 = rows
        .iter()
        .flat_map(|r| &r.edits)
        .map(|&e| e.min(1000) as f64)
        .sum::<f64>()
        / n as f64;
    println!("\n=== {title} ===  ({perfect}/{n} exact, mean edit {mean:.1})");
    for r in rows {
        let max = r.edits.iter().copied().max().unwrap_or(0);
        let mark = if max == 0 {
            "ok "
        } else if max <= 2 {
            "~  "
        } else {
            "XX "
        };
        println!(
            "  {mark} {:<30} truth_len {:>4}  edits {:?}",
            r.label, r.truth_len, r.edits
        );
    }
}

#[test]
#[ignore = "report; run with --ignored --nocapture"]
fn sweep_repeat_alleles() {
    let units = [10usize, 30, 60];
    let errs = [0.0, 0.02, 0.05, 0.08];
    let depths = [5usize, 10, 20];
    let seeds = [1u64, 2, 3];
    let mut rows = Vec::new();
    for (mn, m) in MOTIFS {
        for &u in &units {
            for &d in &depths {
                for &e in &errs {
                    let edits: Vec<usize> = seeds
                        .iter()
                        .map(|&s| eval(&spec(synth::repeat(m, u), m.len(), d, e, 0.0, s)).0)
                        .collect();
                    let tl = synth::LEFT.len() + m.len() * u + synth::RIGHT.len();
                    rows.push(Row {
                        label: format!("{mn}x{u} d{d} e{:.0}%", e * 100.0),
                        edits,
                        truth_len: tl,
                    });
                }
            }
        }
    }
    report("repeat single-allele: consensus vs truth", &rows);
}

#[test]
#[ignore = "report; run with --ignored --nocapture"]
fn sweep_nonrepeat_and_partial() {
    let seeds = [1u64, 2, 3];
    let mut rows = Vec::new();
    // Non-repeat (random) alleles — the "easy" control; must be near-perfect.
    for len in [80usize, 200, 400] {
        for &e in &[0.0, 0.05, 0.08] {
            let edits: Vec<usize> = seeds
                .iter()
                .map(|&s| {
                    let mid = synth::random_seq(len, 9999 + len as u64);
                    eval(&spec(mid, 1, 15, e, 0.0, s)).0
                })
                .collect();
            rows.push(Row {
                label: format!("random{len} d15 e{:.0}%", e * 100.0),
                edits,
                truth_len: synth::LEFT.len() + len + synth::RIGHT.len(),
            });
        }
    }
    // Partial reads (non-spanning) on a repeat.
    for &p in &[0.2f64, 0.4] {
        for &e in &[0.0, 0.05] {
            let edits: Vec<usize> = seeds
                .iter()
                .map(|&s| eval(&spec(synth::repeat(b"CAG", 30), 3, 20, e, p, s)).0)
                .collect();
            rows.push(Row {
                label: format!("CAGx30 d20 partial{:.0}% e{:.0}%", p * 100.0, e * 100.0),
                edits,
                truth_len: synth::LEFT.len() + 90 + synth::RIGHT.len(),
            });
        }
    }
    report("non-repeat control + partial reads", &rows);
}

// ── Must-hold assertions ────────────────────────────────────────────────────────

#[test]
fn clean_repeat_consensus_is_exact() {
    // A clean (0% error), adequately-covered repeat allele must reconstruct
    // exactly, across motifs, lengths, and seeds.
    for (name, m) in MOTIFS {
        for &u in &[10usize, 30, 60] {
            for seed in [1u64, 2, 3] {
                let (edit, clen, tlen) =
                    eval(&spec(synth::repeat(m, u), m.len(), 15, 0.0, 0.0, seed));
                assert_eq!(
                    edit, 0,
                    "{name}x{u} seed{seed}: clean consensus not exact (edit {edit}, len {clen} vs {tlen})"
                );
            }
        }
    }
}

#[test]
fn nonrepeat_consensus_near_exact_with_error() {
    // A non-repetitive allele at ONT error reconstructs near-exactly at depth 15
    // (error averages out; no periodic ambiguity). Small budget catches gross
    // regressions while tolerating a couple of residual high-error positions.
    for len in [80usize, 200] {
        for &e in &[0.0, 0.05, 0.08] {
            for seed in [1u64, 2, 3] {
                let mid = synth::random_seq(len, 777 + len as u64);
                let (edit, _, _) = eval(&spec(mid, 1, 15, e, 0.0, seed));
                let budget = 2 + (len as f64 * 0.02) as usize;
                assert!(
                    edit <= budget,
                    "random{len} e{:.0}% seed{seed}: edit {edit} > {budget}",
                    e * 100.0
                );
            }
        }
    }
}

#[test]
fn noisy_repeat_consensus_is_close() {
    // A repeat allele at 5% error, depth 15, must be within a small edit budget
    // (a few residual errors), not wildly wrong.
    for (name, m) in MOTIFS {
        for &u in &[30usize] {
            for seed in [1u64, 2, 3] {
                let (edit, _, tlen) =
                    eval(&spec(synth::repeat(m, u), m.len(), 15, 0.05, 0.0, seed));
                let budget = (tlen as f64 * 0.05) as usize + 3;
                assert!(
                    edit <= budget,
                    "{name}x{u} e5% seed{seed}: edit {edit} > budget {budget}"
                );
            }
        }
    }
}

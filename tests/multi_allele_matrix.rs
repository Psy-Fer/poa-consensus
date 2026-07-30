//! Multi-allele robustness matrix.
//!
//! Instead of one hand-picked test per scenario (which passes or fails by luck
//! and invites tuning to a fragile point), each CASE is swept across every axis
//! — motif, allele size/gap, depth, error, frequency, partial fraction, and
//! multiple RNG seeds — so we see the *gradient*: where a case is solid, where it
//! degrades, and why. The `*_sweep` tests (run with `--ignored --nocapture`)
//! print a per-axis report; the plain `#[test]`s assert must-hold points.
//!
//! Capture the full baseline of an engine (for before/after comparison when the
//! legacy machinery is removed):
//!   cargo test --release --features cli --test multi_allele_matrix -- --ignored --nocapture
//!   BASELINE_ENGINE=poa2 cargo test --release --features cli --test multi_allele_matrix -- --ignored --nocapture
//!
//! `BASELINE_ENGINE` selects the splitter: `legacy` (default; `consensus_multi`
//! → graph.rs) or `poa2` (`poa2::consensus_multi`, the WIP clean engine).

mod synth;

use poa_consensus::PoaConfig;
use synth::{AlleleSpec, Outcome, Read, Spec};

#[derive(Clone, Copy, PartialEq)]
enum Engine {
    Legacy,
    Poa2,
    Linkage,
    Hybrid,
}

fn engine() -> Engine {
    match std::env::var("BASELINE_ENGINE").as_deref() {
        Ok("poa2") => Engine::Poa2,
        Ok("linkage") => Engine::Linkage,
        Ok("hybrid") => Engine::Hybrid,
        _ => Engine::Legacy,
    }
}

fn engine_name() -> &'static str {
    match engine() {
        Engine::Legacy => "legacy",
        Engine::Poa2 => "poa2",
        Engine::Linkage => "linkage",
        Engine::Hybrid => "hybrid",
    }
}

fn cfg() -> PoaConfig {
    PoaConfig {
        min_reads: 3,
        min_allele_freq: 0.2,
        phasing_bubble_min_span: 10,
        ..PoaConfig::default()
    }
}

/// Run a read-set through the selected engine and return (per-group read indices,
/// per-group consensus length). Empty on error (scored as 0 alleles).
fn run_multi(reads: &[Read], config: &PoaConfig) -> (Vec<Vec<usize>>, Vec<usize>) {
    let seqs = Read::seqs(reads);
    let result = match engine() {
        Engine::Legacy => poa_consensus::consensus_multi(&seqs, 0, config),
        Engine::Poa2 => poa_consensus::poa2::consensus_multi(&seqs, config),
        Engine::Linkage => poa_consensus::poa2::linkage_consensus_multi(&seqs, config),
        Engine::Hybrid => poa_consensus::poa2::hybrid_consensus_multi(&seqs, config),
    };
    match result {
        Ok(alleles) => (
            alleles.iter().map(|a| a.read_indices.clone()).collect(),
            alleles.iter().map(|a| a.sequence.len()).collect(),
        ),
        Err(_) => (vec![], vec![]),
    }
}

fn eval(spec: &Spec) -> Outcome {
    let reads = synth::generate(spec);
    let (groups, lens) = run_multi(&reads, &cfg());
    synth::score_with_lens(&reads, spec, &groups, &lens)
}

/// Average over `seeds`, returning (all-seeds-correct, mean-misassign, got@seed0).
fn eval_seeds(mut make: impl FnMut(u64) -> Spec, seeds: &[u64]) -> (bool, f64, usize, bool) {
    let mut mis = 0.0;
    let mut all_ok = true;
    let mut any_split = false;
    for &s in seeds {
        let o = eval(&make(s));
        mis += o.misassign_rate;
        if !o.count_correct {
            all_ok = false;
        }
        if o.got_alleles > 1 {
            any_split = true;
        }
    }
    let got0 = eval(&make(seeds[0])).got_alleles;
    (all_ok, mis / seeds.len() as f64, got0, any_split)
}

// ── Spec builders (small alleles → fast builds; phasing logic still exercised) ─

const MOTIFS: &[(&str, &[u8])] = &[
    ("CAG", b"CAG"),
    ("GAA", b"GAA"),
    ("AAAAG", b"AAAAG"),
    ("CATCAT", b"CATCAT"),
];
const ERRS: &[f64] = &[0.0, 0.03, 0.06];
const SEEDS: &[u64] = &[1, 2, 3];

fn base_spec(alleles: Vec<AlleleSpec>, err: f64, unit_len: usize, partial: f64, seed: u64) -> Spec {
    Spec {
        left: synth::LEFT.to_vec(),
        right: synth::RIGHT.to_vec(),
        alleles,
        sub: err,
        ins: err,
        del: err,
        partial_frac: partial,
        unit_jitter: true,
        unit_len,
        seed,
    }
}

fn diploid_length(
    motif: &[u8],
    ua: usize,
    ub: usize,
    depth: usize,
    err: f64,
    partial: f64,
    seed: u64,
) -> Spec {
    base_spec(
        vec![
            AlleleSpec {
                mid: synth::repeat(motif, ua),
                n_reads: depth,
            },
            AlleleSpec {
                mid: synth::repeat(motif, ub),
                n_reads: depth,
            },
        ],
        err,
        motif.len(),
        partial,
        seed,
    )
}

fn single(motif: &[u8], units: usize, depth: usize, err: f64, seed: u64) -> Spec {
    base_spec(
        vec![AlleleSpec {
            mid: synth::repeat(motif, units),
            n_reads: depth,
        }],
        err,
        motif.len(),
        0.0,
        seed,
    )
}

/// Two SAME-length alleles differing only by substitutions (interruptions) at
/// fixed positions — the substitution-haplotype case (no length signal at all).
fn diploid_substitution(
    motif: &[u8],
    units: usize,
    n_subs: usize,
    depth: usize,
    err: f64,
    seed: u64,
) -> Spec {
    let a = synth::repeat(motif, units);
    let mut b = a.clone();
    // Space the interruptions across the repeat; flip the base deterministically.
    for k in 0..n_subs {
        let pos = (k + 1) * b.len() / (n_subs + 1);
        b[pos] = if b[pos] == b'A' { b'C' } else { b'A' };
    }
    base_spec(
        vec![
            AlleleSpec {
                mid: a,
                n_reads: depth,
            },
            AlleleSpec {
                mid: b,
                n_reads: depth,
            },
        ],
        err,
        motif.len(),
        0.0,
        seed,
    )
}

/// Majority allele + a minority (mosaic/subclonal) allele at `minor_frac`.
fn mosaic(
    motif: &[u8],
    ua: usize,
    ub: usize,
    total: usize,
    minor_frac: f64,
    err: f64,
    seed: u64,
) -> Spec {
    let minor = ((total as f64) * minor_frac).round() as usize;
    let major = total - minor;
    base_spec(
        vec![
            AlleleSpec {
                mid: synth::repeat(motif, ua),
                n_reads: major,
            },
            AlleleSpec {
                mid: synth::repeat(motif, ub),
                n_reads: minor,
            },
        ],
        err,
        motif.len(),
        0.0,
        seed,
    )
}

fn triallelic(motif: &[u8], u: (usize, usize, usize), depth: usize, err: f64, seed: u64) -> Spec {
    base_spec(
        vec![
            AlleleSpec {
                mid: synth::repeat(motif, u.0),
                n_reads: depth,
            },
            AlleleSpec {
                mid: synth::repeat(motif, u.1),
                n_reads: depth,
            },
            AlleleSpec {
                mid: synth::repeat(motif, u.2),
                n_reads: depth,
            },
        ],
        err,
        motif.len(),
        0.0,
        seed,
    )
}

// ── Report machinery ──────────────────────────────────────────────────────────

struct Row {
    label: String,
    got: usize,
    expected: usize,
    correct: bool,
    misassign: f64,
    note: String,
}

fn report(title: &str, rows: &[Row]) {
    let n = rows.len();
    let ok = rows.iter().filter(|r| r.correct).count();
    let mean_mis = rows.iter().map(|r| r.misassign).sum::<f64>() / n.max(1) as f64;
    println!(
        "\n=== [{}] {title} ===  ({ok}/{n} correct, mean misassign {:.1}%)",
        engine_name(),
        mean_mis * 100.0
    );
    for r in rows {
        let mark = if r.correct { "ok " } else { "XX " };
        println!(
            "  {mark} {:<34} got {}/{}  misassign {:>4.0}%  {}",
            r.label,
            r.got,
            r.expected,
            r.misassign * 100.0,
            r.note
        );
    }
}

// ── Sweeps (one per case) ─────────────────────────────────────────────────────

#[test]
#[ignore = "baseline report; run with --ignored --nocapture"]
fn case_single_allele_no_split() {
    let depths = [10usize, 20, 40];
    let mut rows = Vec::new();
    for (mn, m) in MOTIFS {
        for &err in ERRS {
            for &d in &depths {
                let gots: Vec<usize> = SEEDS
                    .iter()
                    .map(|&s| eval(&single(m, 30, d, err, s)).got_alleles)
                    .collect();
                let splits = gots.iter().filter(|&&g| g > 1).count();
                let errors = gots.iter().filter(|&&g| g == 0).count();
                // Correct = exactly one allele on every seed (not split, not errored).
                let correct = gots.iter().all(|&g| g == 1);
                rows.push(Row {
                    label: format!("{mn}x30 d{d} e{:.0}%", err * 100.0),
                    got: gots[0],
                    expected: 1,
                    correct,
                    misassign: 0.0,
                    note: format!("{splits}/3 over-split, {errors}/3 errored"),
                });
            }
        }
    }
    report(
        "single allele — must be exactly 1 (never over-split)",
        &rows,
    );
}

#[test]
#[ignore = "baseline report; run with --ignored --nocapture"]
fn case_diploid_length() {
    let gaps: &[(usize, usize)] = &[(15, 25), (20, 40), (30, 60), (40, 45)];
    let depths = [8usize, 20];
    let mut rows = Vec::new();
    for (mn, m) in MOTIFS {
        for &(ua, ub) in gaps {
            for &err in ERRS {
                for &d in &depths {
                    let (ok, mis, got, _) =
                        eval_seeds(|s| diploid_length(m, ua, ub, d, err, 0.0, s), SEEDS);
                    rows.push(Row {
                        label: format!("{mn} {ua}v{ub} d{d} e{:.0}%", err * 100.0),
                        got,
                        expected: 2,
                        correct: ok,
                        misassign: mis,
                        note: String::new(),
                    });
                }
            }
        }
    }
    report("diploid length variant", &rows);
}

#[test]
#[ignore = "baseline report; run with --ignored --nocapture"]
fn case_diploid_substitution() {
    let subs = [2usize, 4];
    let depths = [10usize, 20];
    let mut rows = Vec::new();
    for (mn, m) in MOTIFS {
        for &ns in &subs {
            for &err in ERRS {
                for &d in &depths {
                    let (ok, mis, got, _) =
                        eval_seeds(|s| diploid_substitution(m, 20, ns, d, err, s), SEEDS);
                    rows.push(Row {
                        label: format!("{mn}x20 {ns}subs d{d} e{:.0}%", err * 100.0),
                        got,
                        expected: 2,
                        correct: ok,
                        misassign: mis,
                        note: String::new(),
                    });
                }
            }
        }
    }
    report("diploid substitution (same length, SNP haplotypes)", &rows);
}

#[test]
#[ignore = "baseline report; run with --ignored --nocapture"]
fn case_mosaic() {
    let fracs = [0.1f64, 0.2, 0.3];
    let mut rows = Vec::new();
    for (mn, m) in &MOTIFS[..3] {
        for &f in &fracs {
            for &err in &[0.0f64, 0.03] {
                let (ok, mis, got, _) = eval_seeds(|s| mosaic(m, 20, 40, 30, f, err, s), SEEDS);
                rows.push(Row {
                    label: format!("{mn} 20/40 minor{:.0}% e{:.0}%", f * 100.0, err * 100.0),
                    got,
                    expected: 2,
                    correct: ok,
                    misassign: mis,
                    note: "minor allele recovery".into(),
                });
            }
        }
    }
    report("mosaic / subclonal minority allele", &rows);
}

#[test]
#[ignore = "baseline report; run with --ignored --nocapture"]
fn case_triallelic() {
    let sets = [(10usize, 20, 30), (15, 20, 25)];
    let mut rows = Vec::new();
    for (mn, m) in &MOTIFS[..3] {
        for &u in &sets {
            for &err in ERRS {
                let (ok, mis, got, _) = eval_seeds(|s| triallelic(m, u, 8, err, s), SEEDS);
                rows.push(Row {
                    label: format!("{mn} {}/{}/{} d8 e{:.0}%", u.0, u.1, u.2, err * 100.0),
                    got,
                    expected: 3,
                    correct: ok,
                    misassign: mis,
                    note: String::new(),
                });
            }
        }
    }
    report("triallelic (short/medium/long)", &rows);
}

// ── Must-hold CI guards (NOT ignored) ─────────────────────────────────────────
// A small, robust subset of the sweeps above, asserted as hard bounds so the
// headline multi-allele API (`consensus_multi` → hybrid) is actually exercised in
// CI. The full sweeps stay `#[ignore]`d reports; these are the regression floor.
// Run under the default engine (production hybrid via `poa_consensus::consensus_multi`).

#[test]
fn must_single_allele_never_oversplits() {
    // A single true allele must return exactly one consensus — never a spurious
    // second allele — across motifs, seeds, and clean/moderate error. This guards
    // the single-allele safety the whole phasing pipeline is gated on.
    for (mn, m) in MOTIFS {
        for &err in &[0.0f64, 0.03] {
            for &s in SEEDS {
                let o = eval(&single(m, 30, 20, err, s));
                assert_eq!(
                    o.got_alleles, 1,
                    "{mn}×30 single allele over-split at err={err} seed={s}: got {} alleles",
                    o.got_alleles
                );
            }
        }
    }
}

#[test]
fn must_clean_diploid_length_splits_cleanly() {
    // A clean, well-separated, well-covered length diploid (20 vs 40 units, 0%
    // error, depth 20) must split into exactly two alleles with no read
    // misassignment. This guards the core diploid-length calling path.
    for (mn, m) in MOTIFS {
        for &s in SEEDS {
            let o = eval(&diploid_length(m, 20, 40, 20, 0.0, 0.0, s));
            // Count is the headline must-hold; a few % read misassignment under
            // the per-read ±1-unit jitter is tolerable, so the misassign bound is
            // a gross-failure guard, not an exact-0.
            assert_eq!(
                o.got_alleles, 2,
                "{mn} 20v40 clean diploid should split into 2 (seed={s}): got {}",
                o.got_alleles
            );
            assert!(
                o.misassign_rate < 0.15,
                "{mn} 20v40 clean diploid grossly misphased: misassign {:.3} (seed={s})",
                o.misassign_rate
            );
        }
    }
}

#[test]
#[ignore = "baseline report; run with --ignored --nocapture"]
fn case_partial_reads() {
    let partials = [0.2f64, 0.4];
    let mut rows = Vec::new();
    for (mn, m) in &MOTIFS[..3] {
        for &p in &partials {
            for &err in &[0.0f64, 0.03] {
                let (ok, mis, got, _) =
                    eval_seeds(|s| diploid_length(m, 20, 40, 20, err, p, s), SEEDS);
                rows.push(Row {
                    label: format!(
                        "{mn} 20v40 d20 partial{:.0}% e{:.0}%",
                        p * 100.0,
                        err * 100.0
                    ),
                    got,
                    expected: 2,
                    correct: ok,
                    misassign: mis,
                    note: String::new(),
                });
            }
        }
    }
    report("diploid length + partial reads", &rows);
}

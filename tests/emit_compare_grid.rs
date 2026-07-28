//! Emits controlled single-allele read sets + ground truth across the HARD grid
//! (long repeats, high error, low depth, partials, all motifs + non-repeat), for
//! the 3-way caller comparison in `bench/compare_grid.py`. Reuses the one synth
//! generator so poa-consensus / abPOA / SPOA are scored on byte-identical reads.
//!
//! Run: COMPARE_GRID_DIR=/path cargo test --release --features cli \
//!        --test emit_compare_grid -- --ignored --nocapture

mod synth;

use std::fs;
use std::io::Write;
use synth::{AlleleSpec, Read, Spec};

fn spec(mid: Vec<u8>, unit_len: usize, depth: usize, err: f64, partial: f64, seed: u64) -> Spec {
    Spec {
        left: synth::LEFT.to_vec(),
        right: synth::RIGHT.to_vec(),
        alleles: vec![AlleleSpec {
            mid,
            n_reads: depth,
        }],
        // total error split ONT-style across sub/ins/del (they apply cumulatively).
        sub: err * 0.34,
        ins: err * 0.33,
        del: err * 0.33,
        partial_frac: partial,
        unit_jitter: false,
        unit_len,
        seed,
    }
}

fn truth(s: &Spec) -> Vec<u8> {
    let mut t = s.left.clone();
    t.extend_from_slice(&s.alleles[0].mid);
    t.extend_from_slice(&s.right);
    t
}

fn emit(dir: &str, name: &str, s: &Spec) {
    let reads = synth::generate(s);
    let mut fa = String::new();
    for (i, r) in reads.iter().enumerate() {
        fa.push_str(&format!(">r{i}\n{}\n", String::from_utf8_lossy(&r.seq)));
    }
    fs::write(format!("{dir}/{name}.reads.fa"), fa).unwrap();
    fs::write(format!("{dir}/{name}.truth"), truth(s)).unwrap();
}

#[test]
#[ignore = "emitter; run with COMPARE_GRID_DIR set and --ignored"]
fn emit_grid() {
    let dir = std::env::var("COMPARE_GRID_DIR").unwrap_or_else(|_| "/tmp/poa_compare_grid".into());
    fs::create_dir_all(&dir).unwrap();

    let motifs: &[(&str, &[u8])] = &[
        ("CAG", b"CAG"),
        ("GAA", b"GAA"),
        ("AAAAG", b"AAAAG"),
        ("CATCAT", b"CATCAT"),
    ];
    let units = [30usize, 60, 100];
    let depths = [5usize, 10, 20];
    let errs = [0.0f64, 0.05, 0.08, 0.12];
    let seeds = [1u64, 2];

    let mut manifest = String::new();
    for (mn, m) in motifs {
        for &u in &units {
            for &d in &depths {
                for &e in &errs {
                    for &s in &seeds {
                        let name = format!("{mn}_u{u}_d{d}_e{:.0}_s{s}", e * 100.0);
                        emit(
                            &dir,
                            &name,
                            &spec(synth::repeat(m, u), m.len(), d, e, 0.0, s),
                        );
                        manifest.push_str(&name);
                        manifest.push('\n');
                    }
                }
            }
        }
    }
    // Non-repeat controls.
    for &len in &[200usize, 400] {
        for &e in &[0.0f64, 0.08, 0.12] {
            for &s in &seeds {
                let name = format!("random{len}_d15_e{:.0}_s{s}", e * 100.0);
                let mid = synth::random_seq(len, 4242 + len as u64);
                emit(&dir, &name, &spec(mid, 1, 15, e, 0.0, s));
                manifest.push_str(&name);
                manifest.push('\n');
            }
        }
    }
    // Partial-read stress on a repeat.
    for &p in &[0.3f64] {
        for &e in &[0.05f64, 0.08] {
            for &s in &seeds {
                let name = format!("CAG_u60_d20_partial{:.0}_e{:.0}_s{s}", p * 100.0, e * 100.0);
                emit(
                    &dir,
                    &name,
                    &spec(synth::repeat(b"CAG", 60), 3, 20, e, p, s),
                );
                manifest.push_str(&name);
                manifest.push('\n');
            }
        }
    }
    let mut mf = fs::File::create(format!("{dir}/manifest.txt")).unwrap();
    mf.write_all(manifest.as_bytes()).unwrap();
    let n = manifest.lines().count();
    eprintln!("emitted {n} configs to {dir}");
    let _ = Read::seqs; // silence unused-import lints across cfgs
}

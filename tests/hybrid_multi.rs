//! Curated must-hold regression tests for the hybrid multi-allele engine
//! (`poa2::hybrid_consensus_multi` = structural proposal + linkage confirmation).
//! These are named per case so a failure pinpoints the exact scenario; the broad
//! gradient lives in `multi_allele_matrix.rs`.

mod synth;

use poa_consensus::PoaConfig;
use poa_consensus::poa2::hybrid_consensus_multi;
use synth::{AlleleSpec, Read, Spec};

fn cfg() -> PoaConfig {
    PoaConfig {
        min_reads: 3,
        min_allele_freq: 0.2,
        phasing_bubble_min_span: 10,
        ..PoaConfig::default()
    }
}

fn n_alleles(spec: &Spec) -> usize {
    let reads = synth::generate(spec);
    let seqs = Read::seqs(&reads);
    hybrid_consensus_multi(&seqs, &cfg())
        .map(|a| a.len())
        .unwrap_or(0)
}

fn diploid(motif: &[u8], ua: usize, ub: usize, depth: usize, err: f64, seed: u64) -> Spec {
    Spec {
        left: synth::LEFT.to_vec(),
        right: synth::RIGHT.to_vec(),
        alleles: vec![
            AlleleSpec {
                mid: synth::repeat(motif, ua),
                n_reads: depth,
            },
            AlleleSpec {
                mid: synth::repeat(motif, ub),
                n_reads: depth,
            },
        ],
        sub: err,
        ins: err,
        del: err,
        partial_frac: 0.0,
        unit_jitter: true,
        unit_len: motif.len(),
        seed,
    }
}

fn single(motif: &[u8], units: usize, depth: usize, err: f64, seed: u64) -> Spec {
    Spec {
        left: synth::LEFT.to_vec(),
        right: synth::RIGHT.to_vec(),
        alleles: vec![AlleleSpec {
            mid: synth::repeat(motif, units),
            n_reads: depth,
        }],
        sub: err,
        ins: err,
        del: err,
        partial_frac: 0.0,
        unit_jitter: true,
        unit_len: motif.len(),
        seed,
    }
}

#[test]
fn hybrid_clean_single_allele_never_splits() {
    // Safety: a clean single allele must stay one allele across motifs/seeds.
    for (name, m) in [
        ("CAG", &b"CAG"[..]),
        ("GAA", &b"GAA"[..]),
        ("AAAAG", &b"AAAAG"[..]),
    ] {
        for seed in [1u64, 2, 3] {
            let k = n_alleles(&single(m, 30, 20, 0.0, seed));
            assert_eq!(
                k, 1,
                "{name} seed{seed}: clean single allele split into {k}"
            );
        }
    }
}

#[test]
fn hybrid_noisy_single_allele_rarely_splits() {
    // Safety under error: a single allele at 3% error should almost never be
    // fabricated into two (the confirmation gate rejects noise sub-splits).
    let mut splits = 0;
    let mut total = 0;
    for (_, m) in [
        ("CAG", &b"CAG"[..]),
        ("GAA", &b"GAA"[..]),
        ("AAAAG", &b"AAAAG"[..]),
    ] {
        for seed in [1u64, 2, 3, 4] {
            total += 1;
            if n_alleles(&single(m, 30, 20, 0.03, seed)) > 1 {
                splits += 1;
            }
        }
    }
    assert!(
        splits <= 1,
        "single allele over-split in {splits}/{total} noisy configs"
    );
}

#[test]
fn hybrid_clean_diploid_length_splits_into_two() {
    // Sensitivity: a clean, well-separated length variant must split into two.
    for (name, m) in [("CAG", &b"CAG"[..]), ("GAA", &b"GAA"[..])] {
        for seed in [1u64, 2, 3] {
            let k = n_alleles(&diploid(m, 15, 30, 12, 0.0, seed));
            assert_eq!(k, 2, "{name} seed{seed}: clean diploid gave {k} alleles");
        }
    }
}

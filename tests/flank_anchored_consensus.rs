//! Regression tests for the flank-anchoring API (`consensus_flanked` /
//! `consensus_multi_flanked`). The win is on partial-read-heavy loci; on
//! fully-spanning loci it must be neutral. Output is always the repeat region
//! BETWEEN the flanks (flanks excluded).

mod synth;

use poa_consensus::{PoaConfig, consensus_flanked, consensus_multi_flanked};
use synth::{AlleleSpec, Spec};

const L: &[u8] = synth::LEFT;
const R: &[u8] = synth::RIGHT;

fn cfg() -> PoaConfig {
    PoaConfig {
        min_reads: 3,
        min_allele_freq: 0.2,
        phasing_bubble_min_span: 10,
        ..PoaConfig::default()
    }
}

fn spec(alleles: Vec<AlleleSpec>, err: f64, partial: f64, seed: u64) -> Spec {
    Spec {
        left: L.to_vec(),
        right: R.to_vec(),
        alleles,
        sub: err,
        ins: err,
        del: err,
        partial_frac: partial,
        unit_jitter: false,
        unit_len: 3,
        seed,
    }
}

fn contains(hay: &[u8], needle: &[u8]) -> bool {
    hay.windows(needle.len()).any(|w| w == needle)
}

#[test]
fn flanked_output_excludes_flanks_on_spanning_reads() {
    let mid = synth::repeat(b"CAG", 20); // 60 bp truth
    let reads = synth::generate(&spec(
        vec![AlleleSpec {
            mid: mid.clone(),
            n_reads: 12,
        }],
        0.01,
        0.0,
        1,
    ));
    let seqs = synth::Read::seqs(&reads);
    let c = consensus_flanked(&seqs, L, R, &cfg()).expect("spanning consensus");

    // Output is the repeat region only: near 60 bp, and it must NOT contain the flanks.
    assert!(
        (57..=63).contains(&c.sequence.len()),
        "repeat-only length ~60; got {}",
        c.sequence.len()
    );
    assert!(
        !contains(&c.sequence, L),
        "left flank must be trimmed from the output"
    );
    assert!(
        !contains(&c.sequence, R),
        "right flank must be trimmed from the output"
    );
    // Per-position fields stay aligned to the sliced sequence.
    assert_eq!(c.coverage.len(), c.sequence.len());
    assert_eq!(c.path_weights.len(), c.sequence.len());
    // Close to the true repeat.
    assert!(
        synth::edit_distance(&c.sequence, &mid) <= 4,
        "edit {} to CAG x20",
        synth::edit_distance(&c.sequence, &mid)
    );
}

#[test]
fn flanked_recovers_repeat_on_partial_heavy_reads() {
    // 40% of reads are partial (truncated, missing a flank) — the case that
    // distorts a raw consensus. Anchoring drops them and recovers the repeat.
    let mid = synth::repeat(b"CAG", 20);
    let mut ok = 0;
    for seed in 1..=5u64 {
        let reads = synth::generate(&spec(
            vec![AlleleSpec {
                mid: mid.clone(),
                n_reads: 20,
            }],
            0.01,
            0.4,
            seed,
        ));
        let seqs = synth::Read::seqs(&reads);
        let c = consensus_flanked(&seqs, L, R, &cfg()).expect("partial-heavy consensus");
        assert!(!contains(&c.sequence, L), "output must be repeat-only");
        if synth::edit_distance(&c.sequence, &mid) <= 4 {
            ok += 1;
        }
    }
    assert!(
        ok >= 4,
        "flank-anchoring should recover CAG x20 on partial-heavy sets ({ok}/5 seeds)"
    );
}

#[test]
fn flanked_falls_back_without_enough_spanning_reads() {
    // Almost everything is partial → too few reads span to anchor. Must not panic
    // and must still return a consensus (raw fallback).
    let mid = synth::repeat(b"CAG", 20);
    let reads = synth::generate(&spec(vec![AlleleSpec { mid, n_reads: 6 }], 0.01, 0.95, 7));
    let seqs = synth::Read::seqs(&reads);
    let c = consensus_flanked(&seqs, L, R, &cfg()).expect("fallback consensus");
    assert!(
        !c.sequence.is_empty(),
        "fallback still produces a consensus"
    );
}

#[test]
fn multi_flanked_read_indices_map_to_original_reads() {
    // Diploid length variant (CAG x20 vs x28), partial-heavy. The anchored path
    // runs on a filtered subset; read_indices must be remapped back to indices
    // into the ORIGINAL reads slice.
    let a = synth::repeat(b"CAG", 20);
    let b = synth::repeat(b"CAG", 28);
    let reads = synth::generate(&spec(
        vec![
            AlleleSpec {
                mid: a,
                n_reads: 15,
            },
            AlleleSpec {
                mid: b,
                n_reads: 15,
            },
        ],
        0.01,
        0.3,
        3,
    ));
    let seqs = synth::Read::seqs(&reads);
    let alleles = consensus_multi_flanked(&seqs, L, R, &cfg()).expect("multi flanked");
    assert!(!alleles.is_empty());
    for allele in &alleles {
        assert!(
            !contains(&allele.sequence, L),
            "allele output must be repeat-only"
        );
        for &idx in &allele.read_indices {
            assert!(
                idx < reads.len(),
                "read index {idx} out of range (remap failed)"
            );
        }
    }
}

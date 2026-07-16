//! Regression for multi-allele read-partition contamination on periodic
//! (homogeneous) tandem repeats.
//!
//! `find_structural_bubbles` used to measure an arm's span by walking a
//! strictly unbranched chain, stopping at the first fork or reconvergence no
//! matter how small (a single read's indel error). In a periodic repeat,
//! ordinary per-read indel noise (insertions/deletions, not substitutions --
//! the dominant error mode on ONT/HiFi periodic-repeat data, and the kind
//! that actually creates phase-shift ambiguity: matching-vs-inserting at any
//! one column is close to a coin flip) scatters many such tiny disruptions
//! throughout the true length-differentiating region, so no single arm ever
//! measured long enough to be classified "structural" -- `consensus_multi`
//! fell back to a single coincidental SNP-bubble column to decide which
//! allele each read belonged to. That single column's split does not
//! reliably track a read's true (much longer) haplotype identity in a
//! periodic repeat -- confirmed on real GAA×30/GAA×100 HiFi data via
//! independent ground-truth read tracing: ~29% of reads ended up in the
//! wrong allele's output.
//!
//! This test reproduces the same *shape* of problem synthetically (a
//! periodic repeat, two very different lengths, ordinary scattered per-read
//! indel noise -- no real patient data, same discipline as
//! `tests/band_width_zero_unbanded.rs`) and checks that `consensus_multi`
//! recovers exactly two allele groups with no bulk cross-contamination.
//!
//! Confirmed this reproduces the pre-fix failure directly (not just
//! hypothesized): built against the pre-fix `src/graph.rs` (committed at
//! `72e201d`, before this round's noise-tolerant span measurement,
//! `phasing_groups` redesign, and `validate_and_merge_groups`), this exact
//! read set produces two allele groups that are each ~30% cross-contaminated
//! (14 correct + 6 wrong-origin reads in one group, 6 + 14 in the other --
//! matching the real-data magnitude closely). With the fix restored, the
//! same read set produces one clean group (0% contamination) and one group
//! with a single stray read (3/23 ~= 13%, ordinary noise-absorption, not the
//! bulk bug).

use poa_consensus::{AlignmentMode, ConsensusMode, PoaConfig};

/// Deterministic (no external RNG dependency) scattered per-read indel noise:
/// inserts or deletes single bases in `seq` at positions that vary by
/// `read_idx`, cycling through the repeat region so different reads disrupt
/// different columns -- the same shape of noise (independent, scattered
/// single-base indels) that fragments a periodic repeat's true structural
/// bubble into many tiny, individually-insignificant ones on real
/// sequencing data.
fn with_scattered_noise(seq: &[u8], read_idx: usize, region: std::ops::Range<usize>) -> Vec<u8> {
    let span = region.end - region.start;
    let mut out = seq.to_vec();

    // First indel: alternate delete/insert by read parity, position cycles
    // through the region so different reads disrupt different columns.
    let pos = region.start + (read_idx * 7) % span;
    if read_idx % 2 == 0 {
        if pos < out.len() {
            out.remove(pos);
        }
    } else {
        out.insert(pos.min(out.len()), b'A');
    }

    // Second indel: opposite direction, different cycling stride, so noise
    // accumulates from more than one independent source per read (closer to
    // the layered error profile of real long-read sequencing).
    let pos2 = region.start + (read_idx * 13 + 5) % span;
    if read_idx % 2 == 0 {
        out.insert(pos2.min(out.len()), b'A');
    } else if pos2 < out.len() {
        out.remove(pos2);
    }

    // Third, sparser indel affecting only 1-in-5 reads: enough extra noise to
    // push the pre-fix fallback into bulk (~30%) cross-contamination, matching
    // the real-data magnitude, without which this synthetic case stays too
    // clean even on the pre-fix code.
    let pos3 = region.start + (read_idx * 17 + 3) % span;
    if read_idx % 5 == 0 && pos3 < out.len() {
        out.remove(pos3);
    }

    out
}

#[test]
fn structural_phasing_no_contamination_on_noisy_periodic_repeat() {
    let left = b"ACGTACGTCGATCGATTAGCTAGCGCTAGCTA"; // 32bp unique flank, no GAA
    let right = b"ATCGATCGCGATCGATTAGCTAGCTGCATGCA"; // 32bp unique flank, no GAA

    let short_units = 15usize; // 45bp
    let long_units = 40usize; // 120bp
    let short_mid: Vec<u8> = b"GAA".repeat(short_units);
    let long_mid: Vec<u8> = b"GAA".repeat(long_units);

    let make = |mid: &[u8]| -> Vec<u8> {
        let mut r = left.to_vec();
        r.extend_from_slice(mid);
        r.extend_from_slice(right);
        r
    };
    let short_read = make(&short_mid);
    let long_read = make(&long_mid);

    let short_region = left.len()..(left.len() + short_mid.len());
    let long_region = left.len()..(left.len() + long_mid.len());

    let n_per_allele = 20;
    let mut reads: Vec<Vec<u8>> = Vec::new();
    for i in 0..n_per_allele {
        reads.push(with_scattered_noise(&short_read, i, short_region.clone()));
    }
    for i in 0..n_per_allele {
        reads.push(with_scattered_noise(&long_read, i, long_region.clone()));
    }

    let cfg = PoaConfig {
        band_width: 0,
        adaptive_band: true,
        min_reads: 2,
        alignment_mode: AlignmentMode::SemiGlobal,
        consensus_mode: ConsensusMode::HeaviestPath,
        ..PoaConfig::default()
    };

    let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    let alleles = poa_consensus::consensus_multi(&refs, 0, &cfg).unwrap();

    assert_eq!(
        alleles.len(),
        2,
        "expected exactly 2 allele groups for a genuinely biallelic periodic repeat, got {}",
        alleles.len()
    );

    // Each output allele's read_indices should be overwhelmingly dominated by
    // reads from ONE true origin. A single stray read whose own scattered
    // noise happened to make it genuinely ambiguous (no recorded arm choice,
    // or compatible with 2+ clusters) legitimately folds into the largest
    // group by design -- that is not the bug this test guards against. The
    // bug this guards against is bulk misassignment: on real GAA×30/GAA×100
    // data, the pre-fix failure put ~29% of reads (11-15 of 52) in the wrong
    // allele's group; on this synthetic read set the pre-fix failure measures
    // an almost identical ~30% (6 of 20) per group. A generous 15%
    // cross-origin tolerance per group is well below that bulk-failure
    // magnitude but comfortably above ordinary single-read noise-absorption
    // (post-fix, this read set leaves one clean group and one group with a
    // single stray read, 3 of 23 ~= 13%). (Read indices 0..n_per_allele are
    // short-origin; n_per_allele..2*n_per_allele are long-origin.)
    const MAX_CONTAMINATION_FRACTION: f64 = 0.15;
    for allele in &alleles {
        let n_short_origin = allele
            .read_indices
            .iter()
            .filter(|&&idx| idx < n_per_allele)
            .count();
        let n_long_origin = allele
            .read_indices
            .iter()
            .filter(|&&idx| idx >= n_per_allele)
            .count();
        let total = n_short_origin + n_long_origin;
        let minority = n_short_origin.min(n_long_origin);
        assert!(
            total == 0 || (minority as f64 / total as f64) <= MAX_CONTAMINATION_FRACTION,
            "allele group mixes short-origin ({n_short_origin}) and long-origin \
             ({n_long_origin}) reads far beyond ordinary noise-absorption levels \
             -- looks like the bulk cross-contamination bug, not stray noise"
        );
    }

    let mut lens: Vec<usize> = alleles.iter().map(|c| c.sequence.len()).collect();
    lens.sort_unstable();
    let expected_short = left.len() + short_mid.len() + right.len();
    let expected_long = left.len() + long_mid.len() + right.len();
    // Allow a small tolerance for the scattered indel noise itself: each
    // individual read's length wobbles by a couple of bases, but the
    // consensus's own length is a majority-vote outcome and should still
    // land on the true repeat-unit count.
    assert!(
        lens[0].abs_diff(expected_short) <= 3,
        "short allele length {} far from expected {expected_short}",
        lens[0]
    );
    assert!(
        lens[1].abs_diff(expected_long) <= 3,
        "long allele length {} far from expected {expected_long}",
        lens[1]
    );
}

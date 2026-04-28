//! Comparison of poa-consensus (affine-gap) output against the original
//! linear-gap reference implementation in `ref/poa.rs`.
//!
//! **Expected divergence:** The two implementations use different gap penalty
//! models (affine vs. linear). On indel-heavy inputs they may produce
//! different consensus sequences. On SNP-only inputs they should agree.
//! Tests below document where they agree, where they differ, and why.

// ref/poa.rs tests call `crate::misc::medoid_index`; provide a stub so it
// compiles when included as a path module here.
mod misc {
    pub fn medoid_index(reads: &[Vec<u8>]) -> usize {
        let n = reads.len();
        if n == 0 {
            return 0;
        }
        (0..n)
            .min_by_key(|&i| {
                reads
                    .iter()
                    .map(|r| r.len().abs_diff(reads[i].len()))
                    .sum::<usize>()
            })
            .unwrap_or(0)
    }
}

#[path = "../ref/poa.rs"]
mod ref_poa;

use poa_consensus::{PoaConfig, consensus};

// ── Helpers ───────────────────────────────────────────────────────────────────

fn run_ref(reads: &[&[u8]], seed_idx: usize) -> Vec<u8> {
    let owned: Vec<Vec<u8>> = reads.iter().map(|r| r.to_vec()).collect();
    ref_poa::get_banded_consensus(&owned, seed_idx)
}

fn run_new(reads: &[&[u8]], seed_idx: usize) -> Vec<u8> {
    consensus(reads, seed_idx, &PoaConfig::default())
        .expect("consensus failed")
        .sequence
}

fn ascii(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

// ── Agreement cases ───────────────────────────────────────────────────────────

/// Both implementations must agree on a clean SNP-only dataset.
#[test]
fn agree_snp_only() {
    let reads: &[&[u8]] = &[
        b"ACGTACGTACGTACGT",
        b"ACGTACGTACGTACGT",
        b"ACGTACGTACGTACGT",
        b"ACGTACGTAGGTACGT", // SNP at position 8
        b"ACGTACGTACGTACGT",
    ];
    let ref_out = run_ref(reads, 0);
    let new_out = run_new(reads, 0);
    assert_eq!(
        ascii(&ref_out),
        ascii(&new_out),
        "SNP-only: ref={} new={}",
        ascii(&ref_out),
        ascii(&new_out)
    );
}

/// Identical reads — both return the read unchanged.
#[test]
fn agree_identical_reads() {
    let reads: &[&[u8]] = &[
        b"CATCATCATCAT",
        b"CATCATCATCAT",
        b"CATCATCATCAT",
        b"CATCATCATCAT",
        b"CATCATCATCAT",
    ];
    let ref_out = run_ref(reads, 0);
    let new_out = run_new(reads, 0);
    assert_eq!(
        ascii(&ref_out),
        ascii(&new_out),
        "identical: ref={} new={}",
        ascii(&ref_out),
        ascii(&new_out)
    );
    assert_eq!(ascii(&new_out), "CATCATCATCAT");
}

/// Single outlier with a SNP — majority wins in both.
#[test]
fn agree_single_outlier_snp() {
    let reads: &[&[u8]] = &[
        b"GCTAGCTAGCTA",
        b"GCTAGCTAGCTA",
        b"GCTAGCTAGCTA",
        b"GCTAGCTAGCTA",
        b"GCTAGCAAGCTA", // outlier
    ];
    let ref_out = run_ref(reads, 0);
    let new_out = run_new(reads, 0);
    assert_eq!(ascii(&new_out), "GCTAGCTAGCTA", "new consensus wrong");
    assert_eq!(
        ascii(&ref_out),
        ascii(&new_out),
        "ref={} new={}",
        ascii(&ref_out),
        ascii(&new_out)
    );
}

// ── Documented divergence ─────────────────────────────────────────────────────

/// On length-variable inputs (insertions/deletions), affine-gap penalties
/// change which path is preferred. We verify that BOTH implementations pick
/// the majority length, even if the sequences differ slightly.
#[test]
fn diverge_length_variation_both_pick_majority_length() {
    // 9 reads at 12 bp, 1 outlier at 18 bp
    let reads: &[&[u8]] = &[
        b"CATCATCATCAT",
        b"CATCATCATCAT",
        b"CATCATCATCAT",
        b"CATCATCATCAT",
        b"CATCATCATCAT",
        b"CATCATCATCAT",
        b"CATCATCATCAT",
        b"CATCATCATCAT",
        b"CATCATCATCAT",
        b"CATCATCATCATCATCAT", // outlier: 18 bp
    ];
    let ref_len = run_ref(reads, 0).len();
    let new_len = run_new(reads, 0).len();

    // Both must select the 12-bp majority, not the 18-bp outlier.
    assert_eq!(ref_len, 12, "ref chose wrong length: {ref_len}");
    assert_eq!(new_len, 12, "new chose wrong length: {new_len}");
}

/// With affine gaps, the new implementation handles a single-bp insertion
/// outlier more cleanly than the linear-gap ref. We just assert the new
/// result is correct; the ref result is logged for documentation.
#[test]
fn diverge_single_insertion_outlier_new_is_correct() {
    let reads: &[&[u8]] = &[
        b"ATCGATCG",
        b"ATCGATCG",
        b"ATCGATCG",
        b"ATCGATCG",
        b"ATCGATCG",
        b"ATCGAATCG", // 1-bp insertion outlier
    ];
    let new_out = run_new(reads, 0);
    // New impl with affine gaps should not include the spurious insertion.
    assert_eq!(ascii(&new_out), "ATCGATCG", "new={}", ascii(&new_out));

    let ref_out = run_ref(reads, 0);
    // Document ref behaviour (may or may not include the insertion).
    eprintln!(
        "ref_insertion_outlier: ref={} new={}",
        ascii(&ref_out),
        ascii(&new_out)
    );
}

// ── Ground-truth recovery ─────────────────────────────────────────────────────

/// Verify poa-consensus recovers the ground truth on the simple fixture reads.
/// This is the primary correctness sanity check.
#[test]
fn ground_truth_simple_fixture() {
    let ground_truth = b"ACGTACGTACGTACGTACGTACGTACGT";
    let reads: &[&[u8]] = &[
        b"ACGTACGTACGTACGTACGTACGTACGT",
        b"ACGTACGTACGTACGTACGTACGTACGT",
        b"ACGTACGTACTTACGTACGTACGTACGT",
        b"ACGTACGTACGTACGTACGTACGTACGT",
        b"ACGTACGTACGTACGTACGTACGGACGT",
        b"ACGTACGTACGTACGTACGTACGTACGT",
        b"ACGTAGGTACGTACGTACGTACGTACGT",
        b"ACGTACGTACGTACGTACGTACGTACGT",
        b"ACGTACGTACGTACGTACGTACGTACGT",
        b"ACGTACGTACGTACGTACGTACGTACGT",
    ];
    let out = run_new(reads, 0);
    assert_eq!(
        out.as_slice(),
        ground_truth,
        "expected {}, got {}",
        ascii(ground_truth),
        ascii(&out)
    );
}

/// Verify poa-consensus recovers the ground truth on the noisy fixture reads.
#[test]
fn ground_truth_noisy_fixture() {
    let ground_truth = b"GCTAGCTAGCTAGCTA";
    let reads: &[&[u8]] = &[
        b"GCTAGCTAGCTAGCTA",
        b"GCTAGCTAGCTAGCTA",
        b"GCTAACTAGCTAGCTA",
        b"GCTAGCTAGCTAGCTA",
        b"GCTAGCTAGCTAGCTA",
        b"GCTAGCTAGCTAGCTA",
        b"GCTAGCTAGCTAGCTA",
        b"GCTAGCTAGCTATCTA",
        b"GCTAGCTAGCTAGCTA",
        b"GCTAGCTAGCTAGCTA",
    ];
    let out = run_new(reads, 0);
    assert_eq!(
        out.as_slice(),
        ground_truth,
        "expected {}, got {}",
        ascii(ground_truth),
        ascii(&out)
    );
}

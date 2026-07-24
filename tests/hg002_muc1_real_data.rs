//! Real-data integration tests using the MUC1 flVNTR locus (chr1:154327978-154330908, CHM13).
//!
//! Data provenance (all public — HG002 / GIAB, no patient data):
//!   Reads pulled from the HG002 ONT alignment with `bedpull` over the MUC1
//!   flVNTR locus, then split by haplotype. Ground-truth unit counts come from
//!   the HG002 hifiasm assembly of the same locus.
//!
//!   hap1 -- 65-unit allele, reads ~4130 bp (26 reads)
//!   hap2 -- 77-unit allele, reads ~4868 bp (16 reads)
//!
//! The MUC1 flVNTR is a ~60 bp GC-rich tandem repeat with divergent unit
//! variants. At ONT error rates this produces many small phase-registration
//! bubbles, which is the regime where a length-biased consensus extraction
//! over-calls the repeat length (see design/vntr_overcall_delete_edge_visibility.md).
//!
//! hap2 is the regression guard: the correct consensus length is 4868 bp
//! (independently confirmed — abPOA, SPOA, and poa-consensus 0.2.1 all agree),
//! whereas the length-biased heaviest-path over-calls to ~5168 bp (+5 units).
//! hap1 is a positive control that must stay correct.
//!
//! Run with:
//!   cargo test --test hg002_muc1_real_data -- --nocapture

use poa_consensus::{PoaConfig, SeedSelection, consensus, select_seed};

fn parse_fasta(s: &str) -> Vec<Vec<u8>> {
    let mut reads: Vec<Vec<u8>> = Vec::new();
    let mut current: Vec<u8> = Vec::new();
    for line in s.lines() {
        if line.starts_with('>') {
            if !current.is_empty() {
                reads.push(std::mem::take(&mut current));
            }
        } else {
            current.extend_from_slice(line.trim().as_bytes());
        }
    }
    if !current.is_empty() {
        reads.push(current);
    }
    reads
}

fn median_len(reads: &[Vec<u8>]) -> usize {
    let mut lens: Vec<usize> = reads.iter().map(|r| r.len()).collect();
    lens.sort_unstable();
    lens[lens.len() / 2]
}

/// Build a single-allele consensus with the default (CLI-equivalent) config and
/// the auto-selected seed, returning its length.
fn consensus_len(reads: &[Vec<u8>]) -> usize {
    let slices: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
    let seed = select_seed(&slices, &SeedSelection::Auto).expect("seed selection");
    let cons = consensus(&slices, seed, &PoaConfig::default()).expect("consensus");
    cons.sequence.len()
}

// MUC1 flVNTR unit is ~60 bp; allow ~2 units of tolerance for ONT noise at the
// repeat/flank boundaries. The baseline over-call is +5 units (300 bp), so this
// tolerance decisively separates correct from over-called.
const UNIT: usize = 60;
const TOL: usize = 2 * UNIT;

// ── hap2: 77-unit allele — the over-call regression guard ─────────────────────

#[test]
fn hg002_muc1_hap2_length_not_overcalled() {
    let reads = parse_fasta(include_str!("fixtures/hg002_muc1_hap2.fa"));
    assert_eq!(reads.len(), 16, "fixture read count");

    // Ground truth: 77 units; correct consensus length 4868 bp (abPOA / SPOA /
    // poa-consensus 0.2.1 all agree). Median read length is the same.
    const TRUTH: usize = 4868;
    let median = median_len(&reads);
    let len = consensus_len(&reads);

    assert!(
        len.abs_diff(TRUTH) <= TOL,
        "hap2 consensus length {len} bp is not within {TOL} bp of the 77-unit \
         truth {TRUTH} bp (median read {median}); a large positive deviation is \
         the length-biased heaviest-path over-call this test guards against",
    );
}

// ── hap1: 65-unit allele — positive control (must stay correct) ───────────────

#[test]
fn hg002_muc1_hap1_length_correct() {
    let reads = parse_fasta(include_str!("fixtures/hg002_muc1_hap1.fa"));
    assert_eq!(reads.len(), 26, "fixture read count");

    // Ground truth: 65 units; correct consensus length ~4130 bp. This allele is
    // called correctly by the baseline too — the guard is that the fix must not
    // regress it.
    const TRUTH: usize = 4130;
    let median = median_len(&reads);
    let len = consensus_len(&reads);

    assert!(
        len.abs_diff(TRUTH) <= TOL,
        "hap1 consensus length {len} bp is not within {TOL} bp of the 65-unit \
         truth {TRUTH} bp (median read {median})",
    );
}

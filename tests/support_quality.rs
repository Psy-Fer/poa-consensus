//! Per-base support quality (`Consensus::support_quality`): one value per
//! consensus base on every entry point, high where the reads agree and low
//! where they split.

use poa_consensus::analysis::low_support_regions;
use poa_consensus::{
    ConsensusMode, PoaConfig, bridged_consensus, consensus, consensus_flanked, consensus_multi,
};

// The flanks must not continue the CA repeat: a left flank ending in A makes
// deleting "AC" across the boundary equivalent to deleting a "CA" unit, so the
// ambiguity (and low support) would rightly extend one base into the flank.
const LEFT: &str = "GATCGTACGTTAGCATGCTG";
const RIGHT: &str = "TTGCAGTCAGGATCCATGCA";

fn locus(units: usize) -> Vec<u8> {
    format!("{LEFT}{}{RIGHT}", "CA".repeat(units)).into_bytes()
}

fn slices(v: &[Vec<u8>]) -> Vec<&[u8]> {
    v.iter().map(|r| r.as_slice()).collect()
}

#[test]
fn identical_reads_are_solid_everywhere() {
    let reads: Vec<Vec<u8>> = (0..10).map(|_| locus(20)).collect();
    let c = consensus(&slices(&reads), 0, &PoaConfig::default()).unwrap();
    assert_eq!(c.sequence, locus(20));
    assert_eq!(c.support_quality.len(), c.sequence.len());
    assert!(
        c.support_quality.iter().all(|&q| q == 33),
        "{:?}",
        c.support_quality
    );
    assert!(low_support_regions(&c, 10).is_empty());
}

#[test]
fn split_repeat_length_is_marked_and_flanks_stay_solid() {
    // 6 reads say (CA)20, 4 say (CA)19: the call is the majority length, but
    // the reads disagree on it, and that should show inside the repeat only.
    let mut reads: Vec<Vec<u8>> = (0..6).map(|_| locus(20)).collect();
    reads.extend((0..4).map(|_| locus(19)));
    let c = consensus(&slices(&reads), 0, &PoaConfig::default()).unwrap();
    assert_eq!(c.sequence, locus(20));
    let low = low_support_regions(&c, 10);
    assert!(!low.is_empty(), "a 6-4 length split must be marked");
    let repeat = LEFT.len()..LEFT.len() + 40;
    for r in &low {
        assert!(
            r.start >= repeat.start && r.end <= repeat.end,
            "low support {r:?} outside the repeat {repeat:?}"
        );
    }
    let q = &c.support_quality;
    assert!(q[..LEFT.len()].iter().all(|&x| x >= 30));
    assert!(q[q.len() - RIGHT.len()..].iter().all(|&x| x >= 30));
}

#[test]
fn every_entry_point_returns_one_quality_per_base() {
    let mut reads: Vec<Vec<u8>> = (0..7).map(|_| locus(12)).collect();
    reads.push(locus(11));
    reads.push(locus(13));
    let r = slices(&reads);
    let check = |c: &poa_consensus::Consensus, what: &str| {
        assert_eq!(c.support_quality.len(), c.sequence.len(), "{what}");
    };

    check(&consensus(&r, 0, &PoaConfig::default()).unwrap(), "default");
    let plain = PoaConfig {
        homopolymer_min_run: 0,
        homopolymer_gap: None,
        ..PoaConfig::default()
    };
    check(&consensus(&r, 0, &plain).unwrap(), "plain engine");
    let mf = PoaConfig {
        consensus_mode: ConsensusMode::MajorityFrequency,
        ..PoaConfig::default()
    };
    check(&consensus(&r, 0, &mf).unwrap(), "majority frequency");
    let flanked = consensus_flanked(&r, LEFT.as_bytes(), RIGHT.as_bytes(), &plain).unwrap();
    check(&flanked, "flanked");
    assert_eq!(flanked.sequence, "CA".repeat(12).into_bytes());
    for a in consensus_multi(&r, 0, &PoaConfig::default()).unwrap() {
        check(&a, "multi");
    }
    let left: Vec<Vec<u8>> = (0..5).map(|_| locus(12)[..30].to_vec()).collect();
    let right: Vec<Vec<u8>> = (0..5).map(|_| locus(12)[34..].to_vec()).collect();
    let b = bridged_consensus(&slices(&left), 0, &slices(&right), 0, &plain).unwrap();
    check(&b, "bridged");
}

#[test]
fn homopolymer_splice_is_rescored() {
    // Length-variable runs around a single G: the homopolymer candidates can
    // replace the plain call, which must then be rescored on the raw reads.
    let mk = |a: usize, b: usize| format!("{LEFT}{}G{}{RIGHT}", "A".repeat(a), "A".repeat(b));
    let reads: Vec<Vec<u8>> = [
        (12, 11),
        (11, 12),
        (12, 11),
        (13, 10),
        (12, 11),
        (11, 12),
        (12, 11),
        (12, 12),
    ]
    .iter()
    .map(|&(a, b)| mk(a, b).into_bytes())
    .collect();
    let c = consensus(&slices(&reads), 0, &PoaConfig::default()).unwrap();
    assert_eq!(c.support_quality.len(), c.sequence.len());
}

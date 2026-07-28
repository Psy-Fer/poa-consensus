//! Regression test for the homogeneous-VNTR flank-fabrication bug
//! (Known Bug #3 scoring degeneracy; see `design/flank_fabrication_bug3.md`).
//!
//! In a homogeneous tandem repeat, matching a base to a phase-shifted repeat
//! node scores identically to matching the correct node, so a narrow banded
//! aligner can fold a read onto the wrong unit and fabricate bases present in
//! no read. `PoaConfig::path_score_bias` fixes this by adding an abPOA-style
//! log-odds read-support term to the alignment DP.
//!
//! The input is a SYNTHETIC read set (`tests/fixtures/synthetic_vntr_fold.fa`,
//! from `bench/make_synthetic_vntr.py --units 18 --depth 24 --seed 1`); no
//! private data. Fabrication is measured minimap2-free as the number of
//! consensus k-mers (k=14) absent from every read in both orientations -- a
//! fabricated stretch necessarily introduces novel k-mers.

use std::collections::HashSet;

use poa_consensus::{AlignmentMode, PoaConfig, PoaGraph, SeedSelection, auto_orient, select_seed};

fn parse_fasta(text: &str) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    let mut cur = Vec::new();
    for line in text.lines() {
        if line.starts_with('>') || line.starts_with(';') {
            if !cur.is_empty() {
                out.push(std::mem::take(&mut cur));
            }
        } else {
            cur.extend_from_slice(line.trim().as_bytes());
        }
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            x => x,
        })
        .collect()
}

fn read_kmers(reads: &[Vec<u8>], k: usize) -> HashSet<Vec<u8>> {
    let mut set = HashSet::new();
    for r in reads {
        for w in r.windows(k) {
            set.insert(w.to_vec());
        }
        for w in revcomp(r).windows(k) {
            set.insert(w.to_vec());
        }
    }
    set
}

fn consensus_novel_kmers(reads: &[Vec<u8>], path_score_bias: bool, k: usize) -> usize {
    let slices: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
    let seed = select_seed(&slices, &SeedSelection::Auto).unwrap();
    let oriented = auto_orient(reads, seed);
    let os: Vec<&[u8]> = oriented.iter().map(|r| r.as_ref()).collect();

    let cfg = PoaConfig {
        band_width: 50,
        adaptive_band: true,
        alignment_mode: AlignmentMode::SemiGlobal,
        path_score_bias,
        ..PoaConfig::default()
    };
    let mut g = PoaGraph::new(os[seed], cfg).unwrap();
    for (i, r) in os.iter().enumerate() {
        if i != seed {
            g.add_read(r).unwrap();
        }
    }
    let cons = g.consensus().unwrap().sequence;

    let rk = read_kmers(reads, k);
    cons.windows(k).filter(|w| !rk.contains(*w)).count()
}

#[test]
fn path_score_bias_removes_vntr_flank_fabrication() {
    let text = include_str!("fixtures/synthetic_vntr_fold.fa");
    let reads = parse_fasta(text);
    assert_eq!(reads.len(), 24, "fixture should have 24 synthetic reads");

    const K: usize = 14;
    let off = consensus_novel_kmers(&reads, false, K);
    let on = consensus_novel_kmers(&reads, true, K);

    // Non-vacuous: without the fix the consensus fabricates a large run of
    // k-mers present in no read (observed 92 at authoring time). The exact
    // count is not asserted -- only that fabrication is clearly present -- so
    // unrelated algorithm changes don't make this brittle.
    assert!(
        off >= 20,
        "expected substantial fabrication WITHOUT path_score_bias (novel {K}-mers), got {off}"
    );
    // The fix: path_score_bias eliminates the fabrication (observed 0).
    assert!(
        on <= 2,
        "path_score_bias should remove the VNTR flank fabrication: novel {K}-mers went {off} -> {on}"
    );
    assert!(
        on * 5 < off,
        "path_score_bias should sharply reduce fabrication: {off} -> {on}"
    );
}

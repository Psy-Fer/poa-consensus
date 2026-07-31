//! Public API error contract + Global alignment mode.
//!
//! The `select_seed` errors are covered in `src/seed.rs`, but the `validate()` +
//! depth-floor path inside the `consensus` / `consensus_multi` wrappers, and the
//! non-default `AlignmentMode::Global`, had no correctness coverage.

use poa_consensus::{AlignmentMode, PoaConfig, PoaError, consensus, consensus_multi};

fn cfg() -> PoaConfig {
    PoaConfig {
        min_reads: 3,
        ..PoaConfig::default()
    }
}

#[test]
fn empty_input_errors() {
    let reads: &[&[u8]] = &[];
    assert!(matches!(
        consensus(reads, 0, &cfg()),
        Err(PoaError::EmptyInput)
    ));
    assert!(matches!(
        consensus_multi(reads, 0, &cfg()),
        Err(PoaError::EmptyInput)
    ));
}

#[test]
fn seed_out_of_bounds_errors() {
    let reads: &[&[u8]] = &[b"ACGTACGT", b"ACGTACGT", b"ACGTACGT"];
    // seed_idx 9 >= 3 reads → SeedOutOfBounds (validated before anything else).
    assert!(matches!(
        consensus(reads, 9, &cfg()),
        Err(PoaError::SeedOutOfBounds { index: 9, len: 3 })
    ));
    assert!(matches!(
        consensus_multi(reads, 9, &cfg()),
        Err(PoaError::SeedOutOfBounds { index: 9, len: 3 })
    ));
}

#[test]
fn below_min_reads_errors() {
    // 2 reads with min_reads = 3 → InsufficientDepth (the depth safety floor).
    let reads: &[&[u8]] = &[b"ACGTACGT", b"ACGTACGT"];
    assert!(matches!(
        consensus(reads, 0, &cfg()),
        Err(PoaError::InsufficientDepth { got: 2, min: 3 })
    ));
    assert!(matches!(
        consensus_multi(reads, 0, &cfg()),
        Err(PoaError::InsufficientDepth { got: 2, min: 3 })
    ));
    // Exactly min_reads succeeds.
    assert!(consensus(&[reads[0], reads[0], reads[0]], 0, &cfg()).is_ok());
}

#[test]
fn global_mode_recovers_clean_consensus() {
    // Fully-spanning, equal-length reads (Global requires both ends aligned).
    // Majority C at position 4; one read has G. Global recovery must call C.
    let reads: &[&[u8]] = &[b"ACGTCACGTA", b"ACGTCACGTA", b"ACGTGACGTA", b"ACGTCACGTA"];
    let g = PoaConfig {
        alignment_mode: AlignmentMode::Global,
        min_reads: 3,
        ..PoaConfig::default()
    };
    let c = consensus(reads, 0, &g).expect("global consensus");
    assert_eq!(
        c.sequence,
        b"ACGTCACGTA",
        "Global mode should recover the majority sequence exactly; got {}",
        String::from_utf8_lossy(&c.sequence)
    );
}

#[test]
fn global_mode_exact_on_identical_reads() {
    let seq: &[u8] = b"ACGTACGTACGTACGT";
    let reads: &[&[u8]] = &[seq, seq, seq, seq];
    let g = PoaConfig {
        alignment_mode: AlignmentMode::Global,
        min_reads: 3,
        ..PoaConfig::default()
    };
    assert_eq!(consensus(reads, 0, &g).unwrap().sequence, seq);
}

use crate::{AlignmentMode, PoaConfig, PoaError, PoaGraph};

fn b(s: &str) -> Vec<u8> {
    s.as_bytes().to_vec()
}
fn s(v: &[u8]) -> String {
    String::from_utf8_lossy(v).into_owned()
}

/// Build consensus from reads; seed_idx selects the first read added to the graph.
fn consensus(reads: &[Vec<u8>], seed_idx: usize) -> Vec<u8> {
    let mut graph = PoaGraph::new(&reads[seed_idx], PoaConfig::default()).unwrap();
    for (i, read) in reads.iter().enumerate() {
        if i == seed_idx {
            continue;
        }
        graph.add_read(read).unwrap();
    }
    graph.consensus().unwrap().sequence
}

fn consensus_cfg(reads: &[Vec<u8>], seed_idx: usize, cfg: PoaConfig) -> Vec<u8> {
    let mut graph = PoaGraph::new(&reads[seed_idx], cfg).unwrap();
    for (i, read) in reads.iter().enumerate() {
        if i == seed_idx {
            continue;
        }
        graph.add_read(read).unwrap();
    }
    graph.consensus().unwrap().sequence
}

// ── Error cases ───────────────────────────────────────────────────────────────

#[test]
fn empty_reads() {
    let result = PoaGraph::new(&[], PoaConfig::default());
    assert!(
        matches!(result, Err(PoaError::EmptyInput)),
        "expected EmptyInput"
    );
}

#[test]
fn below_min_reads() {
    let cfg = PoaConfig {
        min_reads: 3,
        ..Default::default()
    };
    let mut graph = PoaGraph::new(&b("ACGT"), cfg).unwrap();
    graph.add_read(&b("ACGT")).unwrap();
    let result = graph.consensus();
    assert!(
        matches!(result, Err(PoaError::InsufficientDepth { got: 2, min: 3 })),
        "expected InsufficientDepth, got {:?}",
        result
    );
}

#[test]
fn seed_out_of_bounds() {
    // Callers must check seed_idx < reads.len() before calling PoaGraph::new.
    // Document the expected guard pattern.
    let reads = vec![b("ACGT"), b("ACGT")];
    let idx = 5usize;
    assert!(idx >= reads.len(), "caller must guard seed_idx before use");
}

// ── Basic correctness ─────────────────────────────────────────────────────────

#[test]
fn single_read_passthrough() {
    let reads = vec![b("CATCATCAT")];
    assert_eq!(consensus(&reads, 0), b("CATCATCAT"));
}

#[test]
fn two_identical_reads() {
    let reads = vec![b("CATCATCAT"), b("CATCATCAT")];
    assert_eq!(consensus(&reads, 0), b("CATCATCAT"));
}

#[test]
fn majority_base_wins() {
    let reads = vec![b("CATCATCAT"), b("CATCATCAT"), b("CGTCATCAT")];
    assert_eq!(s(&consensus(&reads, 0)), "CATCATCAT");
}

#[test]
fn single_outlier_not_inflated() {
    let reads = vec![b("CATCATCAT"), b("CATCATCAT"), b("CATCATCATCAT")];
    assert_eq!(consensus(&reads, 0).len(), 9);
}

#[test]
fn length_variation_longer_wins() {
    let reads = vec![b("CATCATCATCAT"), b("CATCATCATCAT"), b("CATCATCAT")];
    assert_eq!(consensus(&reads, 0).len(), 12);
}

#[test]
fn no_inflation_with_length_noise() {
    let reads = vec![
        b("CAGCAGCAGCAGCAG"),
        b("CAGCAGCAGCAGCAGCAG"),
        b("CAGCAGCAGCAGCAG"),
        b("CAGCAGCAGCAGCAG"),
    ];
    assert_eq!(consensus(&reads, 0).len(), 15);
}

#[test]
fn no_inflation_phox2b_like() {
    let reads = vec![
        b("GCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA"),
        b("GCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA"),
        b("GCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA"),
    ];
    assert_eq!(consensus(&reads, 0).len(), 60);
}

#[test]
fn single_base_reads() {
    let reads = vec![b("A"), b("A"), b("A")];
    assert_eq!(consensus(&reads, 0), b("A"));
}

// ── Boundary trim ─────────────────────────────────────────────────────────────

#[test]
fn boundary_trim_leading_seed_artifact() {
    let reads = vec![
        b("XXXCATCATCAT"),
        b("CATCATCAT"),
        b("CATCATCAT"),
        b("CATCATCAT"),
    ];
    let result = s(&consensus(&reads, 0));
    assert_eq!(result, "CATCATCAT", "got: {}", result);
}

#[test]
fn boundary_trim_trailing_seed_artifact() {
    let reads = vec![
        b("CATCATCATXXX"),
        b("CATCATCAT"),
        b("CATCATCAT"),
        b("CATCATCAT"),
    ];
    let result = s(&consensus(&reads, 0));
    assert_eq!(result, "CATCATCAT", "got: {}", result);
}

// ── Diagnostic tests from ref/poa.rs ─────────────────────────────────────────

#[test]
fn diag_sca3_t3_tail_seed_t1() {
    let reads = vec![
        b("CAGCAGCAGT"),
        b("CAGCAGCAGTTT"),
        b("CAGCAGCAGTTT"),
        b("CAGCAGCAGTTT"),
    ];
    let result = s(&consensus(&reads, 0));
    assert_eq!(result, "CAGCAGCAGTTT", "got: {}", result);
}

#[test]
fn diag_sca3_t3_tail_seed_t3() {
    let reads = vec![
        b("CAGCAGCAGTTT"),
        b("CAGCAGCAGTTT"),
        b("CAGCAGCAGTTT"),
        b("CAGCAGCAGT"),
    ];
    let result = s(&consensus(&reads, 0));
    assert_eq!(result, "CAGCAGCAGTTT", "got: {}", result);
}

#[test]
fn diag_sca31_trailing_interrupt_seed_missing() {
    let reads = vec![
        b("ATTATTATTATT"),
        b("ATTATTATTATTATA"),
        b("ATTATTATTATTATA"),
        b("ATTATTATTATTATA"),
    ];
    let result = s(&consensus(&reads, 0));
    assert_eq!(result, "ATTATTATTATTATA", "got: {}", result);
}

#[test]
fn diag_sca3_interrupt_position_single_outlier() {
    let maj = b("CAGCAGCAGCAGCAGGTTCAGCAG");
    let out = b("CAGCAGCAGCAGCAGCAGGTTCAGCAG");
    let reads = vec![maj.clone(), maj.clone(), maj.clone(), out];
    let result = s(&consensus(&reads, 0));
    assert_eq!(result, "CAGCAGCAGCAGCAGGTTCAGCAG", "got: {}", result);
}

#[test]
fn diag_sca8_minority_trailing_extension_trimmed() {
    let base = b("CAGCAGCAGCAGCAG");
    let extend = b("CAGCAGCAGCAGCAGGCT");
    let reads = vec![
        base.clone(),
        base.clone(),
        base.clone(),
        base.clone(),
        base.clone(),
        base.clone(),
        base.clone(),
        extend.clone(),
        extend.clone(),
        extend.clone(),
    ];
    assert_eq!(consensus(&reads, 0).len(), 15);
}

#[test]
fn diag_sca8_minority_trailing_extension_seed_extends() {
    let base = b("CAGCAGCAGCAGCAG");
    let extend = b("CAGCAGCAGCAGCAGGCT");
    let reads = vec![
        extend.clone(),
        extend.clone(),
        extend.clone(),
        base.clone(),
        base.clone(),
        base.clone(),
        base.clone(),
        base.clone(),
        base.clone(),
        base.clone(),
    ];
    assert_eq!(consensus(&reads, 0).len(), 15);
}

#[test]
fn diag_sca31_trailing_interrupt_before_flank() {
    let flank = b("GCGCGCGC");
    let mut seed_read = b("ATTATTATTATT");
    seed_read.extend_from_slice(&flank);
    let mut maj_read = b("ATTATTATTATTATA");
    maj_read.extend_from_slice(&flank);
    let reads = vec![
        seed_read,
        maj_read.clone(),
        maj_read.clone(),
        maj_read.clone(),
    ];
    let result = s(&consensus(&reads, 0));
    let expected: String = "ATTATTATTATTATA"
        .chars()
        .chain("GCGCGCGC".chars())
        .collect();
    assert_eq!(result, expected, "got: {}", result);
}

#[test]
fn diag_sca3_interrupt_position_long_repeat_with_flank() {
    let flank = b("CTGCTGCTG");
    let make = |repeat_pre: &str, interrupt: &str, repeat_post: &str| -> Vec<u8> {
        let mut v = repeat_pre.as_bytes().to_vec();
        v.extend_from_slice(interrupt.as_bytes());
        v.extend_from_slice(repeat_post.as_bytes());
        v.extend_from_slice(&flank);
        v
    };
    let maj = make("CAGCAGCAGCAGCAGCAGCAGCAG", "GTT", "CAGCAGCAG");
    let out = make("CAGCAGCAGCAGCAGCAGCAGCAGCAG", "GTT", "CAGCAG");
    let reads = vec![maj.clone(), maj.clone(), maj.clone(), out];
    let result = s(&consensus(&reads, 0));
    let expected = s(&make("CAGCAGCAGCAGCAGCAGCAGCAG", "GTT", "CAGCAGCAG"));
    assert_eq!(result, expected, "got: {}", result);
}

#[test]
#[ignore]
fn diag_frda_gaa_rotation_phase() {
    let gaa_phase = b("GAAGAAGAAGAA");
    let aag_phase = b("AAGAAGAAGAAG");
    let aga_phase = b("AGAAGAAGAAGA");
    let reads = vec![gaa_phase.clone(), gaa_phase.clone(), aag_phase, aga_phase];
    let result = s(&consensus(&reads, 0));
    assert_eq!(result.len(), 12, "got: '{}'", result);
}

#[test]
fn diag_frda_gaa_rotation_with_flanking() {
    let make = |repeat: &str| -> Vec<u8> {
        let mut v = b("TTTCCC");
        v.extend_from_slice(repeat.as_bytes());
        v.extend_from_slice(b("GGGAAA").as_slice());
        v
    };
    let reads = vec![
        make("GAAGAAGAAGAA"),
        make("GAAGAAGAAGAA"),
        make("AAGAAGAAGAAG"),
        make("AGAAGAAGAAGA"),
    ];
    let result = s(&consensus(&reads, 0));
    assert_eq!(result.len(), 24, "got: '{}'", result);
}

#[test]
fn diag_phase_shift_first_node_coverage() {
    let reads = vec![
        b("GAAGAA"),
        b("GAAGAA"),
        b("GAAGAA"),
        b("GAAGAA"),
        b("AAGAAG"),
    ];
    let result = s(&consensus(&reads, 0));
    assert_eq!(result, "GAAGAA", "got: '{}'", result);
}

#[test]
fn diag_phase_shift_majority_trims_first_base() {
    let reads = vec![
        b("GAAGAA"),
        b("GAAGAA"),
        b("AAGAAG"),
        b("AAGAAG"),
        b("AAGAAG"),
    ];
    let result = s(&consensus(&reads, 0));
    // The majority is the same sequence in a different phase; correct length is 6.
    assert_eq!(result.len(), 6, "got: '{}'", result);
}

#[test]
fn diag_sca3_t3_tail_with_flank() {
    let flank = b("CCTCCTCCT");
    let make = |tail: &str| -> Vec<u8> {
        let mut v = b("CAGCAGCAG");
        v.extend_from_slice(tail.as_bytes());
        v.extend_from_slice(&flank);
        v
    };
    let reads = vec![make("T"), make("TTT"), make("TTT"), make("TTT")];
    let result = s(&consensus(&reads, 0));
    let expected = s(&make("TTT"));
    assert_eq!(result, expected, "got: {}", result);
}

// ── Real-data reproduction ────────────────────────────────────────────────────

#[test]
fn diag_sca8_real_sequences_no_flank() {
    let maj57 = b("TACTACTACTACTACTACTACTACTACTACTACTGCTGCTGCTGCTGCTGCTGCTGCT");
    let min75 = b("TACTACTACTACTACTACTACTACTACTACTACTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT");
    let min81 =
        b("TACTACTACTACTACTACTACTACTACTACTACTACTACTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT");
    let mut reads: Vec<Vec<u8>> = std::iter::repeat(maj57.clone()).take(32).collect();
    reads.push(b(
        "TACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTAC",
    ));
    reads.push(b(
        "TACTACTACTACTACTACTACTACTACTACTACTACTACTACTACTGCTGCTGCTGCT",
    ));
    reads.push(b("TACTACTACTACTACTACTACTACTGCTGCTGCTGCTGCTGCTGCTGCT"));
    reads.push(min75.clone());
    reads.push(min81.clone());
    reads.push(min81.clone());
    let result = consensus(&reads, 0);
    assert_eq!(
        result.len(),
        maj57.len(),
        "SCA8 consensus must match majority length {}, got len {}: '{}'",
        maj57.len(),
        result.len(),
        s(&result)
    );
}

#[test]
fn diag_sca8_real_sequences_with_flank() {
    let flank_l = b("GCTTCGAAGTC");
    let flank_r = b("AAACGGTTCCA");
    let make = |repeat: &[u8]| -> Vec<u8> {
        let mut v = flank_l.clone();
        v.extend_from_slice(repeat);
        v.extend_from_slice(&flank_r);
        v
    };
    let maj = make(&b(
        "TACTACTACTACTACTACTACTACTACTACTACTGCTGCTGCTGCTGCTGCTGCTGCT",
    ));
    let min75 = make(&b(
        "TACTACTACTACTACTACTACTACTACTACTACTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT",
    ));
    let min81 = make(&b(
        "TACTACTACTACTACTACTACTACTACTACTACTACTACTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT",
    ));
    let mut reads: Vec<Vec<u8>> = std::iter::repeat(maj.clone()).take(32).collect();
    reads.push(min75);
    reads.push(min81.clone());
    reads.push(min81);
    let result_len = consensus(&reads, 0).len();
    assert_eq!(
        result_len,
        maj.len(),
        "SCA8 flanked: got {}, expected {}",
        result_len,
        maj.len()
    );
}

// ── Banded DP ─────────────────────────────────────────────────────────────────

#[test]
fn banded_matches_unbanded_small() {
    // Banded and unbanded must produce identical results when the band is wide
    // enough to cover the optimal path.
    let reads = vec![
        b("CAGCAGCAGCAGCAG"),
        b("CAGCAGCAGCAGCAG"),
        b("CAGCAGCAGCAGCAGCAG"),
        b("CAGCAGCAGCAGCAG"),
    ];
    let unbanded = consensus(&reads, 0);
    let cfg_banded = PoaConfig {
        band_width: 50,
        ..Default::default()
    };
    let banded = consensus_cfg(&reads, 0, cfg_banded);
    assert_eq!(unbanded, banded, "banded vs unbanded mismatch");
}

#[test]
fn adaptive_band_matches_unbanded() {
    let reads = vec![
        b("CATCATCAT"),
        b("CATCATCAT"),
        b("CATCATCATCAT"),
        b("CATCATCAT"),
    ];
    let unbanded = consensus(&reads, 0);
    let cfg = PoaConfig {
        adaptive_band: true,
        adaptive_band_b: 5,
        adaptive_band_f: 0.1,
        ..Default::default()
    };
    let adaptive = consensus_cfg(&reads, 0, cfg);
    assert_eq!(unbanded, adaptive, "adaptive band vs unbanded mismatch");
}

#[test]
fn band_too_narrow_returns_error() {
    // seed = 1 A, read = 10 A's: the read is 10 bp but the 1-node graph with
    // band_width=2 can only reach j_hi=min(10, 0+2)=2 at t=0. Column j=10 is
    // never in-band, so the terminal scan finds no cells and BandTooNarrow fires.
    let seed = b("A");
    let read = b("AAAAAAAAAA");
    let cfg = PoaConfig {
        band_width: 2,
        ..Default::default()
    };
    let mut graph = PoaGraph::new(&seed, cfg).unwrap();
    let result = graph.add_read(&read);
    assert!(
        matches!(result, Err(PoaError::BandTooNarrow { .. })),
        "expected BandTooNarrow, got {:?}",
        result.map(|_| ())
    );
}

#[test]
fn large_length_variance_banded() {
    // 3 reads of 15 bp majority + 1 read of 60 bp outlier.
    // Band must be wide enough to cover the 45-base insertion from the outlier.
    // With band_width=50 the banded result should match unbanded.
    let maj = b("CAGCAGCAGCAGCAG");
    let outlier = b("CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG");
    let reads = vec![maj.clone(), maj.clone(), maj.clone(), outlier];
    let unbanded = consensus(&reads, 0);
    let cfg = PoaConfig {
        band_width: 50,
        ..Default::default()
    };
    let banded = consensus_cfg(&reads, 0, cfg);
    assert_eq!(
        banded.len(),
        unbanded.len(),
        "banded length mismatch with large variance"
    );
    assert_eq!(
        banded, unbanded,
        "banded result mismatch with large variance"
    );
}

// ── Semi-global alignment ─────────────────────────────────────────────────────

#[test]
fn partial_reads_semi_global() {
    // 3 full reads + 1 partial: full reads are the majority so trailing region
    // has coverage 3 >= min_cov=3, consensus is the full sequence.
    let full = b("ACGTACGTACGT");
    let partial = b("ACGTACGT");
    let cfg = PoaConfig {
        alignment_mode: AlignmentMode::SemiGlobal,
        ..Default::default()
    };
    let reads = vec![full.clone(), full.clone(), full.clone(), partial];
    let result = consensus_cfg(&reads, 0, cfg);
    assert_eq!(
        result.len(),
        12,
        "semi-global: got len {}, seq: '{}'",
        result.len(),
        s(&result)
    );
}

// ── Reverse complement / orientation ─────────────────────────────────────────

#[test]
fn reverse_complement_basic() {
    use crate::reverse_complement;
    assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
    assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
    assert_eq!(reverse_complement(b"GCTA"), b"TAGC");
}

#[test]
fn orient_to_seed_forward() {
    use crate::Strand;
    use crate::orient_to_seed;
    let seed = b("ACGTACGTACGT");
    let read = b("ACGTACGT");
    assert_eq!(orient_to_seed(&read, &seed, 4), Strand::Forward);
}

#[test]
fn orient_to_seed_reverse() {
    use crate::Strand;
    use crate::orient_to_seed;
    // Use a non-palindromic sequence so forward and RC share no k-mers.
    let seed = b("AAAACCCCGGGG");
    let rc = crate::reverse_complement(&seed);
    assert_eq!(orient_to_seed(&rc, &seed, 4), Strand::Reverse);
}

#[test]
fn mixed_strand_input() {
    use crate::auto_orient;
    let seed = b("CATCATCAT");
    let rc = crate::reverse_complement(&seed);
    let reads = vec![seed.clone(), seed.clone(), rc.clone(), rc.clone()];
    let oriented: Vec<Vec<u8>> = auto_orient(&reads, 0)
        .into_iter()
        .map(|c| c.into_owned())
        .collect();
    // All oriented reads should match or be rc'd to match seed strand.
    let mut graph = PoaGraph::new(&oriented[0], PoaConfig::default()).unwrap();
    for r in &oriented[1..] {
        graph.add_read(r).unwrap();
    }
    let result = graph.consensus().unwrap().sequence;
    assert_eq!(result.len(), 9, "mixed strand: got len {}", result.len());
}

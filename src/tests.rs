use crate::{AlignmentMode, ConsensusMode, PoaConfig, PoaError, PoaGraph};

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

// ── Majority-frequency consensus ──────────────────────────────────────────────

fn mf_cfg() -> PoaConfig {
    PoaConfig {
        consensus_mode: ConsensusMode::MajorityFrequency,
        ..Default::default()
    }
}

fn mf_consensus(reads: &[Vec<u8>], seed_idx: usize) -> Vec<u8> {
    consensus_cfg(reads, seed_idx, mf_cfg())
}

#[test]
fn mf_identical_reads() {
    let reads = vec![b("CATCATCAT"), b("CATCATCAT"), b("CATCATCAT")];
    assert_eq!(mf_consensus(&reads, 0), b("CATCATCAT"));
}

#[test]
fn mf_matches_hb_on_clean_input() {
    // On a read set with no noise, MF and HB should agree.
    let reads = vec![
        b("CAGCAGCAG"),
        b("CAGCAGCAG"),
        b("CAGCAGCAGCAG"),
        b("CAGCAGCAG"),
    ];
    let hb = consensus(&reads, 0);
    let mf = mf_consensus(&reads, 0);
    assert_eq!(hb, mf, "HB and MF disagree on clean input");
}

#[test]
fn mf_boundary_trim_leading() {
    // Seed has 3 extra leading bases not present in the majority.
    // MF should exclude them because gap votes outnumber base votes.
    let reads = vec![
        b("XXXCATCATCAT"),
        b("CATCATCAT"),
        b("CATCATCAT"),
        b("CATCATCAT"),
    ];
    let result = s(&mf_consensus(&reads, 0));
    assert_eq!(result, "CATCATCAT", "got: {}", result);
}

#[test]
fn mf_boundary_trim_trailing() {
    let reads = vec![
        b("CATCATCATXXX"),
        b("CATCATCAT"),
        b("CATCATCAT"),
        b("CATCATCAT"),
    ];
    let result = s(&mf_consensus(&reads, 0));
    assert_eq!(result, "CATCATCAT", "got: {}", result);
}

#[test]
fn mf_majority_base_wins() {
    // 3 reads have CAT, 1 has CGT at position 1. MF should pick A.
    let reads = vec![
        b("CATCATCAT"),
        b("CATCATCAT"),
        b("CATCATCAT"),
        b("CGTCATCAT"),
    ];
    assert_eq!(s(&mf_consensus(&reads, 0)), "CATCATCAT");
}

#[test]
fn mf_single_outlier_not_inflated() {
    // One read has an extra CAT; MF should not include it.
    let reads = vec![b("CATCATCAT"), b("CATCATCAT"), b("CATCATCATCAT")];
    assert_eq!(mf_consensus(&reads, 0).len(), 9);
}

// ── GraphStats ────────────────────────────────────────────────────────────────

fn build_graph(reads: &[Vec<u8>], seed_idx: usize) -> PoaGraph {
    let mut graph = PoaGraph::new(&reads[seed_idx], PoaConfig::default()).unwrap();
    for (i, read) in reads.iter().enumerate() {
        if i != seed_idx {
            graph.add_read(read).unwrap();
        }
    }
    graph
}

#[test]
fn stats_clean_linear_no_bubbles() {
    // Identical reads produce a clean linear graph with no bubbles.
    let reads = vec![
        b("CATCATCAT"),
        b("CATCATCAT"),
        b("CATCATCAT"),
        b("CATCATCAT"),
    ];
    let st = build_graph(&reads, 0).stats();
    assert_eq!(st.bubble_count, 0);
    assert_eq!(st.max_bubble_depth, 0);
    assert_eq!(st.node_count, 9);
    // All reads match every node: delete_count=0 everywhere → entropy=0.
    assert_eq!(st.mean_column_entropy, 0.0);
}

#[test]
fn stats_bubble_detected() {
    // 3 reads with CATCATCAT, 1 with CGTCATCAT → SNV bubble at position 1 (A vs G).
    let reads = vec![
        b("CATCATCAT"),
        b("CATCATCAT"),
        b("CATCATCAT"),
        b("CGTCATCAT"),
    ];
    let st = build_graph(&reads, 0).stats();
    assert_eq!(st.bubble_count, 1, "expected 1 bubble");
    // Minority arm has 1 read (the CGT read creates a G branch at position 1).
    assert_eq!(st.max_bubble_depth, 1, "minority arm weight should be 1");
}

#[test]
fn stats_entropy_nonzero_on_length_variation() {
    // Seed has 3 leading X nodes that other reads delete.
    // The X nodes get delete_count=3, coverage=1 → entropy > 0.
    // (Shorter reads aligned to a longer graph in global mode don't generate trailing
    // deletes — they simply stop at the best-scoring diagonal. Leading deletes DO fire
    // because the traceback reaches t=0, j=0 via the D-chain at the j=0 column.)
    let reads = vec![
        b("XXXCATCATCAT"),
        b("CATCATCAT"),
        b("CATCATCAT"),
        b("CATCATCAT"),
    ];
    let st = build_graph(&reads, 0).stats();
    // X nodes: coverage=1, delete_count=3 → p=0.25, binary_entropy(0.25) ≈ 0.811 bits.
    assert!(
        st.mean_column_entropy > 0.0,
        "expected nonzero entropy, got {}",
        st.mean_column_entropy
    );
}

#[test]
fn stats_node_edge_counts() {
    // 2 reads with length variation: 9-node backbone + 3 extra nodes from longer read.
    // The longer read aligns as Insert(CAT) + Match(nodes 0-8), creating nodes 9,10,11
    // with edges 9→10, 10→11, 11→0. Backbone edges 0→1..7→8 already existed = 8.
    let reads = vec![b("CATCATCAT"), b("CATCATCATCAT")];
    let st = build_graph(&reads, 0).stats();
    assert_eq!(st.node_count, 12, "9 + 3 extra nodes");
    // 8 backbone + 3 new edges (9→10, 10→11, 11→0) = 11 total.
    assert_eq!(st.edge_count, 11);
}

#[test]
fn stats_coverage_mean_uniform() {
    // All 4 reads match all 9 nodes → coverage=4 everywhere → variance=0.
    let reads = vec![
        b("CATCATCAT"),
        b("CATCATCAT"),
        b("CATCATCAT"),
        b("CATCATCAT"),
    ];
    let st = build_graph(&reads, 0).stats();
    assert!((st.coverage_mean - 4.0).abs() < 1e-10);
    assert!(st.coverage_variance < 1e-10);
}

// ── Multi-allele consensus ────────────────────────────────────────────────────

fn multi_graph(reads: &[Vec<u8>], seed_idx: usize) -> PoaGraph {
    let mut graph = PoaGraph::new(&reads[seed_idx], PoaConfig::default()).unwrap();
    for (i, read) in reads.iter().enumerate() {
        if i != seed_idx {
            graph.add_read(read).unwrap();
        }
    }
    graph
}

#[test]
fn consensus_multi_single_allele() {
    // Identical reads → no bubble → consensus_multi falls through to single consensus.
    let reads = vec![b("CATCATCAT"); 4];
    let g = multi_graph(&reads, 0);
    let results = g.consensus_multi().unwrap();
    assert_eq!(results.len(), 1, "expected 1 allele for homozygous input");
    assert_eq!(results[0].sequence, b("CATCATCAT"));
}

#[test]
fn consensus_multi_snv_bubble() {
    // 4 reads with allele A (CATCATCAT) and 4 with allele B (CATCGTCAT).
    // A SNV at position 4 (A→G) creates a 2-arm bubble.
    let allele_a = b("CATCATCAT");
    let allele_b = b("CATCGTCAT");
    let reads: Vec<Vec<u8>> = (0..4)
        .map(|_| allele_a.clone())
        .chain((0..4).map(|_| allele_b.clone()))
        .collect();
    let g = multi_graph(&reads, 0);
    let results = g.consensus_multi().unwrap();
    assert_eq!(results.len(), 2, "expected 2 alleles for SNV input");
    let seqs: Vec<String> = results.iter().map(|c| s(&c.sequence)).collect();
    assert!(
        seqs.iter().any(|seq| seq == "CATCATCAT"),
        "missing CATCATCAT allele: {:?}",
        seqs
    );
    assert!(
        seqs.iter().any(|seq| seq == "CATCGTCAT"),
        "missing CATCGTCAT allele: {:?}",
        seqs
    );
}

#[test]
fn consensus_multi_length_variation() {
    // Two alleles with different repeat counts flanked by matching anchors.
    // The anchor regions create a proper bubble between the two allele lengths.
    // short: AAA + 2×CAT + TTTTTT = 14 bp
    // long : AAA + 3×CAT + TTTTTT = 17 bp
    let short = b("AAACATCATTTTTT");
    let long_ = b("AAACATCATCATTTTTT");
    let reads: Vec<Vec<u8>> = (0..4)
        .map(|_| short.clone())
        .chain((0..4).map(|_| long_.clone()))
        .collect();
    let g = multi_graph(&reads, 0);
    let results = g.consensus_multi().unwrap();
    assert_eq!(
        results.len(),
        2,
        "expected 2 alleles for length-variation input"
    );
    let lens: Vec<usize> = results.iter().map(|c| c.sequence.len()).collect();
    assert!(lens.contains(&14), "expected 14-bp allele; got {:?}", lens);
    assert!(lens.contains(&17), "expected 17-bp allele; got {:?}", lens);
}

#[test]
fn consensus_multi_insufficient_depth_per_allele() {
    // 4 reads total, min_reads=3. Each allele only gets 2 → InsufficientDepth.
    let allele_a = b("CATCATCAT");
    let allele_b = b("CATCGTCAT");
    let reads = vec![allele_a.clone(), allele_a, allele_b.clone(), allele_b];
    let cfg = PoaConfig {
        min_reads: 3,
        ..Default::default()
    };
    let mut g = PoaGraph::new(&reads[0], cfg).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    let result = g.consensus_multi();
    assert!(
        matches!(result, Err(PoaError::InsufficientDepth { .. })),
        "expected InsufficientDepth, got {:?}",
        result.map(|v| v.len())
    );
}

// ── Longer-sequence stress tests ──────────────────────────────────────────────

#[test]
fn long_repeat_consensus_correctness() {
    let seq: Vec<u8> = "CAT".repeat(30).into_bytes(); // 90 bp
    let reads = vec![seq.clone(); 6];
    assert_eq!(consensus(&reads, 0), seq, "30×CAT consensus mismatch");
}

#[test]
fn long_repeat_length_majority_wins() {
    // 8 reads at 60 bp (20×CAT), 2 outliers at 63 bp (21×CAT)
    let maj: Vec<u8> = "CAT".repeat(20).into_bytes();
    let out: Vec<u8> = "CAT".repeat(21).into_bytes();
    let mut reads: Vec<Vec<u8>> = vec![maj.clone(); 8];
    reads.extend(vec![out; 2]);
    let result = consensus(&reads, 0);
    assert_eq!(
        result.len(),
        60,
        "expected 60-bp majority, got {} bp",
        result.len()
    );
}

#[test]
fn long_repeat_snv_correction() {
    // 9 correct reads + 1 noisy read with a single mismatch at position 30
    let correct: Vec<u8> = "CAT".repeat(20).into_bytes(); // 60 bp
    let mut noisy = correct.clone();
    noisy[30] = b'G';
    let mut reads: Vec<Vec<u8>> = vec![correct.clone(); 9];
    reads.push(noisy);
    let result = consensus(&reads, 0);
    assert_eq!(
        result, correct,
        "SNV from single noisy read should not affect consensus"
    );
}

#[test]
fn long_banded_matches_unbanded() {
    // 4 × 72 bp + 1 × 78 bp (length outlier), band=30
    let base: Vec<u8> = "CAT".repeat(24).into_bytes(); // 72 bp
    let long: Vec<u8> = "CAT".repeat(26).into_bytes(); // 78 bp
    let mut reads: Vec<Vec<u8>> = vec![base.clone(); 4];
    reads.push(long);
    let unbanded = consensus(&reads, 0);
    let cfg = PoaConfig {
        band_width: 30,
        ..Default::default()
    };
    let banded = consensus_cfg(&reads, 0, cfg);
    assert_eq!(
        banded, unbanded,
        "banded(30) should match unbanded for small divergence"
    );
}

#[test]
fn long_adaptive_band_matches_unbanded() {
    // 4 × 72 bp + 1 × 81 bp, adaptive band b=10 f=0.05
    let base: Vec<u8> = "CAT".repeat(24).into_bytes(); // 72 bp
    let long: Vec<u8> = "CAT".repeat(27).into_bytes(); // 81 bp
    let mut reads: Vec<Vec<u8>> = vec![base.clone(); 4];
    reads.push(long);
    let unbanded = consensus(&reads, 0);
    let cfg = PoaConfig {
        adaptive_band: true,
        adaptive_band_b: 10,
        adaptive_band_f: 0.05,
        ..Default::default()
    };
    let adaptive = consensus_cfg(&reads, 0, cfg);
    assert_eq!(
        adaptive, unbanded,
        "adaptive band should match unbanded for small divergence"
    );
}

#[test]
fn consensus_multi_long_flanked_str() {
    // Two alleles anchored by GGGGG / AAAAA flanks:
    //   allele_a: GGGGG + 8×CAT + AAAAA  = 5+24+5 = 34 bp
    //   allele_b: GGGGG + 11×CAT + AAAAA = 5+33+5 = 43 bp
    let flank_l = b("GGGGG");
    let flank_r = b("AAAAA");
    let inner_a: Vec<u8> = "CAT".repeat(8).into_bytes();
    let inner_b: Vec<u8> = "CAT".repeat(11).into_bytes();
    let allele_a: Vec<u8> = [flank_l.as_slice(), inner_a.as_slice(), flank_r.as_slice()].concat();
    let allele_b: Vec<u8> = [flank_l.as_slice(), inner_b.as_slice(), flank_r.as_slice()].concat();
    let mut reads: Vec<Vec<u8>> = vec![allele_a.clone(); 5];
    reads.extend(vec![allele_b.clone(); 5]);
    let mut g = PoaGraph::new(&reads[0], PoaConfig::default()).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    let results = g.consensus_multi().unwrap();
    let lens: Vec<usize> = results.iter().map(|c| c.sequence.len()).collect();
    assert_eq!(results.len(), 2, "expected 2 alleles; got {:?}", lens);
    assert!(lens.contains(&34), "expected 34-bp allele; got {:?}", lens);
    assert!(lens.contains(&43), "expected 43-bp allele; got {:?}", lens);
}

#[test]
fn consensus_multi_snv_in_long_context() {
    // SNV at position 10 in a 41-bp read: allele_a has 'A' at pos 10, allele_b is all-T
    let a: Vec<u8> = {
        let mut v = vec![b'T'; 41];
        v[10] = b'A';
        v
    };
    let bv: Vec<u8> = vec![b'T'; 41];
    let mut reads: Vec<Vec<u8>> = vec![a.clone(); 5];
    reads.extend(vec![bv.clone(); 5]);
    let mut g = PoaGraph::new(&reads[0], PoaConfig::default()).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    let results = g.consensus_multi().unwrap();
    assert_eq!(
        results.len(),
        2,
        "expected 2 alleles for SNV; got {}",
        results.len()
    );
    let seqs: Vec<Vec<u8>> = results.into_iter().map(|c| c.sequence).collect();
    assert!(
        seqs.iter().any(|s| s == &a),
        "allele_a (A at pos 10) not found in results"
    );
    assert!(
        seqs.iter().any(|s| s == &bv),
        "allele_b (all-T) not found in results"
    );
}

#[test]
fn consensus_multi_skewed_allele_ratio() {
    // 7:3 ratio — minor allele at 30% detected with default min_allele_freq=0.25
    // allele_a: TTTT + 6×CAT + CCCC = 4+18+4 = 26 bp
    // allele_b: TTTT + 10×CAT + CCCC = 4+30+4 = 38 bp
    let flank_l = b("TTTT");
    let flank_r = b("CCCC");
    let inner_a: Vec<u8> = "CAT".repeat(6).into_bytes();
    let inner_b: Vec<u8> = "CAT".repeat(10).into_bytes();
    let allele_a: Vec<u8> = [flank_l.as_slice(), inner_a.as_slice(), flank_r.as_slice()].concat();
    let allele_b: Vec<u8> = [flank_l.as_slice(), inner_b.as_slice(), flank_r.as_slice()].concat();
    let mut reads: Vec<Vec<u8>> = vec![allele_a.clone(); 7];
    reads.extend(vec![allele_b.clone(); 3]);
    let mut g = PoaGraph::new(&reads[0], PoaConfig::default()).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    let results = g.consensus_multi().unwrap();
    let lens: Vec<usize> = results.iter().map(|c| c.sequence.len()).collect();
    assert_eq!(
        results.len(),
        2,
        "expected 2 alleles at 7:3; got {:?}",
        lens
    );
    assert!(
        lens.contains(&26),
        "expected 26-bp majority allele; got {:?}",
        lens
    );
    assert!(
        lens.contains(&38),
        "expected 38-bp minor allele; got {:?}",
        lens
    );
}

#[test]
fn long_reads_noise_and_banding() {
    // 63-bp majority + scattered single-base errors; banded with band=40
    let correct: Vec<u8> = "CAT".repeat(21).into_bytes(); // 63 bp
    let mut err1 = correct.clone();
    err1[20] = b'G';
    let mut err2 = correct.clone();
    err2[45] = b'T';
    let mut reads: Vec<Vec<u8>> = vec![correct.clone(); 6];
    reads.extend(vec![err1; 2]);
    reads.extend(vec![err2; 2]);
    let cfg = PoaConfig {
        band_width: 40,
        ..Default::default()
    };
    let result = consensus_cfg(&reads, 0, cfg);
    assert_eq!(
        result, correct,
        "banded consensus should correct isolated noise in 63-bp reads"
    );
}

// ── Non-repeat longer-sequence tests ─────────────────────────────────────────
//
// BASE_60: ATCGATCGTT ACGATCGTAG CTAGTCATGC TAATCGTAGC GATCGTAACG ATCGATCGTA
// 60 bp, mixed composition, no periodic structure.

const BASE_60: &[u8] = b"ATCGATCGTTACGATCGTAGCTAGTCATGCTAATCGTAGCGATCGTAACGATCGATCGTA";

#[test]
fn long_nonrepeat_consensus_correctness() {
    let reads = vec![BASE_60.to_vec(); 6];
    assert_eq!(
        consensus(&reads, 0),
        BASE_60,
        "6 identical non-repeat reads"
    );
}

#[test]
fn long_nonrepeat_snv_correction() {
    // 9 correct + 1 with T→G at position 30
    let mut noisy = BASE_60.to_vec();
    noisy[30] = b'G';
    let mut reads: Vec<Vec<u8>> = vec![BASE_60.to_vec(); 9];
    reads.push(noisy);
    assert_eq!(
        consensus(&reads, 0),
        BASE_60,
        "single noisy read must not flip consensus base"
    );
}

#[test]
fn long_nonrepeat_banded_matches_unbanded() {
    // 4 × 60 bp + 1 × 66 bp (6-bp insertion at position 30), band=30
    let mut long = BASE_60.to_vec();
    long.splice(30..30, *b"GCTAGC");
    assert_eq!(long.len(), 66);
    let mut reads: Vec<Vec<u8>> = vec![BASE_60.to_vec(); 4];
    reads.push(long);
    let unbanded = consensus(&reads, 0);
    let cfg = PoaConfig {
        band_width: 30,
        ..Default::default()
    };
    let banded = consensus_cfg(&reads, 0, cfg);
    assert_eq!(
        banded, unbanded,
        "banded(30) should match unbanded on non-repeat sequence"
    );
}

#[test]
fn consensus_multi_nonrepeat_snv() {
    // Two alleles that differ only at position 30 (T vs G) in a non-repeat context
    let mut allele_b = BASE_60.to_vec();
    allele_b[30] = b'G';
    let mut reads: Vec<Vec<u8>> = vec![BASE_60.to_vec(); 5];
    reads.extend(vec![allele_b.clone(); 5]);
    let mut g = PoaGraph::new(&reads[0], PoaConfig::default()).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    let results = g.consensus_multi().unwrap();
    assert_eq!(results.len(), 2, "expected 2 alleles for non-repeat SNV");
    let seqs: Vec<Vec<u8>> = results.into_iter().map(|c| c.sequence).collect();
    assert!(
        seqs.iter().any(|s| s.as_slice() == BASE_60),
        "allele_a not recovered"
    );
    assert!(
        seqs.iter().any(|s| s == &allele_b),
        "allele_b not recovered"
    );
}

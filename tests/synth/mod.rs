//! Deterministic synthetic read generator for the multi-allele robustness matrix.
//!
//! Zero external deps (tiny inline SplitMix64 PRNG) so tests are self-contained
//! and reproducible: the same `Spec` always yields the same reads. Every read
//! carries a `truth` allele label so tests can score read misassignment, not
//! just allele count/length. Shared across integration tests via `mod synth;`.
#![allow(dead_code)]

/// SplitMix64 — small, fast, good-enough PRNG. Seedable and portable.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        // Avoid the all-zero state producing a degenerate first draw.
        Rng(seed ^ 0x9E37_79B9_7F4A_7C15)
    }
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// Uniform in [0, 1).
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    pub fn below(&mut self, p: f64) -> bool {
        self.next_f64() < p
    }
    pub fn range(&mut self, lo: usize, hi: usize) -> usize {
        if hi <= lo {
            return lo;
        }
        lo + (self.next_u64() as usize) % (hi - lo)
    }
    pub fn base_other_than(&mut self, b: u8) -> u8 {
        let bases = [b'A', b'C', b'G', b'T'];
        loop {
            let x = bases[(self.next_u64() % 4) as usize];
            if x != b {
                return x;
            }
        }
    }
    pub fn base(&mut self) -> u8 {
        [b'A', b'C', b'G', b'T'][(self.next_u64() % 4) as usize]
    }
}

/// One allele: the mid (repeat/variant) sequence between the shared flanks, plus
/// how many reads are drawn from it (frequency = count / total).
#[derive(Clone)]
pub struct AlleleSpec {
    pub mid: Vec<u8>,
    pub n_reads: usize,
}

/// A full read-set specification. Deterministic given `seed`.
#[derive(Clone)]
pub struct Spec {
    pub left: Vec<u8>,
    pub right: Vec<u8>,
    pub alleles: Vec<AlleleSpec>,
    /// Per-base substitution / insertion / deletion error rates.
    pub sub: f64,
    pub ins: f64,
    pub del: f64,
    /// Fraction of reads truncated to a random contiguous sub-span (partial reads).
    pub partial_frac: f64,
    /// Per-read +/-1 unit-copy jitter for repeat mids (set false for exact).
    pub unit_jitter: bool,
    /// Optional unit length for jitter (bytes to add/drop at the mid boundary).
    pub unit_len: usize,
    pub seed: u64,
}

/// A generated read plus its ground-truth allele index (into `Spec::alleles`).
#[derive(Clone)]
pub struct Read {
    pub seq: Vec<u8>,
    pub truth: usize,
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            x => x,
        })
        .collect()
}

fn mutate(seq: &[u8], sub: f64, ins: f64, del: f64, rng: &mut Rng) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq {
        let r = rng.next_f64();
        if r < del {
            continue;
        } else if r < del + ins {
            out.push(rng.base());
            out.push(b);
        } else if r < del + ins + sub {
            out.push(rng.base_other_than(b));
        } else {
            out.push(b);
        }
    }
    out
}

/// Generate the read set. Reads are emitted allele-block by allele-block, each
/// read = left + mid(+/-jitter) + right, then per-base error, then optional
/// truncation. Strand is NOT randomized here (the library/CLI auto-orients; keep
/// tests focused on phasing, not orientation).
pub fn generate(spec: &Spec) -> Vec<Read> {
    let mut rng = Rng::new(spec.seed);
    let mut reads = Vec::new();
    for (ai, allele) in spec.alleles.iter().enumerate() {
        for _ in 0..allele.n_reads {
            let mut mid = allele.mid.clone();
            if spec.unit_jitter && spec.unit_len > 0 && mid.len() >= spec.unit_len {
                match rng.range(0, 5) {
                    0 => mid.truncate(mid.len() - spec.unit_len), // drop a unit
                    1 => mid.extend_from_within(mid.len() - spec.unit_len..), // add a unit
                    _ => {}
                }
            }
            let mut seq = spec.left.clone();
            seq.extend_from_slice(&mid);
            seq.extend_from_slice(&spec.right);
            seq = mutate(&seq, spec.sub, spec.ins, spec.del, &mut rng);
            if spec.partial_frac > 0.0 && rng.below(spec.partial_frac) {
                let l = seq.len();
                if l > 4 {
                    let a = rng.range(0, l / 2);
                    let b = rng.range(l / 2, l);
                    seq = seq[a..b.max(a + 1)].to_vec();
                }
            }
            reads.push(Read { seq, truth: ai });
        }
    }
    reads
}

/// Convenience: `motif` repeated `units` times.
pub fn repeat(motif: &[u8], units: usize) -> Vec<u8> {
    motif
        .iter()
        .copied()
        .cycle()
        .take(motif.len() * units)
        .collect()
}

/// A pseudo-random non-periodic sequence of length `len` (deterministic per seed).
pub fn random_seq(len: usize, seed: u64) -> Vec<u8> {
    let mut rng = Rng::new(seed);
    (0..len).map(|_| rng.base()).collect()
}

// Shared flanks used across the matrix (unique, GC-balanced, not repeat-like).
pub const LEFT: &[u8] = b"ACGTACGTACGTAGGCTTACAGGCTAGATCCAGT";
pub const RIGHT: &[u8] = b"TTGGATCCAGTAGCTTGGCAGATCGATTCAGGCA";

impl Read {
    pub fn seqs(reads: &[Read]) -> Vec<&[u8]> {
        reads.iter().map(|r| r.seq.as_slice()).collect()
    }
}

// ── Scoring ───────────────────────────────────────────────────────────────────

/// Outcome of one generated read-set vs the multi-allele output. Measures the
/// three things that matter separately, so the matrix shows *which* degrades:
/// did we call the right number of alleles, did reads land in the right group,
/// and is each allele's length right.
#[derive(Clone, Debug)]
pub struct Outcome {
    pub expected_alleles: usize,
    pub got_alleles: usize,
    pub count_correct: bool,
    /// Fraction of phased reads whose truth label != their group's assigned label
    /// (optimal 1-1 group↔label assignment when counts match; greedy otherwise).
    pub misassign_rate: f64,
    /// Per truth-allele: |best-matching consensus length − true full length|.
    /// Empty entry (usize::MAX) if no output group mapped to that truth label.
    pub len_abs_err: Vec<usize>,
}

fn permutations(n: usize) -> Vec<Vec<usize>> {
    if n == 0 {
        return vec![vec![]];
    }
    let mut out = Vec::new();
    for i in 0..n {
        for mut p in permutations(n - 1) {
            // insert i, shifting >= i up (build a permutation of 0..n)
            for v in p.iter_mut() {
                if *v >= i {
                    *v += 1;
                }
            }
            p.insert(0, i);
            out.push(p);
        }
    }
    out
}

/// Score multi-allele output against ground truth. `groups[g]` = truth-label
/// histogram of the reads assigned to output allele g (indexed via read_indices).
pub fn score(reads: &[Read], spec: &Spec, group_read_indices: &[Vec<usize>]) -> Outcome {
    let k_truth = spec.alleles.len();
    let g = group_read_indices.len();

    // Per-group truth histogram.
    let hist: Vec<Vec<usize>> = group_read_indices
        .iter()
        .map(|idxs| {
            let mut h = vec![0usize; k_truth];
            for &ri in idxs {
                if ri < reads.len() {
                    h[reads[ri].truth] += 1;
                }
            }
            h
        })
        .collect();

    let total_phased: usize = hist.iter().flatten().sum();

    // Best group→label assignment (maximize correctly-labeled reads).
    let (correct, assign): (usize, Vec<usize>) = if g == k_truth && g <= 6 {
        let mut best = (0usize, (0..g).collect::<Vec<_>>());
        for perm in permutations(g) {
            let c: usize = (0..g).map(|gi| hist[gi][perm[gi]]).sum();
            if c > best.0 {
                best = (c, perm.clone());
            }
        }
        best
    } else {
        // Greedy: each group takes its majority label (labels may repeat).
        let assign: Vec<usize> = hist
            .iter()
            .map(|h| {
                h.iter()
                    .enumerate()
                    .max_by_key(|&(_, &c)| c)
                    .map(|(i, _)| i)
                    .unwrap_or(0)
            })
            .collect();
        let correct: usize = (0..g).map(|gi| hist[gi][assign[gi]]).sum();
        (correct, assign)
    };

    let misassign_rate = if total_phased == 0 {
        1.0
    } else {
        (total_phased - correct) as f64 / total_phased as f64
    };

    // Length is scored by `score_with_lens` (needs the per-group consensus
    // lengths); this entry point leaves it unset.
    let _ = assign;
    Outcome {
        expected_alleles: k_truth,
        got_alleles: g,
        count_correct: g == k_truth,
        misassign_rate,
        len_abs_err: vec![usize::MAX; k_truth],
    }
}

/// Like [`score`] but also fills `len_abs_err` from the per-group consensus
/// lengths (same order as `group_read_indices`).
pub fn score_with_lens(
    reads: &[Read],
    spec: &Spec,
    group_read_indices: &[Vec<usize>],
    group_lens: &[usize],
) -> Outcome {
    let k_truth = spec.alleles.len();
    let g = group_read_indices.len();
    let hist: Vec<Vec<usize>> = group_read_indices
        .iter()
        .map(|idxs| {
            let mut h = vec![0usize; k_truth];
            for &ri in idxs {
                if ri < reads.len() {
                    h[reads[ri].truth] += 1;
                }
            }
            h
        })
        .collect();
    let assign: Vec<usize> = if g == k_truth && g <= 6 {
        let mut best = (0usize, (0..g).collect::<Vec<_>>());
        for perm in permutations(g) {
            let c: usize = (0..g).map(|gi| hist[gi][perm[gi]]).sum();
            if c > best.0 {
                best = (c, perm.clone());
            }
        }
        best.1
    } else {
        hist.iter()
            .map(|h| {
                h.iter()
                    .enumerate()
                    .max_by_key(|&(_, &c)| c)
                    .map(|(i, _)| i)
                    .unwrap_or(0)
            })
            .collect()
    };
    let mut out = score(reads, spec, group_read_indices);
    let full_len = |ai: usize| spec.left.len() + spec.alleles[ai].mid.len() + spec.right.len();
    let mut len_abs_err = vec![usize::MAX; k_truth];
    for (gi, &label) in assign.iter().enumerate() {
        if gi < group_lens.len() {
            let e = group_lens[gi].abs_diff(full_len(label));
            if len_abs_err[label] == usize::MAX {
                len_abs_err[label] = e;
            } else {
                len_abs_err[label] = len_abs_err[label].min(e);
            }
        }
    }
    out.len_abs_err = len_abs_err;
    out
}

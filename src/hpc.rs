//! Homopolymer-compressed (HPC) consensus.
//!
//! Each read is rewritten in an extended alphabet before POA: a run of the
//! same base of length `>= min_run` becomes a single "poly" symbol (e.g.
//! `polyA`), shorter runs stay as plain bases. With `min_run = 1` this is pure
//! homopolymer compression (every run is one symbol); with `min_run = 2` it is
//! the 8-letter `A,C,G,T,polyA,polyC,polyG,polyT` vocabulary.
//!
//! The point: in raw space, a single base between two long, length-variable
//! runs (`A{12}GA{11}` with reads carrying `A{11}GA{12}`, ...) has no fixed
//! position, and with harsh affine gaps the DP prefers two G/A mismatches over
//! a deletion + insertion, so node fusion splits the one biological base across
//! two MSA columns. In compressed space every such read reads `polyA G polyA`,
//! so the base has exactly one place to go.
//!
//! Run lengths are recovered after consensus: every read is aligned to the
//! compressed consensus, and each consensus poly symbol takes the median run
//! length of the reads that matched it with the same symbol.

use crate::config::{ConsensusMode, PoaConfig};
use crate::poa2::{Op, Poa};
use crate::types::Consensus;

/// High bit marks a poly symbol; the low bits keep the base.
const POLY: u8 = 0x80;

fn is_poly(s: u8) -> bool {
    s & POLY != 0
}

fn base_of(s: u8) -> u8 {
    s & !POLY
}

/// Compress `read`: returns the symbol string and, per symbol, the raw run
/// length it stands for (1 for plain bases).
pub(crate) fn encode(read: &[u8], min_run: usize) -> (Vec<u8>, Vec<u32>) {
    let min_run = min_run.max(1);
    let mut syms = Vec::with_capacity(read.len());
    let mut lens = Vec::with_capacity(read.len());
    let mut i = 0;
    while i < read.len() {
        let b = read[i];
        let mut j = i + 1;
        while j < read.len() && read[j] == b {
            j += 1;
        }
        let run = j - i;
        if run >= min_run {
            syms.push(b | POLY);
            lens.push(run as u32);
        } else {
            for _ in 0..run {
                syms.push(b);
                lens.push(1);
            }
        }
        i = j;
    }
    (syms, lens)
}

/// Median of a non-empty slice (lower median for even counts, so a 50/50 tie
/// between n and n+1 calls the shorter run).
fn median(v: &mut [u32]) -> u32 {
    v.sort_unstable();
    v[(v.len() - 1) / 2]
}

/// Per consensus symbol, the raw run length to expand it to. Plain symbols are
/// always 1; a poly symbol takes the median run length over reads whose own
/// symbol at that position is the same poly symbol.
///
/// Observations are "clean" when the read's neighbouring symbols also match the
/// consensus neighbours. A substitution inside a run splits it (`A5 C A6`), and
/// either fragment would otherwise vote a too-short length, so clean
/// observations are preferred; all matching observations are the fallback.
fn estimate_run_lengths(
    cons: &[u8],
    enc_reads: &[(Vec<u8>, Vec<u32>)],
    cfg: &PoaConfig,
) -> Vec<u32> {
    let g = Poa::new(cons, cfg.clone());
    let (topo, rank) = g.topo_order();
    let n = cons.len();
    let mut clean: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut all: Vec<Vec<u32>> = vec![Vec::new(); n];
    for (syms, lens) in enc_reads {
        // matched[rpos] = consensus position the read symbol matched exactly.
        let mut matched: Vec<Option<usize>> = vec![None; syms.len()];
        let mut rpos = 0usize;
        for op in g.align(syms, &topo, &rank) {
            match op {
                Op::Match(node) => {
                    // The seed chain's nodes are 0..cons.len(), one per symbol.
                    if node < n && syms[rpos] == cons[node] {
                        matched[rpos] = Some(node);
                    }
                    rpos += 1;
                }
                Op::Ins(_) => rpos += 1,
                Op::Del(_) => {}
            }
        }
        for (r, m) in matched.iter().enumerate() {
            let Some(node) = *m else { continue };
            if !is_poly(cons[node]) {
                continue;
            }
            all[node].push(lens[r]);
            // A neighbour that doesn't exist on either side (read or consensus
            // end) can't disagree.
            let left_ok = node == 0 || r == 0 || matched[r - 1] == Some(node - 1);
            let right_ok = node + 1 == n || r + 1 == syms.len() || matched[r + 1] == Some(node + 1);
            if left_ok && right_ok {
                clean[node].push(lens[r]);
            }
        }
    }
    (0..n)
        .map(|i| {
            if !is_poly(cons[i]) {
                1
            } else if !clean[i].is_empty() {
                median(&mut clean[i])
            } else if !all[i].is_empty() {
                median(&mut all[i])
            } else {
                // No read matched it exactly (should not happen for a consensus
                // symbol); fall back to the poly floor.
                1
            }
        })
        .collect()
}

/// Repeat each per-position value across its run.
fn repeat_runs<T: Copy>(v: &[T], runs: &[u32]) -> Vec<T> {
    v.iter()
        .zip(runs)
        .flat_map(|(&x, &r)| std::iter::repeat_n(x, r as usize))
        .collect()
}

/// Expand a compressed `Consensus` back to raw bases using per-symbol run
/// lengths. Per-position fields are repeated across each run; positions in
/// gaps and bubble sites are remapped. Bubble arm sequences are decoded with a
/// run length of 1 per poly symbol (their own run lengths are not estimated).
fn expand(mut c: Consensus, runs: &[u32], raw_reads: &[&[u8]]) -> Consensus {
    let n = c.sequence.len();
    debug_assert_eq!(runs.len(), n);
    // starts[i] = raw position of compressed position i; starts[n] = raw length.
    let mut starts = Vec::with_capacity(n + 1);
    let mut acc = 0usize;
    for &r in runs {
        starts.push(acc);
        acc += r as usize;
    }
    starts.push(acc);
    let seq: Vec<u8> = c
        .sequence
        .iter()
        .zip(runs)
        .flat_map(|(&s, &r)| std::iter::repeat_n(base_of(s), r as usize))
        .collect();
    c.coverage = repeat_runs(&c.coverage, runs);
    c.path_weights = repeat_runs(&c.path_weights, runs);
    c.sequence = seq;
    let map = |p: usize| starts[p.min(n)];
    for gap in &mut c.gaps {
        gap.start = map(gap.start);
        gap.end = map(gap.end);
    }
    for site in &mut c.bubble_sites {
        site.consensus_pos = map(site.consensus_pos);
        for arm in &mut site.arm_sequences {
            for s in arm.iter_mut() {
                *s = base_of(*s);
            }
        }
    }
    let mut lens: Vec<usize> = raw_reads.iter().map(|r| r.len()).collect();
    lens.sort_unstable();
    c.graph_stats.median_input_read_len = lens[lens.len() / 2];
    c
}

/// Single-allele consensus in compressed space, expanded back to raw bases.
/// Mirrors `crate::consensus`'s engine path (median-length backbone, then
/// `consensus_mode`), with the compression applied around it.
pub(crate) fn consensus_compressed(reads: &[&[u8]], config: &PoaConfig) -> Consensus {
    let min_run = config.homopolymer_min_run;
    let enc: Vec<(Vec<u8>, Vec<u32>)> = reads.iter().map(|r| encode(r, min_run)).collect();
    let syms: Vec<&[u8]> = enc.iter().map(|(s, _)| s.as_slice()).collect();
    let mut order: Vec<usize> = (0..syms.len()).collect();
    order.sort_by_key(|&i| syms[i].len());
    let med = order[order.len() / 2];
    let mut g = Poa::new(syms[med], config.clone());
    for (i, s) in syms.iter().enumerate() {
        if i != med {
            g.add_read(s);
        }
    }
    let c = match config.consensus_mode {
        ConsensusMode::MajorityFrequency => g.consensus_full_majority(),
        ConsensusMode::BestFit => g.consensus_full_best_fit(&syms),
    };
    let runs = estimate_run_lengths(&c.sequence, &enc, config);
    expand(c, &runs, reads)
}

/// Raw-space misfit of `candidate` to `reads`: mismatches + inserted +
/// deleted bases per read base (lower is better). Unlike the engine's
/// indel-only `fit_score`, mismatches count: a G shifted within a homopolymer
/// aligns as two G/A mismatches under harsh affine gaps, not as indels.
fn raw_misfit(candidate: &[u8], reads: &[&[u8]], cfg: &PoaConfig) -> f64 {
    if candidate.is_empty() {
        return f64::MAX;
    }
    let g = Poa::new(candidate, cfg.clone());
    let (topo, rank) = g.topo_order();
    let (mut cost, mut bases) = (0usize, 0usize);
    for r in reads {
        let mut rpos = 0usize;
        for op in g.align(r, &topo, &rank) {
            match op {
                Op::Match(node) => {
                    if node < candidate.len() && candidate[node] != r[rpos] {
                        cost += 1;
                    }
                    rpos += 1;
                }
                Op::Ins(_) => {
                    cost += 1;
                    rpos += 1;
                }
                Op::Del(_) => cost += 1,
            }
        }
        bases += r.len();
    }
    cost as f64 / bases.max(1) as f64
}

/// Length of the identical-base run each position belongs to.
fn run_len_at(seq: &[u8]) -> Vec<usize> {
    let mut out = vec![0; seq.len()];
    let mut i = 0;
    while i < seq.len() {
        let mut j = i + 1;
        while j < seq.len() && seq[j] == seq[i] {
            j += 1;
        }
        out[i..j].fill(j - i);
        i = j;
    }
    out
}

/// Splice the two consensuses: align `hpc` to `raw`, and for each block where
/// they differ take the compressed call only if the block touches a run of
/// `>= min_run` bases in either sequence; elsewhere keep the raw call.
/// Compression only claims to help inside homopolymers, and outside them the
/// raw engine is the better-tested answer (e.g. a periodic repeat's unit
/// count). Per-position fields come from whichever side each base came from;
/// gaps and bubble sites are the raw call's, remapped.
fn splice(raw: &Consensus, hpc: &Consensus, min_run: usize, config: &PoaConfig) -> Consensus {
    let (r, h) = (&raw.sequence, &hpc.sequence);
    let min_run = min_run.max(1);
    let (rr, hr) = (run_len_at(r), run_len_at(h));
    let touches = |runs: &[usize], a: usize, b: usize| {
        let lo = a.saturating_sub(1);
        let hi = (b + 1).min(runs.len());
        runs[lo..hi.max(lo)].iter().any(|&x| x >= min_run)
    };
    let cfg = PoaConfig {
        alignment_mode: crate::config::AlignmentMode::Global,
        ..config.clone()
    };
    let g = Poa::new(r, cfg);
    let (topo, rank) = g.topo_order();
    let ops = g.align(h, &topo, &rank);

    // (from_hpc, index) per emitted base.
    let mut src: Vec<(bool, usize)> = Vec::with_capacity(r.len().max(h.len()));
    let (mut i, mut j) = (0usize, 0usize);
    let mut block: Option<(usize, usize)> = None; // (i1, j1) of an open diff block
    let flush =
        |block: &mut Option<(usize, usize)>, i2: usize, j2: usize, src: &mut Vec<(bool, usize)>| {
            if let Some((i1, j1)) = block.take() {
                if touches(&rr, i1, i2) || touches(&hr, j1, j2) {
                    src.extend((j1..j2).map(|x| (true, x)));
                } else {
                    src.extend((i1..i2).map(|x| (false, x)));
                }
            }
        };
    for op in ops {
        match op {
            Op::Match(node) if node == i && h[j] == r[i] => {
                flush(&mut block, i, j, &mut src);
                src.push((false, i));
                i += 1;
                j += 1;
            }
            Op::Match(_) => {
                block.get_or_insert((i, j));
                i += 1;
                j += 1;
            }
            Op::Ins(_) => {
                block.get_or_insert((i, j));
                j += 1;
            }
            Op::Del(_) => {
                block.get_or_insert((i, j));
                i += 1;
            }
        }
    }
    flush(&mut block, i, j, &mut src);

    let mut out = raw.clone();
    out.sequence = src
        .iter()
        .map(|&(fh, x)| if fh { h[x] } else { r[x] })
        .collect();
    out.coverage = src
        .iter()
        .map(|&(fh, x)| if fh { hpc.coverage[x] } else { raw.coverage[x] })
        .collect();
    out.path_weights = src
        .iter()
        .map(|&(fh, x)| {
            if fh {
                hpc.path_weights[x]
            } else {
                raw.path_weights[x]
            }
        })
        .collect();
    // raw index -> first merged index at or after it.
    let mut map = vec![out.sequence.len(); r.len() + 1];
    for (k, &(fh, x)) in src.iter().enumerate().rev() {
        if !fh {
            map[x] = k;
        }
    }
    for x in (0..r.len()).rev() {
        map[x] = map[x].min(map[x + 1]);
    }
    for gap in &mut out.gaps {
        gap.start = map[gap.start.min(r.len())];
        gap.end = map[gap.end.min(r.len())];
    }
    for site in &mut out.bubble_sites {
        site.consensus_pos = map[site.consensus_pos.min(r.len())];
    }
    out
}

/// Homopolymer-aware candidates combined with the plain engine call `raw`
/// (built with `base`: both homopolymer options off). Candidates: the
/// compressed consensus when `homopolymer_min_run > 0`, and the engine run with
/// homopolymer gap costs when `homopolymer_gap` is set (gap costs are not
/// applied inside the compressed build: compressed reads have no runs). Each
/// candidate is spliced into `raw` only within homopolymer runs (`splice`), and
/// the best-fitting splice is kept only if it fits the raw reads strictly
/// better than `raw` (judged under `base` scoring).
///
/// Compression wins where run-length noise scatters an interrupting base; the
/// fit check rejects it where substitutions inside runs pile up into fabricated
/// compressed structure, and the splice keeps both options from changing calls
/// outside homopolymers (e.g. a CAG repeat's unit count).
pub(crate) fn consensus_homopolymer(
    reads: &[&[u8]],
    config: &PoaConfig,
    base: &PoaConfig,
    raw: Consensus,
) -> Consensus {
    let k = config.homopolymer_min_run;
    let mut cands: Vec<(Consensus, usize)> = Vec::new();
    if k > 0 {
        let cfg = PoaConfig {
            homopolymer_gap: None,
            ..config.clone()
        };
        cands.push((consensus_compressed(reads, &cfg), k));
    }
    if config.homopolymer_gap.is_some() {
        let cfg = PoaConfig {
            homopolymer_min_run: 0,
            ..config.clone()
        };
        cands.push((
            crate::engine_consensus(reads, &cfg),
            if k > 0 { k } else { 2 },
        ));
    }
    let fr = raw_misfit(&raw.sequence, reads, base);
    let mut best: Option<(f64, Consensus)> = None;
    for (cand, min_run) in cands {
        if cand.sequence == raw.sequence {
            continue;
        }
        let merged = splice(&raw, &cand, min_run, base);
        if merged.sequence == raw.sequence {
            continue;
        }
        let fm = raw_misfit(&merged.sequence, reads, base);
        if fm < fr && best.as_ref().is_none_or(|(b, _)| fm < *b) {
            best = Some((fm, merged));
        }
    }
    best.map_or(raw, |(_, c)| c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_min_run_1_compresses_every_run() {
        let (s, l) = encode(b"AAAGCCT", 1);
        assert_eq!(s, vec![b'A' | POLY, b'G' | POLY, b'C' | POLY, b'T' | POLY]);
        assert_eq!(l, vec![3, 1, 2, 1]);
    }

    #[test]
    fn encode_min_run_3_keeps_short_runs_raw() {
        let (s, l) = encode(b"AAAGCCT", 3);
        assert_eq!(s, vec![b'A' | POLY, b'G', b'C', b'C', b'T']);
        assert_eq!(l, vec![3, 1, 1, 1, 1]);
    }

    #[test]
    fn interrupt_between_variable_runs_lands_once() {
        // Same compressed string for every read regardless of run lengths.
        let a = encode(b"CTCAAAAAAAAAAAAGAAAAAAAAAAAGTT", 2).0;
        let b = encode(b"CTCAAAAAAAAAAAGAAAAAAAAAAAAGTT", 2).0;
        assert_eq!(a, b);
    }

    #[test]
    fn roundtrip_recovers_median_run_lengths() {
        let mk = |x: usize, y: usize| {
            let mut r = b"TTCGCATGCC".to_vec();
            r.extend(std::iter::repeat_n(b'A', x));
            r.push(b'G');
            r.extend(std::iter::repeat_n(b'A', y));
            r.extend_from_slice(b"CGTCCATGCT");
            r
        };
        let reads = [
            mk(12, 11),
            mk(11, 12),
            mk(12, 11),
            mk(13, 11),
            mk(12, 10),
            mk(12, 11),
        ];
        let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
        for min_run in [1, 2, 4] {
            let cfg = PoaConfig {
                homopolymer_min_run: min_run,
                ..PoaConfig::default()
            };
            let c = consensus_compressed(&refs, &cfg);
            assert_eq!(c.sequence, mk(12, 11), "min_run={min_run}");
            assert_eq!(c.coverage.len(), c.sequence.len());
        }
    }
}

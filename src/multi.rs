//! Multi-allele consensus pipeline (extracted from `poa2.rs` — F22).
//!
//! Structural-bubble discovery PROPOSES splits, linkage discovery refines each,
//! and `confirm_and_merge` rejects over-splits. This layer operates ON the
//! `poa2::Poa` engine (it builds a per-partition sub-graph per candidate allele);
//! the engine itself carries no multi-allele logic. Not part of the public API —
//! reached only via the crate-root `consensus_multi` / `bridged_consensus`.

use crate::config::PoaConfig;
use crate::poa2::{MatchedInc, MatchedOut, Op, Poa};
use std::collections::HashMap;

const MIN_LENGTH_GAP_BP: f64 = 8.0;
const LENGTH_SEPARATION_MADS: f64 = 3.0;
const MIN_SPREAD_FLOOR_BP: f64 = 3.0;

fn median_usize(vals: &[usize]) -> f64 {
    if vals.is_empty() {
        return 0.0;
    }
    let mut v: Vec<f64> = vals.iter().map(|&x| x as f64).collect();
    v.sort_by(|a, b| a.total_cmp(b)); // total_cmp: NaN-safe, no unwrap panic path
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

fn mad_usize(vals: &[usize], median: f64) -> f64 {
    let dev: Vec<usize> = vals
        .iter()
        .map(|&x| (x as f64 - median).abs().round() as usize)
        .collect();
    median_usize(&dev)
}

/// True when two groups' read-length distributions are genuinely bimodally
/// separated (real length-differing alleles) rather than ordinary scatter
/// around one shared mode.
fn length_separated(lens_a: &[usize], lens_b: &[usize]) -> bool {
    if lens_a.is_empty() || lens_b.is_empty() {
        return false;
    }
    let med_a = median_usize(lens_a);
    let med_b = median_usize(lens_b);
    let gap = (med_a - med_b).abs();
    if gap < MIN_LENGTH_GAP_BP {
        return false;
    }
    let spread_a = mad_usize(lens_a, med_a).max(MIN_SPREAD_FLOOR_BP);
    let spread_b = mad_usize(lens_b, med_b).max(MIN_SPREAD_FLOOR_BP);
    gap >= LENGTH_SEPARATION_MADS * spread_a.max(spread_b)
}

/// Number of `node`'s in-edges (matched-view) whose own source-side weight
/// clears `sig_threshold` (a genuine reconvergence vs one real path plus noise).
fn real_in_edge_count(
    out: &MatchedOut,
    inc: &MatchedInc,
    node: usize,
    sig_threshold: i32,
) -> usize {
    inc[node]
        .iter()
        .filter(|&&p| {
            out[p]
                .iter()
                .find(|&&(to, _)| to == node)
                .is_some_and(|&(_, w)| w >= sig_threshold)
        })
        .count()
}

/// Noise-tolerant arm-span walk over the matched-view: a fork or reconvergence
/// ends the arm only when it is a *real* one (2+ edges each clearing
/// `sig_threshold`). See the legacy `materialize_arm_len_tolerant` for why safe.
fn materialize_arm_len_tolerant(
    out: &MatchedOut,
    inc: &MatchedInc,
    start: usize,
    max_depth: usize,
    sig_threshold: i32,
) -> usize {
    let mut len = 0;
    let mut cur = start;
    for _ in 0..max_depth {
        if len >= 1 && real_in_edge_count(out, inc, cur, sig_threshold) > 1 {
            break;
        }
        len += 1;
        match out[cur].as_slice() {
            [] => break,
            [(next, _)] => cur = *next,
            edges => {
                let real_arms = edges.iter().filter(|&&(_, w)| w >= sig_threshold).count();
                if real_arms > 1 {
                    break; // a genuine nested fork -- this arm ends here.
                }
                let &(next, _) = edges.iter().max_by_key(|&&(_, w)| w).unwrap();
                cur = next;
            }
        }
    }
    len
}

/// Bubbles where at least one arm has structural size (span ≥
/// `phasing_bubble_min_span`): genuine length variants / SVs, not SNPs.
fn find_structural_bubbles(
    out: &MatchedOut,
    inc: &MatchedInc,
    topo: &[usize],
    n_reads: usize,
    cfg: &PoaConfig,
) -> Vec<(usize, Vec<usize>)> {
    let threshold = ((n_reads as f64 * cfg.min_allele_freq).ceil() as i32).max(1);
    let max_check = cfg.phasing_bubble_min_span.saturating_add(1);

    topo.iter()
        .copied()
        .filter_map(|entry_node| {
            let arms: Vec<usize> = out[entry_node]
                .iter()
                .filter(|&&(_, w)| w >= threshold)
                .map(|&(to, _)| to)
                .collect();
            if arms.len() < 2 {
                return None;
            }
            let max_span = arms
                .iter()
                .map(|&start| {
                    if real_in_edge_count(out, inc, start, threshold) > 1 {
                        0
                    } else {
                        materialize_arm_len_tolerant(out, inc, start, max_check, threshold)
                    }
                })
                .max()
                .unwrap_or(0);
            if max_span >= cfg.phasing_bubble_min_span {
                Some((entry_node, arms))
            } else {
                None
            }
        })
        .collect()
}

/// Bubbles with 2+ arms each clearing the `min_allele_freq` vote threshold
/// (any span; the SNP-fallback finder) over the matched-view.
fn find_bubbles(
    out: &MatchedOut,
    topo: &[usize],
    n_reads: usize,
    min_allele_freq: f64,
) -> Vec<(usize, Vec<usize>)> {
    let threshold = ((n_reads as f64 * min_allele_freq).ceil() as i32).max(1);
    topo.iter()
        .copied()
        .filter_map(|node_idx| {
            let arms: Vec<usize> = out[node_idx]
                .iter()
                .filter(|&&(_, w)| w >= threshold)
                .map(|&(to, _)| to)
                .collect();
            if arms.len() >= 2 {
                Some((node_idx, arms))
            } else {
                None
            }
        })
        .collect()
}

/// Groups reads by arm-choice compatibility across all structural bubbles.
/// Incremental, most-informative-first clustering (NOT plain union-find); bridge
/// candidates and unassigned reads resolved by length with an 85% plausibility
/// gate. `lens` is internal-read-indexed.
fn phasing_groups(
    edge_reads: &HashMap<(usize, usize), Vec<u32>>,
    bubbles: &[(usize, Vec<usize>)],
    n_reads: usize,
    min_reads: usize,
    lens: &[usize],
) -> Vec<Vec<usize>> {
    let n_bubbles = bubbles.len();

    let mut sig: Vec<Vec<Option<usize>>> = vec![vec![None; n_bubbles]; n_reads];
    for (b, (entry, arm_starts)) in bubbles.iter().enumerate() {
        for (arm_idx, &arm_start) in arm_starts.iter().enumerate() {
            if let Some(rs) = edge_reads.get(&(*entry, arm_start)) {
                for &r in rs {
                    let r = r as usize;
                    if r < n_reads {
                        sig[r][b] = Some(arm_idx);
                    }
                }
            }
        }
    }

    let mut assigned: Vec<usize> = Vec::new();
    let mut unassigned: Vec<usize> = Vec::new();
    for (r, row) in sig.iter().enumerate() {
        if row.iter().any(|s| s.is_some()) {
            assigned.push(r);
        } else {
            unassigned.push(r);
        }
    }

    let mut order = assigned.clone();
    order.sort_by_key(|&r| {
        let known = sig[r].iter().filter(|s| s.is_some()).count();
        (std::cmp::Reverse(known), r)
    });

    let mut clusters: Vec<Vec<Option<usize>>> = Vec::new();
    let mut cluster_members: Vec<Vec<usize>> = Vec::new();
    let mut bridge_candidates: Vec<(usize, Vec<usize>)> = Vec::new();

    for r in order {
        let row = &sig[r];
        let compatible_clusters: Vec<usize> = clusters
            .iter()
            .enumerate()
            .filter(|(_, known)| {
                (0..n_bubbles).all(|b| match (row[b], known[b]) {
                    (Some(a), Some(bv)) => a == bv,
                    _ => true,
                })
            })
            .map(|(ci, _)| ci)
            .collect();

        match compatible_clusters.as_slice() {
            [] => {
                clusters.push(row.clone());
                cluster_members.push(vec![r]);
            }
            [only] => {
                cluster_members[*only].push(r);
                for b in 0..n_bubbles {
                    if clusters[*only][b].is_none() {
                        clusters[*only][b] = row[b];
                    }
                }
            }
            _ => bridge_candidates.push((r, compatible_clusters)),
        }
    }

    const PLAUSIBLE_LEN_FRACTION: f64 = 0.85;

    if !cluster_members.is_empty() {
        let snapshot_lens: Vec<f64> = cluster_members
            .iter()
            .map(|members| median_usize(&members.iter().map(|&r| lens[r]).collect::<Vec<_>>()))
            .collect();
        let snapshot_sizes: Vec<usize> = cluster_members.iter().map(|m| m.len()).collect();
        let min_full_len = (0..cluster_members.len())
            .filter(|&c| snapshot_sizes[c] >= min_reads)
            .map(|c| snapshot_lens[c])
            .fold(f64::INFINITY, f64::min);

        let nearest_cluster_among = |read_idx: usize, candidates: &[usize]| -> usize {
            let len = lens[read_idx] as f64;
            candidates
                .iter()
                .copied()
                .min_by(|&a, &b| {
                    let da = (snapshot_lens[a] - len).abs();
                    let db = (snapshot_lens[b] - len).abs();
                    da.total_cmp(&db)
                        .then_with(|| snapshot_sizes[b].cmp(&snapshot_sizes[a]))
                        .then_with(|| a.cmp(&b))
                })
                .unwrap()
        };
        let largest_overall = (0..cluster_members.len())
            .max_by_key(|&c| snapshot_sizes[c])
            .unwrap();
        let all_clusters: Vec<usize> = (0..cluster_members.len()).collect();
        let largest_compatible = |candidates: &[usize]| -> usize {
            *candidates
                .iter()
                .max_by_key(|&&c| snapshot_sizes[c])
                .unwrap()
        };

        let mut resolved: Vec<(usize, usize)> = Vec::new();
        for (r, compat) in &bridge_candidates {
            let len = lens[*r] as f64;
            let target = if min_full_len.is_finite() && len >= PLAUSIBLE_LEN_FRACTION * min_full_len
            {
                nearest_cluster_among(*r, &all_clusters)
            } else {
                largest_compatible(compat)
            };
            resolved.push((*r, target));
        }
        for &r in &unassigned {
            let len = lens[r] as f64;
            let target = if min_full_len.is_finite() && len >= PLAUSIBLE_LEN_FRACTION * min_full_len
            {
                nearest_cluster_among(r, &all_clusters)
            } else {
                largest_overall
            };
            resolved.push((r, target));
        }
        for (r, target) in resolved {
            cluster_members[target].push(r);
        }
    } else {
        cluster_members.push(Vec::new());
        cluster_members[0].extend(bridge_candidates.into_iter().map(|(r, _)| r));
        cluster_members[0].extend(unassigned);
    }

    let mut groups: Vec<Vec<usize>> = cluster_members;
    groups.sort_unstable_by_key(|g| std::cmp::Reverse(g.len()));
    if groups.is_empty() {
        groups.push(Vec::new());
    }
    let mut i = 1;
    while i < groups.len() {
        if groups[i].len() < min_reads {
            let g = groups.remove(i);
            groups[0].extend(g);
        } else {
            i += 1;
        }
    }
    groups.retain(|g| !g.is_empty());
    groups
}

/// Safety net over `phasing_groups`: a provisional cluster is trusted as its
/// own allele only if it is length-separated from AND has a clean distinguishing
/// bubble vs every confirmed group (both gated on 2+ structural bubbles).
/// Rejected candidates route per-member to the length-nearest confirmed group.
fn validate_and_merge_groups(
    groups: Vec<Vec<usize>>,
    lens: &[usize],
    min_reads: usize,
    bubbles: &[(usize, Vec<usize>)],
    edge_reads: &HashMap<(usize, usize), Vec<u32>>,
) -> Vec<Vec<usize>> {
    if groups.len() < 2 {
        return groups;
    }
    let mut groups = groups;
    groups.sort_unstable_by_key(|g| std::cmp::Reverse(g.len()));

    let n_bubbles = bubbles.len();
    let require_length_separation = n_bubbles >= 2;

    const CLEAR_MAJORITY: f64 = 0.60;
    const MIN_ARM_COV: usize = 2;

    let arm_of: Vec<HashMap<usize, usize>> = bubbles
        .iter()
        .map(|(entry, arm_starts)| {
            let mut m: HashMap<usize, usize> = HashMap::new();
            for (k, &a) in arm_starts.iter().enumerate() {
                if let Some(rs) = edge_reads.get(&(*entry, a)) {
                    for &r in rs {
                        m.entry(r as usize).or_insert(k);
                    }
                }
            }
            m
        })
        .collect();

    let group_bubble_stats = |grp: &[usize], b: usize| -> (Option<usize>, f64, usize) {
        let mut counts: HashMap<usize, usize> = HashMap::new();
        for &r in grp {
            if let Some(&k) = arm_of[b].get(&r) {
                *counts.entry(k).or_default() += 1;
            }
        }
        let cov: usize = counts.values().sum();
        match counts.into_iter().max_by_key(|&(_, c)| c) {
            Some((arm, c)) => (Some(arm), c as f64 / cov as f64, cov),
            None => (None, 0.0, 0),
        }
    };

    let clean_distinguishing = |ga: &[usize], gb: &[usize]| -> bool {
        (0..n_bubbles).any(|b| {
            let (ma, fa, ca) = group_bubble_stats(ga, b);
            let (mb, fb, cb) = group_bubble_stats(gb, b);
            ma.is_some()
                && mb.is_some()
                && ma != mb
                && fa >= CLEAR_MAJORITY
                && fb >= CLEAR_MAJORITY
                && ca >= MIN_ARM_COV
                && cb >= MIN_ARM_COV
        })
    };

    let mut confirmed: Vec<Vec<usize>> = vec![groups[0].clone()];
    let mut confirmed_lens: Vec<Vec<usize>> = vec![groups[0].iter().map(|&r| lens[r]).collect()];

    for cand in groups.into_iter().skip(1) {
        let cand_lens: Vec<usize> = cand.iter().map(|&r| lens[r]).collect();

        let significant = cand.len() >= min_reads;
        let separated_from_all = !require_length_separation
            || confirmed_lens
                .iter()
                .all(|existing| length_separated(&cand_lens, existing));

        if significant && separated_from_all {
            let structurally_distinct = !require_length_separation
                || confirmed.iter().all(|c| clean_distinguishing(&cand, c));
            if structurally_distinct {
                confirmed.push(cand);
                confirmed_lens.push(cand_lens);
            } else {
                let cand_med = median_usize(&cand_lens);
                let target = (0..confirmed.len())
                    .filter(|&c| !clean_distinguishing(&cand, &confirmed[c]))
                    .min_by(|&a, &b| {
                        let da = (median_usize(&confirmed_lens[a]) - cand_med).abs();
                        let db = (median_usize(&confirmed_lens[b]) - cand_med).abs();
                        da.total_cmp(&db)
                    })
                    .expect("structurally_distinct false => at least one indistinguishable group");
                confirmed_lens[target].extend(cand_lens);
                confirmed[target].extend(cand);
            }
        } else if confirmed.len() == 1 {
            confirmed_lens[0].extend(cand_lens);
            confirmed[0].extend(cand);
        } else {
            let snapshot_medians: Vec<f64> =
                confirmed_lens.iter().map(|l| median_usize(l)).collect();
            let mut per_read_targets: Vec<(usize, usize)> = Vec::with_capacity(cand.len());
            for (&r, &len) in cand.iter().zip(cand_lens.iter()) {
                let len = len as f64;
                let target = snapshot_medians
                    .iter()
                    .enumerate()
                    .min_by(|&(_, &a), &(_, &b)| {
                        let da = (a - len).abs();
                        let db = (b - len).abs();
                        da.total_cmp(&db)
                    })
                    .map(|(i, _)| i)
                    .unwrap_or(0);
                per_read_targets.push((r, target));
            }
            for (r, target) in per_read_targets {
                confirmed_lens[target].push(lens[r]);
                confirmed[target].push(r);
            }
        }
    }

    confirmed
}

/// SNP-fallback partition: assign each read to the first of the top-2 arms
/// (by read count) it appears in; unassigned reads go to the largest group.
fn partition_reads_by_bubble(
    edge_reads: &HashMap<(usize, usize), Vec<u32>>,
    entry_node: usize,
    arm_starts: &[usize],
    n_reads: usize,
) -> Vec<Vec<usize>> {
    let mut arm_order: Vec<usize> = (0..arm_starts.len()).collect();
    arm_order.sort_unstable_by(|&a, &b| {
        let wa = edge_reads
            .get(&(entry_node, arm_starts[a]))
            .map_or(0, |v| v.len());
        let wb = edge_reads
            .get(&(entry_node, arm_starts[b]))
            .map_or(0, |v| v.len());
        wb.cmp(&wa)
    });
    let n_arms = arm_order.len().min(2);

    let mut groups: Vec<Vec<usize>> = vec![Vec::new(); n_arms];
    let mut assigned = vec![false; n_reads];

    for (slot, &arm_idx) in arm_order[..n_arms].iter().enumerate() {
        let arm_start = arm_starts[arm_idx];
        if let Some(rs) = edge_reads.get(&(entry_node, arm_start)) {
            for &r in rs {
                let r = r as usize;
                if r < n_reads && !assigned[r] {
                    groups[slot].push(r);
                    assigned[r] = true;
                }
            }
        }
    }

    let largest = groups
        .iter()
        .enumerate()
        .max_by_key(|(_, g)| g.len())
        .map(|(i, _)| i)
        .unwrap_or(0);
    for (r, &done) in assigned.iter().enumerate() {
        if !done {
            groups[largest].push(r);
        }
    }

    groups.retain(|g| !g.is_empty());
    groups
}

/// Build a poa2 graph from `reads`, seeding on the median-length (non-empty)
/// read. Returns `(graph, ext)` where `ext[internal] = external index into
/// reads`; `graph.n_reads == ext.len()` (empty reads are skipped consistently).
fn build_median(reads: &[&[u8]], cfg: &PoaConfig) -> (Poa, Vec<usize>) {
    let mut order: Vec<usize> = (0..reads.len()).filter(|&i| !reads[i].is_empty()).collect();
    order.sort_by_key(|&i| reads[i].len());
    let med = order[order.len() / 2];
    let mut g = Poa::new(reads[med], cfg.clone());
    let mut ext = vec![med];
    for (i, r) in reads.iter().enumerate() {
        if i != med && !r.is_empty() {
            g.add_read(r);
            ext.push(i);
        }
    }
    (g, ext)
}

/// Classic O(nm) edit distance (small consensus sequences, exact).
fn edit_distance(a: &[u8], b: &[u8]) -> usize {
    let (n, m) = (a.len(), b.len());
    let mut prev: Vec<usize> = (0..=m).collect();
    let mut cur = vec![0usize; m + 1];
    for i in 1..=n {
        cur[0] = i;
        for j in 1..=m {
            let cost = usize::from(a[i - 1] != b[j - 1]);
            cur[j] = (prev[j] + 1).min(cur[j - 1] + 1).min(prev[j - 1] + cost);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev[m]
}

/// Linkage-based multi-allele consensus (design/multi_allele_phasing_design.md).
///
/// Phases reads by **linkage** — a genuine allele is a set of differences that
/// co-occur on the same reads across many het sites — rather than by any single
/// bubble's arm span (a folded-graph artifact). The dominant read bipartition is
/// only *trusted* when its two sides yield genuinely different consensus
/// sequences (the consensus-difference confirmation), so phase-shift noise in a
/// homogeneous repeat — which splits reads but yields the SAME consensus on both
/// sides — collapses back to one allele. Recurses for 3+ alleles. Each allele's
/// `read_indices` are the external (input) read indices.
pub fn linkage_consensus_multi(
    reads: &[&[u8]],
    config: &PoaConfig,
) -> Result<Vec<crate::types::Consensus>, crate::error::PoaError> {
    use crate::error::PoaError;
    // Enforce the depth floor on the count of reads that carry sequence. Empty
    // reads contribute nothing, and an all-empty input would otherwise panic the
    // median-seed selection in `build_median` (empty `order`, out-of-bounds
    // median index) rather than returning a clean error.
    let nonempty = reads.iter().filter(|r| !r.is_empty()).count();
    if nonempty < config.min_reads {
        return Err(PoaError::InsufficientDepth {
            got: nonempty,
            min: config.min_reads,
        });
    }
    let all: Vec<usize> = (0..reads.len()).filter(|&i| !reads[i].is_empty()).collect();
    let mut groups = Vec::new();
    linkage_partition_into(reads, &all, config, 0, &mut groups);
    groups_to_consensuses(reads, &groups, config)
}

/// For each read, align it to a fresh graph seeded on `candidate` and return
/// its `(n_insert, n_delete)` op counts. Public hook so the analysis layer
/// (`analysis::consensus_fit`) can score a candidate against reads without
/// depending on any engine internals.
pub fn align_indel_counts(
    candidate: &[u8],
    reads: &[&[u8]],
    cfg: &PoaConfig,
) -> Vec<(usize, usize)> {
    if candidate.is_empty() {
        return Vec::new();
    }
    let g = Poa::new(candidate, cfg.clone());
    let (topo, rank) = g.topo_order();
    reads
        .iter()
        .map(|r| {
            let (mut ins, mut del) = (0usize, 0usize);
            for op in g.align(r, &topo, &rank) {
                match op {
                    Op::Ins(_) => ins += 1,
                    Op::Del(_) => del += 1,
                    Op::Match(_) => {}
                }
            }
            (ins, del)
        })
        .collect()
}

/// Max recursion depth (allele count is bounded by 2^depth; 3 covers up to ~8).
const LINKAGE_MAX_DEPTH: usize = 3;
/// A read bipartition is only trusted as a real allele split when it is this
/// clean (fraction of calls matching their side's per-site majority). Noise
/// sub-clusters of one allele split at lower consistency, so this gate stops the
/// recursion from over-splitting a single allele at higher error rates.
/// Read-support floor for a het-site branch to count as significant: the
/// `min_allele_freq` fraction of reads (min 2). Deliberately NOT maxed with
/// `min_reads` — that is the consensus-*depth* floor (checked separately per
/// group), and using it here would hide a legitimate minority arm below depth
/// (e.g. a 3-read minor allele when min_reads=4), causing the split to be missed.
/// Shared minority-arm support floor: the minimum read count for an arm to be a
/// callable allele at frequency `min_allele_freq`. Used BOTH by the multi-allele
/// splitter (so the split threshold) AND by the analysis layer's competing-allele
/// recommendation (`analysis::has_competing_allele`), so "re-run with --multi" can
/// never fire at a looser threshold than the engine will actually act on.
pub(crate) fn phasing_floor(n_reads: usize, min_allele_freq: f64) -> u32 {
    ((n_reads as f64 * min_allele_freq).ceil() as u32).max(2)
}

/// Minimum median read-length gap (bp) between two groups to accept a LENGTH
/// split. Above a broad allele's own stutter-bisection spacing (a few bp) but
/// well below any genuine repeat-count difference (tens of bp). Rejects the
/// stutter over-split without needing to know the repeat unit size.
const LINKAGE_MIN_LEN_GAP: f64 = 12.0;

/// Consistency bar for confirming a SAME-LENGTH (substitution) split,
/// where length gives no signal: a clean haplotype split (even one SNP) is
/// ~1.0 consistent, error sub-clusters of one allele ~0.7, so this admits the
/// former and rejects the latter.
const LINKAGE_SUBST_CONSISTENCY: f64 = 0.85;

/// Build a consensus for `subset` of `reads`, tagged with the subset's external
/// indices (sorted).
fn subset_consensus(
    reads: &[&[u8]],
    subset: &[usize],
    config: &PoaConfig,
) -> crate::types::Consensus {
    let sub_reads: Vec<&[u8]> = subset.iter().map(|&i| reads[i]).collect();
    let (g, _) = build_median(&sub_reads, config);
    let mut c = g.consensus_full();
    let mut idx = subset.to_vec();
    idx.sort_unstable();
    c.read_indices = idx;
    c
}

/// Recursively split `subset` by the dominant linkage bipartition, keeping a
/// split only when the two sides' consensuses genuinely differ.
/// Recursively partition `subset` by the dominant linkage bipartition, keeping a
/// split only when the two sides are genuinely distinct alleles (`groups_distinct`
/// — length-bimodal, or a clean same-length substitution). Returns read-index
/// groups (recurses for 3+ alleles). This is the linkage *discovery* primitive,
/// used both standalone and to refine the hybrid's structural proposal.
fn linkage_partition_into(
    reads: &[&[u8]],
    subset: &[usize],
    config: &PoaConfig,
    depth: usize,
    out: &mut Vec<Vec<usize>>,
) {
    if subset.len() < 2 * config.min_reads || depth >= LINKAGE_MAX_DEPTH {
        out.push(subset.to_vec());
        return;
    }
    let sub_reads: Vec<&[u8]> = subset.iter().map(|&i| reads[i]).collect();
    let (g, ext) = build_median(&sub_reads, config);
    let bp = g
        .phasing_matrix(phasing_floor(g.n_reads, config.min_allele_freq))
        .dominant_bipartition();

    let mut s0 = Vec::new();
    let mut s1 = Vec::new();
    for internal in 0..g.n_reads {
        let orig = subset[ext[internal]];
        if bp.assign[internal] == 0 {
            s0.push(orig);
        } else {
            s1.push(orig);
        }
    }
    // No `n_informative >= 2` pre-filter here: a clean LENGTH variant can diverge
    // at a single fork, and `groups_distinct` confirms it by length-bimodality
    // regardless of site count. The substitution branch of `groups_distinct` keeps
    // its own `informative >= 2` requirement, so a noise split still can't pass.
    if s0.len() < config.min_reads
        || s1.len() < config.min_reads
        || !groups_distinct(reads, &s0, &s1, config)
    {
        out.push(subset.to_vec());
        return;
    }
    linkage_partition_into(reads, &s0, config, depth + 1, out);
    linkage_partition_into(reads, &s1, config, depth + 1, out);
}

/// Multi-allele consensus over `reads` (WORK IN PROGRESS — not the production
/// path; `crate::consensus_multi` routes through the legacy engine). Structural
/// -bubble phasing (cross-bubble clustering + validate/merge), an unbanded
/// rebuild when no structural bubble is found on a banded graph, then a
/// single-best SNP-bubble fallback; each partition is rebuilt into its own poa2
/// sub-graph for the per-allele consensus.
///
/// KNOWN LIMITATION (does not yet reach legacy real-data parity): poa2's
/// node-fusion folds periodic repeats, so a shorter allele's length difference
/// is smeared across many equal-length phase bubbles instead of surfacing as one
/// high-span structural bubble the way legacy's separate `bypass_edges`
/// representation exposed it. Empirically, no phasing signal tried matches legacy
/// on the real HiFi multi scenarios: structural-only under-splits the clean
/// synthetic 3-allele case and fails multi_cag20_50 / multi_skew_cag20_40 /
/// multi_gaa30_100; all-bubbles over-splits multi_gaa30_100 (best at 19/20);
/// length-refinement over-splits worse (real reads have within-allele length
/// variance). Reaching parity needs flanking-anchor pre-processing
/// (`extract_repeat_segment`), a separate effort. This code carries the ported
/// clustering logic (correct given the right bubble set) for that follow-up.
pub fn consensus_multi(
    reads: &[&[u8]],
    config: &PoaConfig,
) -> Result<Vec<crate::types::Consensus>, crate::error::PoaError> {
    use crate::error::PoaError;
    // Enforce the depth floor on the count of reads that carry sequence. Empty
    // reads contribute nothing, and an all-empty input would otherwise panic the
    // median-seed selection in `build_median` (empty `order`, out-of-bounds
    // median index) rather than returning a clean error.
    let nonempty = reads.iter().filter(|r| !r.is_empty()).count();
    if nonempty < config.min_reads {
        return Err(PoaError::InsufficientDepth {
            got: nonempty,
            min: config.min_reads,
        });
    }
    consensus_multi_impl(reads, config, true)
}

fn consensus_multi_impl(
    reads: &[&[u8]],
    config: &PoaConfig,
    allow_unbanded_retry: bool,
) -> Result<Vec<crate::types::Consensus>, crate::error::PoaError> {
    let groups = structural_read_groups(reads, config, allow_unbanded_retry);
    groups_to_consensuses(reads, &groups, config)
}

/// Structural-phasing read partitions (external read indices): structural
/// bubbles → cross-bubble clustering → validate/merge; an unbanded rebuild when
/// no structural bubble is found on a banded graph; then a single-best
/// SNP-bubble fallback. A single group (all reads) means "one allele". This is
/// the split PROPOSAL the hybrid engine confirms with the linkage gate.
fn structural_read_groups(
    reads: &[&[u8]],
    config: &PoaConfig,
    allow_unbanded_retry: bool,
) -> Vec<Vec<usize>> {
    let (g, ext) = build_median(reads, config);
    let lens: Vec<usize> = ext.iter().map(|&e| reads[e].len()).collect();
    let (topo, _) = g.topo_order();
    let (mout, min) = g.matched_view();

    let structural = find_structural_bubbles(&mout, &min, &topo, g.n_reads, config);
    if structural.is_empty()
        && allow_unbanded_retry
        && (config.band_width > 0 || config.adaptive_band)
    {
        let mut cfg2 = config.clone();
        cfg2.band_width = 0;
        cfg2.adaptive_band = false;
        return structural_read_groups(reads, &cfg2, false);
    }

    let edge_reads = g.edge_reads();
    let internal: Vec<Vec<usize>> = if !structural.is_empty() {
        let g0 = phasing_groups(&edge_reads, &structural, g.n_reads, config.min_reads, &lens);
        validate_and_merge_groups(g0, &lens, config.min_reads, &structural, &edge_reads)
    } else {
        let bubbles = find_bubbles(&mout, &topo, g.n_reads, config.min_allele_freq);
        if bubbles.is_empty() {
            return vec![ext.clone()];
        }
        let (entry, arms) = bubbles
            .iter()
            .max_by_key(|(entry, arms)| {
                arms.iter()
                    .map(|&a| edge_reads.get(&(*entry, a)).map_or(0, |v| v.len()))
                    .min()
                    .unwrap_or(0)
            })
            .unwrap();
        partition_reads_by_bubble(&edge_reads, *entry, arms, g.n_reads)
    };
    // Map internal read indices back to external (input) indices.
    internal
        .into_iter()
        .map(|grp| grp.into_iter().map(|i| ext[i]).collect())
        .collect()
}

/// Build one consensus per read group (external indices), depth-guarded, with a
/// single-allele shortcut when there is fewer than two groups.
fn groups_to_consensuses(
    reads: &[&[u8]],
    groups: &[Vec<usize>],
    config: &PoaConfig,
) -> Result<Vec<crate::types::Consensus>, crate::error::PoaError> {
    use crate::error::PoaError;
    if groups.len() < 2 {
        let all: Vec<usize> = match groups.first() {
            Some(g) => g.clone(),
            None => (0..reads.len()).filter(|&i| !reads[i].is_empty()).collect(),
        };
        return Ok(vec![subset_consensus(reads, &all, config)]);
    }
    let mut out = Vec::with_capacity(groups.len());
    for group in groups {
        if group.len() < config.min_reads {
            return Err(PoaError::InsufficientDepth {
                got: group.len(),
                min: config.min_reads,
            });
        }
        out.push(subset_consensus(reads, group, config));
    }
    Ok(out)
}

/// Are two read groups genuinely distinct alleles (vs one allele a structural
/// proposer over-split)? The linkage confirmation gate: a LENGTH variant is
/// length-bimodal (median gap >= 3×MAD, robust to stutter); a SUBSTITUTION
/// haplotype is same-length but the `gi`-vs-`gj` split is linkage-consistent
/// (reads cleanly take sides at ≥2 het sites) with a real consensus difference.
fn groups_distinct(reads: &[&[u8]], gi: &[usize], gj: &[usize], config: &PoaConfig) -> bool {
    // Frequency floor: neither side can be a *called* allele if it holds fewer
    // reads than `min_allele_freq` of the whole locus. This rejects the classic
    // multi-allele over-split — a small stutter/partial cluster whose median
    // length lands in the VALLEY between two real alleles (so it is length-
    // separated from both) and would otherwise pass the length branch as a
    // spurious third allele. `reads.len()` is the full locus read count (gi/gj
    // index into it), so the floor is a genuine global allele frequency. A real
    // minority allele sits at or above `min_allele_freq` by definition; a valley
    // noise cluster does not. (Lower `min_allele_freq` to call true mosaics.)
    let freq_floor = phasing_floor(reads.len(), config.min_allele_freq) as usize;
    if gi.len().min(gj.len()) < freq_floor {
        return false;
    }
    let lens_i: Vec<usize> = gi.iter().map(|&r| reads[r].len()).collect();
    let lens_j: Vec<usize> = gj.iter().map(|&r| reads[r].len()).collect();
    // Length variant: separated (median gap vs MAD) AND the median gap clears an
    // absolute floor. The floor rejects a broad stutter spread bisected into
    // adjacent sub-clusters (a few bp apart) while keeping real alleles, which
    // differ by many repeat units. (`length_separated` alone is fooled by the
    // bisection; a pure outlier-gap test is fooled the other way — real alleles'
    // stutter/partial tails fill the valley on noisy data.)
    let median_gap = (median_usize(&lens_i) - median_usize(&lens_j)).abs();
    if length_separated(&lens_i, &lens_j) && median_gap >= LINKAGE_MIN_LEN_GAP {
        return true;
    }
    // Same length: confirm the split against the linkage signal on a pooled graph.
    let pooled_reads: Vec<&[u8]> = gi.iter().chain(gj.iter()).map(|&r| reads[r]).collect();
    let (g, ext) = build_median(&pooled_reads, config);
    let floor = phasing_floor(g.n_reads, config.min_allele_freq);
    let m = g.phasing_matrix(floor);
    // side per internal read: pooled input index < gi.len() ⟹ from gi (side 0).
    let side: Vec<i8> = (0..g.n_reads)
        .map(|internal| i8::from(ext[internal] >= gi.len()))
        .collect();
    let bp = m.score_partition(side);
    let ci = subset_consensus(reads, gi, config);
    let cj = subset_consensus(reads, gj, config);
    let same_length = ci.sequence.len().abs_diff(cj.sequence.len()) <= 2;
    let edit = edit_distance(&ci.sequence, &cj.sequence);
    // SUBSTITUTION split (same length; length gave no signal). The discriminator
    // is CONSISTENCY: a real haplotype split — even a single clean SNP —
    // separates reads cleanly across its het site(s) (consistency ~1.0), while
    // two error/stutter sub-clusters of one allele do not (~0.7 at ONT error). A
    // high consistency bar with ≥1 informative site and ≥1 consensus difference
    // admits the clean single-SNP diploid (`two_allele_snp_recovery`) while
    // rejecting noise. Periodic-repeat phase noise — which can fake one spurious
    // single-column fork with consistency ~1.0 — is handled *upstream* by the
    // frequency floor at the top of this function (such a phantom sub-cluster is
    // a small minority, below `min_allele_freq`), not by this branch. The
    // `same_length` gate is load-bearing: a length-jitter sub-split (which
    // `length_separated` already rejected as sub-bimodal) must NOT leak in here
    // as a "substitution".
    same_length
        && bp.consistency >= LINKAGE_SUBST_CONSISTENCY
        && bp.n_informative_sites >= 1
        && edit >= 1
}

/// Merge any pair of proposed groups that is NOT a genuinely distinct allele,
/// to fixpoint. This is the hybrid's safety gate: it rejects the single-allele
/// over-splits that structural phasing produces (two sub-clusters of one allele
/// that are neither length-bimodal nor linkage-consistent).
fn confirm_and_merge(
    reads: &[&[u8]],
    mut groups: Vec<Vec<usize>>,
    config: &PoaConfig,
) -> Vec<Vec<usize>> {
    loop {
        let mut merge: Option<(usize, usize)> = None;
        'outer: for i in 0..groups.len() {
            for j in (i + 1)..groups.len() {
                if !groups_distinct(reads, &groups[i], &groups[j], config) {
                    merge = Some((i, j));
                    break 'outer;
                }
            }
        }
        match merge {
            Some((i, j)) => {
                let gj = groups.remove(j);
                groups[i].extend(gj);
            }
            None => break,
        }
    }
    groups
}

/// HYBRID multi-allele consensus: structural-bubble discovery *proposes* the
/// split (good sensitivity on clean length variants), the linkage
/// consensus-difference + bimodality gate *confirms* it (rejects the
/// single-allele over-splits that structural phasing alone produces). Combines
/// the split-sensitivity of `consensus_multi` with the safety of
/// `linkage_consensus_multi`.
pub fn hybrid_consensus_multi(
    reads: &[&[u8]],
    config: &PoaConfig,
) -> Result<Vec<crate::types::Consensus>, crate::error::PoaError> {
    use crate::error::PoaError;
    // Enforce the depth floor on the count of reads that carry sequence. Empty
    // reads contribute nothing, and an all-empty input would otherwise panic the
    // median-seed selection in `build_median` (empty `order`, out-of-bounds
    // median index) rather than returning a clean error.
    let nonempty = reads.iter().filter(|r| !r.is_empty()).count();
    if nonempty < config.min_reads {
        return Err(PoaError::InsufficientDepth {
            got: nonempty,
            min: config.min_reads,
        });
    }
    // Structural phasing PROPOSES the split (good sensitivity on clean length
    // variants). Then refine each proposed group with linkage discovery — this
    // finds splits structural phasing MISSES entirely (e.g. a folded periodic
    // GAA diploid where no structural bubble surfaces, or a 3rd allele inside a
    // structural group). Finally confirm_and_merge rejects any over-splits.
    let proposed = structural_read_groups(reads, config, true);
    let mut refined = Vec::new();
    for group in &proposed {
        linkage_partition_into(reads, group, config, 0, &mut refined);
    }
    let confirmed = confirm_and_merge(reads, refined, config);
    groups_to_consensuses(reads, &confirmed, config)
}

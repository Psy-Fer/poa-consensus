//! Linkage-based multi-allele phasing signal.
//!
//! See `design/multi_allele_phasing_design.md`. The core idea: a real allele is
//! a set of differences that **co-occur on the same reads across multiple
//! sites** (linkage). Sequencing error is independent per read per site (no
//! linkage). So instead of phasing on any single bubble's arm span — a graph
//! representation artifact that collapses on folded periodic repeats — we build
//! a read × het-site matrix and measure cross-site linkage.
//!
//! This module is the *interpretation* layer: it operates on a `PhasingMatrix`
//! (extracted from a built graph by the engine) and scores linkage, so it is
//! independent of the graph representation. Clustering (MEC / correlation) and
//! the single-mode het/mosaic/error classifier build on top of these scores.

/// A read's call at a site is a branch index, or `MISSING` when the read does
/// not cover the site (partial reads, semi-global ends, or a read that took a
/// below-threshold minor branch).
pub const MISSING: i16 = -1;

/// One het site: a graph fork where reads take ≥2 well-supported branches.
/// Branch `k` is identified by the successor node id `branches[k]`.
#[derive(Clone, Debug)]
pub struct Site {
    /// The fork node (the position at which reads diverge).
    pub node: usize,
    /// Successor node ids of the significant branches (each ≥ the support floor).
    pub branches: Vec<usize>,
    /// Reads on each branch (parallel to `branches`).
    pub branch_support: Vec<u32>,
}

/// The dominant 2-way read bipartition and how well it explains the sites.
#[derive(Clone, Debug)]
pub struct Bipartition {
    /// Side (0 or 1) per read.
    pub assign: Vec<i8>,
    /// Fraction of covered calls matching their side's per-site majority branch
    /// (1.0 = a perfectly clean two-allele split; low = no consistent partition).
    pub consistency: f64,
    /// Smaller side as a fraction of assigned reads (allele-frequency estimate;
    /// ~0.5 balanced het, small = mosaic/minority, ~0 = effectively one allele).
    pub minority_frac: f64,
    /// Sites where the two sides' majority branches actually differ (the sites
    /// that carry the split).
    pub n_informative_sites: usize,
}

/// Read × het-site allele-call matrix — the "fragment matrix" of read-backed
/// phasing / Minimum Error Correction. `calls[r][s]` is read `r`'s branch index
/// at site `s`, or [`MISSING`].
#[derive(Clone, Debug)]
pub struct PhasingMatrix {
    pub n_reads: usize,
    pub sites: Vec<Site>,
    pub calls: Vec<Vec<i16>>,
}

impl PhasingMatrix {
    pub fn n_sites(&self) -> usize {
        self.sites.len()
    }

    /// Reads that carry a non-missing call at site `s`.
    pub fn covered(&self, s: usize) -> usize {
        (0..self.n_reads)
            .filter(|&r| self.calls[r][s] != MISSING)
            .count()
    }

    /// Linkage between two sites in [0, 1]: over reads that cover *both*, how
    /// well the branch at `a` predicts the branch at `b` (a perfect haplotype
    /// pairing → 1; independent/random → the no-information baseline → 0).
    ///
    /// Computed as the best-matching agreement of the branch×branch contingency
    /// table, normalized against the random baseline so that a real allele split
    /// (co-segregating reads) scores high regardless of branch balance, while a
    /// noise fork (random membership) scores near 0.
    pub fn linkage(&self, a: usize, b: usize) -> f64 {
        let (na, nb) = (self.sites[a].branches.len(), self.sites[b].branches.len());
        if na == 0 || nb == 0 {
            return 0.0;
        }
        // Contingency table counts[branch_a][branch_b].
        let mut counts = vec![vec![0u32; nb]; na];
        let mut co = 0u32;
        for r in 0..self.n_reads {
            let (ca, cb) = (self.calls[r][a], self.calls[r][b]);
            if ca != MISSING && cb != MISSING {
                counts[ca as usize][cb as usize] += 1;
                co += 1;
            }
        }
        if co == 0 {
            return 0.0;
        }
        // Predict b from a: each a-branch commits to its most-common b-branch.
        let pred: u32 = counts
            .iter()
            .map(|row| row.iter().copied().max().unwrap_or(0))
            .sum();
        // Baseline: guessing b's single most-common branch overall.
        let mut col_tot = vec![0u32; nb];
        for row in &counts {
            for (j, &c) in row.iter().enumerate() {
                col_tot[j] += c;
            }
        }
        let baseline = col_tot.iter().copied().max().unwrap_or(0);
        if co <= baseline {
            return 0.0;
        }
        // How far above baseline the branch-conditioned prediction gets.
        ((pred.saturating_sub(baseline)) as f64 / (co - baseline) as f64).clamp(0.0, 1.0)
    }

    /// The strongest linkage between any pair of sites (0 if fewer than 2 sites).
    /// A locus with a genuine second allele has at least one strongly-linked
    /// pair; a single noisy allele does not.
    pub fn max_linkage(&self) -> f64 {
        let s = self.n_sites();
        let mut best = 0.0f64;
        for a in 0..s {
            for b in (a + 1)..s {
                best = best.max(self.linkage(a, b));
            }
        }
        best
    }

    /// The dominant 2-way read bipartition (lightweight MEC / Lloyd iteration).
    /// This is the real discriminator: a genuine second allele is explained by a
    /// SINGLE read-partition that is consistent across MANY sites (high
    /// `consistency`); phase-shift noise in a homogeneous repeat produces only
    /// locally-correlated forks that no single global partition explains (low
    /// `consistency`), even though isolated site pairs may look linked.
    pub fn dominant_bipartition(&self) -> Bipartition {
        let ns = self.n_sites();
        if ns == 0 || self.n_reads == 0 {
            return Bipartition {
                assign: vec![0; self.n_reads],
                consistency: 1.0,
                minority_frac: 0.0,
                n_informative_sites: 0,
            };
        }
        // Seed from the best-covered, most-balanced site, then Lloyd-refine.
        // (Multi-restart over several seeds was measured on the robustness matrix
        // to give negligible count improvement and slightly worse single-allele
        // safety — the clustering seed is not the bottleneck, the split
        // confirmation is — so a single seed is kept for simplicity.)
        let seed_site = (0..ns)
            .max_by_key(|&s| {
                let sup = &self.sites[s].branch_support;
                let cov: u32 = sup.iter().sum();
                let mx = *sup.iter().max().unwrap_or(&0);
                (cov.saturating_sub(mx)) * 1000 + cov
            })
            .unwrap();
        self.score_partition(self.lloyd_from_seed(seed_site))
    }

    /// Lloyd-refine a bipartition seeded from `seed_site` (branch 0 → side 0,
    /// any other branch → side 1; uncovered reads → side 0). Returns the
    /// converged `side` assignment.
    #[allow(clippy::needless_range_loop)] // loops cross-index `side` and `calls`
    fn lloyd_from_seed(&self, seed_site: usize) -> Vec<i8> {
        let ns = self.n_sites();
        let mut side: Vec<i8> = (0..self.n_reads)
            .map(|r| match self.calls[r][seed_site] {
                0 => 0,
                MISSING => 0,
                _ => 1,
            })
            .collect();
        let maj = |side: &[i8], want: i8, s: usize| -> i16 {
            let nb = self.sites[s].branches.len();
            let mut c = vec![0u32; nb];
            for r in 0..self.n_reads {
                if side[r] == want {
                    let v = self.calls[r][s];
                    if v != MISSING {
                        c[v as usize] += 1;
                    }
                }
            }
            c.iter()
                .enumerate()
                .max_by_key(|&(_, &n)| n)
                .map(|(i, _)| i as i16)
                .unwrap_or(MISSING)
        };
        for _ in 0..12 {
            let maj0: Vec<i16> = (0..ns).map(|s| maj(&side, 0, s)).collect();
            let maj1: Vec<i16> = (0..ns).map(|s| maj(&side, 1, s)).collect();
            let mut changed = false;
            for r in 0..self.n_reads {
                let (mut a0, mut a1) = (0i32, 0i32);
                for s in 0..ns {
                    let v = self.calls[r][s];
                    if v == MISSING {
                        continue;
                    }
                    if v == maj0[s] {
                        a0 += 1;
                    }
                    if v == maj1[s] {
                        a1 += 1;
                    }
                }
                let new = if a1 > a0 { 1 } else { 0 };
                if new != side[r] {
                    side[r] = new;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        side
    }

    /// Score a GIVEN read bipartition (`side[r] ∈ {0, 1}`): its consistency
    /// (fraction of covered calls matching their side's per-site majority),
    /// minority fraction, and number of informative sites (where the two sides'
    /// majorities differ). Used both by `dominant_bipartition` and to confirm a
    /// split proposed by another engine (structural phasing) against the linkage
    /// signal.
    pub fn score_partition(&self, side: Vec<i8>) -> Bipartition {
        let ns = self.n_sites();
        let members =
            |want: i8| -> Vec<usize> { (0..self.n_reads).filter(|&r| side[r] == want).collect() };
        let (m0, m1) = (members(0), members(1));
        let maj_of = |members: &[usize], s: usize| -> i16 {
            let nb = self.sites[s].branches.len();
            let mut c = vec![0u32; nb];
            for &r in members {
                let v = self.calls[r][s];
                if v != MISSING {
                    c[v as usize] += 1;
                }
            }
            c.iter()
                .enumerate()
                .max_by_key(|&(_, &n)| n)
                .map(|(i, _)| i as i16)
                .unwrap_or(MISSING)
        };
        let mut agree = 0u32;
        let mut total = 0u32;
        let mut informative = 0usize;
        for s in 0..ns {
            let (j0, j1) = (maj_of(&m0, s), maj_of(&m1, s));
            if j0 != MISSING && j1 != MISSING && j0 != j1 {
                informative += 1;
            }
            for &r in &m0 {
                if self.calls[r][s] != MISSING {
                    total += 1;
                    if self.calls[r][s] == j0 {
                        agree += 1;
                    }
                }
            }
            for &r in &m1 {
                if self.calls[r][s] != MISSING {
                    total += 1;
                    if self.calls[r][s] == j1 {
                        agree += 1;
                    }
                }
            }
        }
        let consistency = if total == 0 {
            1.0
        } else {
            agree as f64 / total as f64
        };
        let n0 = m0.len();
        let n1 = m1.len();
        let minority_frac = if n0 + n1 == 0 {
            0.0
        } else {
            n0.min(n1) as f64 / (n0 + n1) as f64
        };
        Bipartition {
            assign: side,
            consistency,
            minority_frac,
            n_informative_sites: informative,
        }
    }

    /// Mean pairwise linkage over site pairs that share ≥ `min_shared` reads.
    /// A summary of how consistently the sites co-segregate.
    pub fn mean_linkage(&self, min_shared: usize) -> f64 {
        let s = self.n_sites();
        let mut sum = 0.0;
        let mut n = 0;
        for a in 0..s {
            for b in (a + 1)..s {
                let shared = (0..self.n_reads)
                    .filter(|&r| self.calls[r][a] != MISSING && self.calls[r][b] != MISSING)
                    .count();
                if shared >= min_shared {
                    sum += self.linkage(a, b);
                    n += 1;
                }
            }
        }
        if n == 0 { 0.0 } else { sum / n as f64 }
    }
}

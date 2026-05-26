//! Analysis helpers for interpreting [`Consensus`] output.
//!
//! These functions are sequence-general: none of them have any concept of a
//! repeat unit or genomic feature.  Domain-specific interpretation (repeat
//! count CI, expansion detection, pathogenicity thresholds) belongs in the
//! calling application.
//!
//! # Typical workflow
//!
//! ```no_run
//! use poa_consensus::{consensus, PoaConfig};
//! use poa_consensus::analysis::{consensus_confidence, should_call_multiallele,
//!                                count_credible_interval};
//!
//! # fn main() -> Result<(), poa_consensus::PoaError> {
//! # let reads: Vec<&[u8]> = vec![b"ACGT", b"ACGT", b"ACGT"];
//! let cfg = PoaConfig::default();
//! let result = consensus(&reads, 0, &cfg)?;
//!
//! // Quick go/no-go.
//! let conf = consensus_confidence(&result, cfg.min_allele_freq);
//! if conf.is_low_confidence() {
//!     eprintln!("low-confidence consensus: {:?}", conf);
//! }
//!
//! // Should we re-run with consensus_multi?
//! if should_call_multiallele(&result, 0.25) {
//!     eprintln!("competing allele detected — consider multi-allele mode");
//! }
//!
//! // Credible interval from per-observation estimates (e.g. per-read counts).
//! let per_read_counts = vec![40.0_f64, 41.0, 39.0, 40.0, 42.0];
//! let (lo, hi) = count_credible_interval(&per_read_counts, 0.95);
//! eprintln!("95% CI: [{lo:.1}, {hi:.1}]");
//! # Ok(())
//! # }
//! ```

use std::ops::Range;

use crate::types::{BubbleSite, Consensus};

// ─── Math utilities ───────────────────────────────────────────────────────────
//
// Zero-dependency implementations of the standard normal CDF and its inverse.
// Accuracy is sufficient for all use-cases here (< 1.5×10⁻⁷ for erf).

/// Abramowitz & Stegun 7.1.27 rational approximation; |error| < 1.5×10⁻⁷.
fn erf(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.327_591_1 * x.abs());
    let poly = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    let r = 1.0 - poly * (-x * x).exp();
    if x >= 0.0 { r } else { -r }
}

/// Standard normal CDF: P(Z ≤ x).
fn normal_cdf(x: f64) -> f64 {
    0.5 * (1.0 + erf(x / std::f64::consts::SQRT_2))
}

/// Inverse standard normal CDF (probit).
///
/// Uses Abramowitz & Stegun 26.2.22 rational approximation; |error| < 4.5×10⁻⁴.
/// Suitable for constructing confidence intervals at any conventional level.
fn probit(p: f64) -> f64 {
    let p = p.clamp(f64::EPSILON, 1.0 - f64::EPSILON);
    let sign = if p >= 0.5 { 1.0_f64 } else { -1.0_f64 };
    let q = p.min(1.0 - p);
    let t = (-2.0 * q.ln()).sqrt();
    let num = 2.515_517 + 0.802_853 * t + 0.010_328 * t * t;
    let den = 1.0 + 1.432_788 * t + 0.189_269 * t * t + 0.001_308 * t * t * t;
    sign * (t - num / den)
}

// ─── Raw data helpers ─────────────────────────────────────────────────────────

/// Minimum per-position read coverage across the consensus sequence.
///
/// Returns 0 for an empty consensus.  Positions with coverage below half the
/// read depth (`n_reads / 2`) are likely seed-only or partial-read artefacts;
/// use [`low_coverage_regions`] to locate them.
///
/// # Examples
/// ```
/// use poa_consensus::{consensus, PoaConfig};
/// use poa_consensus::analysis::min_coverage;
///
/// let reads: &[&[u8]] = &[b"ACGT", b"ACGT", b"ACGT"];
/// let result = consensus(reads, 0, &PoaConfig::default()).unwrap();
/// let m = min_coverage(&result);
/// assert!(m >= 1);
/// ```
pub fn min_coverage(consensus: &Consensus) -> u32 {
    consensus.coverage.iter().copied().min().unwrap_or(0)
}

/// Contiguous runs of positions where per-base coverage is below `threshold`.
///
/// A common threshold is `n_reads / 2` (half-depth), which identifies positions
/// supported by fewer reads than needed for a reliable majority call.  The
/// returned ranges are in ascending position order and never overlap.
///
/// # Examples
/// ```
/// use poa_consensus::analysis::low_coverage_regions;
/// use poa_consensus::{Consensus, GraphStats};
///
/// // Build a synthetic Consensus with a dip in the middle.
/// let cov = vec![5u32, 5, 1, 1, 5, 5];
/// let n = cov.len();
/// let result = Consensus {
///     sequence: b"AAAAAA".to_vec(),
///     coverage: cov,
///     path_weights: vec![5; n],
///     n_reads: 5,
///     graph_stats: GraphStats::default(),
///     gaps: vec![],
///     bubble_sites: vec![],
/// };
///
/// let low = low_coverage_regions(&result, 3);
/// assert_eq!(low, vec![2..4]);
/// ```
pub fn low_coverage_regions(consensus: &Consensus, threshold: u32) -> Vec<Range<usize>> {
    let mut result = Vec::new();
    let mut run_start: Option<usize> = None;
    for (i, &cov) in consensus.coverage.iter().enumerate() {
        if cov < threshold {
            run_start.get_or_insert(i);
        } else if let Some(s) = run_start.take() {
            result.push(s..i);
        }
    }
    if let Some(s) = run_start {
        result.push(s..consensus.coverage.len());
    }
    result
}

/// Per-arm fractions of reads at one bubble site, summing to 1.0.
///
/// The denominator is the total reads recorded across all arms at this site.
/// This may be less than `Consensus::n_reads` when reads predate the bubble or
/// are partial.  Returns all-zeros for a site with no recorded reads.
///
/// # Examples
/// ```
/// use poa_consensus::analysis::allele_fractions;
/// use poa_consensus::BubbleSite;
///
/// let site = BubbleSite {
///     consensus_pos: 10,
///     arm_read_counts: vec![6, 4],
///     arm_sequences: vec![b"CAG".to_vec(), b"CAGCAG".to_vec()],
///     is_structural: true,
/// };
/// let fracs = allele_fractions(&site);
/// assert!((fracs[0] - 0.6).abs() < 1e-9);
/// assert!((fracs[1] - 0.4).abs() < 1e-9);
/// ```
pub fn allele_fractions(site: &BubbleSite) -> Vec<f64> {
    let total: u32 = site.arm_read_counts.iter().sum();
    if total == 0 {
        return vec![0.0; site.arm_read_counts.len()];
    }
    site.arm_read_counts
        .iter()
        .map(|&c| c as f64 / total as f64)
        .collect()
}

// ─── Confidence and credible interval ─────────────────────────────────────────

/// Confidence interval on the mean of a set of independent scalar observations.
///
/// Given `n` observations with sample mean `x̄` and sample standard deviation
/// `s`, returns `(x̄ − z·s/√n, x̄ + z·s/√n)` where `z = Φ⁻¹((1+confidence)/2)`.
///
/// This is the **equivalence zone**: two tools whose outputs both fall inside
/// this interval are statistically indistinguishable given these observations —
/// the data cannot prefer one over the other.
///
/// `confidence` must be in (0, 1).  Typical values: 0.90, 0.95, 0.99.
///
/// Returns `(NaN, NaN)` for an empty slice and `(values[0], values[0])` for a
/// single observation (interval undefined, point estimate returned).
///
/// # Examples
/// ```
/// use poa_consensus::analysis::count_credible_interval;
///
/// // Ten per-read repeat-count estimates with little variance.
/// let counts = vec![40.0, 41.0, 39.0, 40.0, 42.0, 40.0, 41.0, 39.0, 40.0, 41.0];
/// let (lo, hi) = count_credible_interval(&counts, 0.95);
/// // Mean ≈ 40.3, CI ≈ [39.5, 41.1] — both tools in this range are equivalent.
/// assert!(lo < 40.3 && 40.3 < hi);
/// ```
///
/// # Interpretation
/// A result of `(38.5, 42.0)` means: even the theoretically optimal estimator
/// cannot determine the true count more precisely than this range, given these
/// observations.  Any two tools reporting values inside `[38.5, 42.0]` should
/// be considered equally correct.
pub fn count_credible_interval(values: &[f64], confidence: f64) -> (f64, f64) {
    assert!(
        (0.0..1.0).contains(&confidence),
        "confidence must be in (0, 1)"
    );
    match values.len() {
        0 => (f64::NAN, f64::NAN),
        1 => (values[0], values[0]),
        n => {
            let nf = n as f64;
            let mean = values.iter().sum::<f64>() / nf;
            let var = values.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (nf - 1.0);
            let std = var.sqrt();
            if std == 0.0 {
                return (mean, mean);
            }
            let z = probit((1.0 + confidence) / 2.0);
            let hw = z * std / nf.sqrt();
            (mean - hw, mean + hw)
        }
    }
}

/// Theoretical maximum probability of the integer-exact correct estimate.
///
/// Given `n` independent observations each with per-observation standard
/// deviation `sigma_per_obs`, returns:
///
/// ```text
/// P(correct) = 2·Φ(√n / (2·σ)) − 1
/// ```
///
/// This is the **Cramér-Rao ceiling**: no unbiased estimator can exceed this
/// probability of returning the exact correct integer, regardless of algorithm.
///
/// Returns 1.0 when `sigma_per_obs ≤ 0` (exact observations) or 0.0 when
/// `n == 0` (no observations).
///
/// # Examples
/// ```
/// use poa_consensus::analysis::max_achievable_accuracy;
///
/// // CAG×40, ONT R10: σ ≈ 1.37 units/read.
/// let acc_5  = max_achievable_accuracy(5,  1.37);  // ≈ 0.59
/// let acc_20 = max_achievable_accuracy(20, 1.37);  // ≈ 0.90
/// let acc_50 = max_achievable_accuracy(50, 1.37);  // ≈ 0.99
/// assert!(acc_5 < acc_20 && acc_20 < acc_50);
/// assert!(acc_50 > 0.98);
///
/// // Noisy locus: σ = 3.0 units/read, only 20 reads.
/// let hard = max_achievable_accuracy(20, 3.0);  // ≈ 0.55
/// // Even the optimal estimator is wrong ~45% of the time.
/// assert!(hard < 0.60);
/// ```
///
/// # Interpretation
/// Pass the empirical standard deviation of your per-read estimates as
/// `sigma_per_obs`.  The return value is the accuracy ceiling; two tools
/// within the [`count_credible_interval`] of each other have both hit this
/// ceiling and should be treated as equivalent.
pub fn max_achievable_accuracy(n: usize, sigma_per_obs: f64) -> f64 {
    if n == 0 {
        return 0.0;
    }
    if sigma_per_obs <= 0.0 {
        return 1.0;
    }
    let arg = (n as f64).sqrt() / (2.0 * sigma_per_obs);
    (2.0 * normal_cdf(arg) - 1.0).clamp(0.0, 1.0)
}

// ─── Allele / bubble helpers ──────────────────────────────────────────────────

/// Returns the first bubble site where the second-strongest arm has at least
/// `min_freq` of total reads (`Consensus::n_reads`), or `None` if no such site
/// exists.
///
/// A `Some` result is the primary signal to re-run with
/// [`PoaGraph::consensus_multi`]: a meaningful fraction of reads support a
/// different sequence at this position.
///
/// `min_freq` is applied against `n_reads`, consistent with how the graph's
/// `min_allele_freq` threshold is applied during bubble detection.
///
/// # Examples
/// ```
/// use poa_consensus::analysis::has_competing_allele;
/// use poa_consensus::{Consensus, GraphStats, BubbleSite};
///
/// let site = BubbleSite {
///     consensus_pos: 5,
///     arm_read_counts: vec![7, 3],   // 3/10 = 30% on the minority arm
///     arm_sequences: vec![],
///     is_structural: true,
/// };
/// let result = Consensus {
///     sequence: b"ACGT".to_vec(),
///     coverage: vec![10; 4],
///     path_weights: vec![7; 4],
///     n_reads: 10,
///     graph_stats: GraphStats::default(),
///     gaps: vec![],
///     bubble_sites: vec![site],
/// };
///
/// // 30% > 25%: a competing allele is present.
/// assert!(has_competing_allele(&result, 0.25).is_some());
/// // 35% threshold: 30% does not meet it.
/// assert!(has_competing_allele(&result, 0.35).is_none());
/// ```
pub fn has_competing_allele(consensus: &Consensus, min_freq: f64) -> Option<&BubbleSite> {
    if consensus.n_reads == 0 {
        return None;
    }
    let n = consensus.n_reads as f64;
    consensus.bubble_sites.iter().find(|site| {
        // Sort descending; the second entry is the strongest minority arm.
        let mut counts = site.arm_read_counts.clone();
        counts.sort_unstable_by(|a, b| b.cmp(a));
        counts.get(1).is_some_and(|&c| c as f64 / n >= min_freq)
    })
}

/// Returns `true` when at least one bubble site has a minority arm at or above
/// `min_freq`.
///
/// Thin wrapper over [`has_competing_allele`] for callers that only need the
/// boolean signal.  When `true`, re-run with
/// [`PoaGraph::consensus_multi`] to separate the alleles.
///
/// # Examples
/// ```
/// use poa_consensus::analysis::should_call_multiallele;
/// use poa_consensus::{Consensus, GraphStats, BubbleSite};
///
/// let result = Consensus {
///     sequence: b"ACGT".to_vec(),
///     coverage: vec![10; 4],
///     path_weights: vec![10; 4],
///     n_reads: 10,
///     graph_stats: GraphStats::default(),
///     gaps: vec![],
///     bubble_sites: vec![],        // no bubbles
/// };
/// assert!(!should_call_multiallele(&result, 0.25));
/// ```
pub fn should_call_multiallele(consensus: &Consensus, min_freq: f64) -> bool {
    has_competing_allele(consensus, min_freq).is_some()
}

// ─── Summary quality report ───────────────────────────────────────────────────

/// Summary quality indicators for a [`Consensus`].
///
/// Produced by [`consensus_confidence`].  Each field is independently
/// interpretable; `is_low_confidence` combines them into a single go/no-go.
///
/// # Field interpretation
///
/// | Field | High-confidence value | Action when flagged |
/// |---|---|---|
/// | `min_cov` | ≥ 2 | Low end-coverage; check for partial reads |
/// | `mean_cov` | ≥ `n_reads / 2` | Normal |
/// | `has_gaps` | `false` | Reads do not fully overlap the region |
/// | `competing_allele` | `false` | Re-run with `consensus_multi` |
/// | `low_cov_fraction` | < 0.1 | > 10% of positions may be unreliable |
/// | `single_support_fraction` | < 0.15 | High noise; consider stricter `min_reads` |
#[derive(Debug, Clone)]
pub struct ConsensusConfidence {
    /// Minimum per-position coverage across the consensus.
    pub min_cov: u32,
    /// Mean per-position coverage across the consensus.
    pub mean_cov: f64,
    /// Whether any coverage gap was detected (see [`CoverageGap`]).
    pub has_gaps: bool,
    /// Whether any bubble site has a minority arm above the `min_allele_freq`
    /// threshold.  If `true`, re-run with [`PoaGraph::consensus_multi`].
    pub competing_allele: bool,
    /// Fraction of consensus positions with coverage below half the read depth.
    /// Values above ~0.1 indicate widespread partial-read coverage.
    pub low_cov_fraction: f64,
    /// Fraction of graph nodes supported by exactly one read (from
    /// [`GraphStats`]).  Values above ~0.15 suggest the graph is noisy.
    pub single_support_fraction: f64,
}

impl ConsensusConfidence {
    /// Returns `true` when any red-flag indicator is set.
    ///
    /// Red flags: `min_cov < 2`, `has_gaps`, `competing_allele`,
    /// `low_cov_fraction > 0.1`, or `single_support_fraction > 0.15`.
    ///
    /// This is a heuristic; callers may inspect individual fields to apply
    /// domain-specific thresholds.
    pub fn is_low_confidence(&self) -> bool {
        self.min_cov < 2
            || self.has_gaps
            || self.competing_allele
            || self.low_cov_fraction > 0.1
            || self.single_support_fraction > 0.15
    }
}

/// Compute a [`ConsensusConfidence`] summary for a consensus result.
///
/// `min_allele_freq` is forwarded to [`has_competing_allele`] and should match
/// the value used when building the graph (default: 0.25).
///
/// # Examples
/// ```
/// use poa_consensus::{consensus, PoaConfig};
/// use poa_consensus::analysis::consensus_confidence;
///
/// let reads: &[&[u8]] = &[b"ACGTACGT", b"ACGTACGT", b"ACGTACGT",
///                          b"ACGTACGT", b"ACGTACGT"];
/// let cfg = PoaConfig::default();
/// let result = consensus(reads, 0, &cfg).unwrap();
/// let conf = consensus_confidence(&result, cfg.min_allele_freq);
///
/// // Clean identical reads: no gaps, no competing alleles.
/// assert!(!conf.has_gaps);
/// assert!(!conf.competing_allele);
/// assert!(!conf.is_low_confidence());
/// ```
pub fn consensus_confidence(consensus: &Consensus, min_allele_freq: f64) -> ConsensusConfidence {
    let n = consensus.coverage.len();
    let min_cov = min_coverage(consensus);
    let mean_cov = if n == 0 {
        0.0
    } else {
        consensus.coverage.iter().sum::<u32>() as f64 / n as f64
    };
    let half_depth = ((consensus.n_reads as f64 / 2.0).ceil() as u32).max(1);
    let low_cov_fraction = if n == 0 {
        0.0
    } else {
        consensus
            .coverage
            .iter()
            .filter(|&&c| c < half_depth)
            .count() as f64
            / n as f64
    };

    ConsensusConfidence {
        min_cov,
        mean_cov,
        has_gaps: !consensus.gaps.is_empty(),
        competing_allele: should_call_multiallele(consensus, min_allele_freq),
        low_cov_fraction,
        single_support_fraction: consensus.graph_stats.single_support_fraction,
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::GraphStats;

    fn make_consensus(
        coverage: Vec<u32>,
        n_reads: usize,
        bubble_sites: Vec<BubbleSite>,
    ) -> Consensus {
        let n = coverage.len();
        Consensus {
            sequence: vec![b'A'; n],
            coverage,
            path_weights: vec![1; n],
            n_reads,
            graph_stats: GraphStats::default(),
            gaps: vec![],
            bubble_sites,
        }
    }

    fn make_site(counts: Vec<u32>) -> BubbleSite {
        BubbleSite {
            consensus_pos: 0,
            arm_read_counts: counts,
            arm_sequences: vec![],
            is_structural: false,
        }
    }

    // ── erf / normal_cdf / probit ─────────────────────────────────────────────

    #[test]
    fn erf_known_values() {
        assert!((erf(0.0)).abs() < 1e-6);
        assert!((erf(1.0) - 0.842_701).abs() < 1e-5);
        assert!((erf(-1.0) + 0.842_701).abs() < 1e-5);
        assert!((erf(3.0) - 0.999_978).abs() < 1e-5);
    }

    #[test]
    fn probit_known_z_scores() {
        // probit(0.975) ≈ 1.96  (used for 95% CI)
        assert!((probit(0.975) - 1.96).abs() < 5e-3);
        // probit(0.95) ≈ 1.645  (used for 90% CI)
        assert!((probit(0.95) - 1.645).abs() < 5e-3);
        // probit(0.5) = 0  (median)
        assert!(probit(0.5).abs() < 1e-6);
        // symmetry
        assert!((probit(0.025) + probit(0.975)).abs() < 1e-6);
    }

    // ── min_coverage ─────────────────────────────────────────────────────────

    #[test]
    fn min_coverage_empty() {
        let c = make_consensus(vec![], 0, vec![]);
        assert_eq!(min_coverage(&c), 0);
    }

    #[test]
    fn min_coverage_uniform() {
        let c = make_consensus(vec![5, 5, 5], 5, vec![]);
        assert_eq!(min_coverage(&c), 5);
    }

    #[test]
    fn min_coverage_dip() {
        let c = make_consensus(vec![5, 5, 1, 5], 5, vec![]);
        assert_eq!(min_coverage(&c), 1);
    }

    // ── low_coverage_regions ─────────────────────────────────────────────────

    #[test]
    fn low_cov_none_below_threshold() {
        let c = make_consensus(vec![5, 5, 5, 5], 5, vec![]);
        assert!(low_coverage_regions(&c, 3).is_empty());
    }

    #[test]
    fn low_cov_all_below_threshold() {
        let c = make_consensus(vec![1, 1, 1], 5, vec![]);
        assert_eq!(low_coverage_regions(&c, 3), vec![0..3]);
    }

    #[test]
    fn low_cov_middle_dip() {
        let c = make_consensus(vec![5, 5, 1, 1, 5, 5], 5, vec![]);
        assert_eq!(low_coverage_regions(&c, 3), vec![2..4]);
    }

    #[test]
    fn low_cov_at_ends() {
        let c = make_consensus(vec![1, 5, 5, 1], 5, vec![]);
        assert_eq!(low_coverage_regions(&c, 3), vec![0..1, 3..4]);
    }

    #[test]
    fn low_cov_trailing_run() {
        let c = make_consensus(vec![5, 5, 1, 1], 5, vec![]);
        assert_eq!(low_coverage_regions(&c, 3), vec![2..4]);
    }

    // ── allele_fractions ──────────────────────────────────────────────────────

    #[test]
    fn allele_fractions_zero_total() {
        let s = make_site(vec![0, 0]);
        assert_eq!(allele_fractions(&s), vec![0.0, 0.0]);
    }

    #[test]
    fn allele_fractions_equal() {
        let s = make_site(vec![5, 5]);
        let f = allele_fractions(&s);
        assert!((f[0] - 0.5).abs() < 1e-9);
        assert!((f[1] - 0.5).abs() < 1e-9);
    }

    #[test]
    fn allele_fractions_unequal() {
        let s = make_site(vec![6, 4]);
        let f = allele_fractions(&s);
        assert!((f[0] - 0.6).abs() < 1e-9);
        assert!((f[1] - 0.4).abs() < 1e-9);
    }

    #[test]
    fn allele_fractions_sum_to_one() {
        let s = make_site(vec![3, 4, 3]);
        let f = allele_fractions(&s);
        assert!((f.iter().sum::<f64>() - 1.0).abs() < 1e-9);
    }

    // ── count_credible_interval ───────────────────────────────────────────────

    #[test]
    fn ci_empty_is_nan() {
        let (lo, hi) = count_credible_interval(&[], 0.95);
        assert!(lo.is_nan() && hi.is_nan());
    }

    #[test]
    fn ci_single_is_point() {
        let (lo, hi) = count_credible_interval(&[42.0], 0.95);
        assert_eq!(lo, 42.0);
        assert_eq!(hi, 42.0);
    }

    #[test]
    fn ci_identical_values() {
        let vals = vec![10.0; 5];
        let (lo, hi) = count_credible_interval(&vals, 0.95);
        assert_eq!(lo, 10.0);
        assert_eq!(hi, 10.0);
    }

    #[test]
    fn ci_contains_mean() {
        let vals = vec![40.0, 41.0, 39.0, 40.0, 42.0, 40.0, 41.0, 39.0, 40.0, 41.0];
        let (lo, hi) = count_credible_interval(&vals, 0.95);
        let mean = vals.iter().sum::<f64>() / vals.len() as f64;
        assert!(lo < mean && mean < hi);
    }

    #[test]
    fn ci_wider_at_lower_confidence() {
        let vals = vec![40.0, 41.0, 39.0, 42.0, 38.0];
        let (lo90, hi90) = count_credible_interval(&vals, 0.90);
        let (lo99, hi99) = count_credible_interval(&vals, 0.99);
        assert!(hi90 - lo90 < hi99 - lo99);
    }

    #[test]
    fn ci_narrows_with_more_observations() {
        let few: Vec<f64> = vec![40.0, 41.0, 39.0, 42.0, 38.0];
        let many: Vec<f64> = few.iter().cycle().take(50).copied().collect();
        let (lo_few, hi_few) = count_credible_interval(&few, 0.95);
        let (lo_many, hi_many) = count_credible_interval(&many, 0.95);
        assert!(hi_many - lo_many < hi_few - lo_few);
    }

    // ── max_achievable_accuracy ───────────────────────────────────────────────

    #[test]
    fn accuracy_zero_reads() {
        assert_eq!(max_achievable_accuracy(0, 1.0), 0.0);
    }

    #[test]
    fn accuracy_zero_sigma() {
        assert_eq!(max_achievable_accuracy(10, 0.0), 1.0);
    }

    #[test]
    fn accuracy_increases_with_depth() {
        let s = 1.37;
        let a5 = max_achievable_accuracy(5, s);
        let a20 = max_achievable_accuracy(20, s);
        let a50 = max_achievable_accuracy(50, s);
        assert!(a5 < a20 && a20 < a50);
    }

    #[test]
    fn accuracy_decreases_with_sigma() {
        let n = 20;
        assert!(max_achievable_accuracy(n, 0.5) > max_achievable_accuracy(n, 1.5));
    }

    #[test]
    fn accuracy_cag40_50reads() {
        // CAG×40, σ ≈ 1.37: should be close to 99%
        let acc = max_achievable_accuracy(50, 1.37);
        assert!(acc > 0.98 && acc <= 1.0);
    }

    #[test]
    fn accuracy_hard_locus() {
        // σ=3.0, n=20: should be around 55%
        let acc = max_achievable_accuracy(20, 3.0);
        assert!(acc > 0.50 && acc < 0.65);
    }

    // ── has_competing_allele / should_call_multiallele ────────────────────────

    #[test]
    fn no_bubbles_no_competing_allele() {
        let c = make_consensus(vec![10; 4], 10, vec![]);
        assert!(has_competing_allele(&c, 0.25).is_none());
        assert!(!should_call_multiallele(&c, 0.25));
    }

    #[test]
    fn minority_arm_below_threshold() {
        // 1/10 = 10% < 25% threshold
        let c = make_consensus(vec![10; 4], 10, vec![make_site(vec![9, 1])]);
        assert!(has_competing_allele(&c, 0.25).is_none());
    }

    #[test]
    fn minority_arm_above_threshold() {
        // 4/10 = 40% >= 25% threshold
        let c = make_consensus(vec![10; 4], 10, vec![make_site(vec![6, 4])]);
        assert!(has_competing_allele(&c, 0.25).is_some());
        assert!(should_call_multiallele(&c, 0.25));
    }

    #[test]
    fn returns_first_qualifying_site() {
        let low_site = make_site(vec![9, 1]); // 10% — below threshold
        let high_site = make_site(vec![6, 4]); // 40% — above threshold
        let c = make_consensus(vec![10; 4], 10, vec![low_site, high_site]);
        let found = has_competing_allele(&c, 0.25).unwrap();
        assert_eq!(found.arm_read_counts, vec![6, 4]);
    }

    #[test]
    fn zero_reads_returns_none() {
        let c = make_consensus(vec![], 0, vec![make_site(vec![0, 0])]);
        assert!(has_competing_allele(&c, 0.25).is_none());
    }

    // ── consensus_confidence ─────────────────────────────────────────────────

    #[test]
    fn confidence_clean_reads() {
        let c = make_consensus(vec![10; 8], 10, vec![]);
        let conf = consensus_confidence(&c, 0.25);
        assert_eq!(conf.min_cov, 10);
        assert!((conf.mean_cov - 10.0).abs() < 1e-9);
        assert!(!conf.has_gaps);
        assert!(!conf.competing_allele);
        assert_eq!(conf.low_cov_fraction, 0.0);
        assert!(!conf.is_low_confidence());
    }

    #[test]
    fn confidence_flags_competing_allele() {
        let c = make_consensus(vec![10; 4], 10, vec![make_site(vec![6, 4])]);
        let conf = consensus_confidence(&c, 0.25);
        assert!(conf.competing_allele);
        assert!(conf.is_low_confidence());
    }

    #[test]
    fn confidence_flags_low_coverage() {
        // More than 10% of positions below half-depth.
        let mut cov = vec![10u32; 20];
        cov[0] = 1;
        cov[1] = 1;
        cov[2] = 1; // 3/20 = 15% below threshold
        let c = make_consensus(cov, 10, vec![]);
        let conf = consensus_confidence(&c, 0.25);
        assert!(conf.low_cov_fraction > 0.1);
        assert!(conf.is_low_confidence());
    }
}

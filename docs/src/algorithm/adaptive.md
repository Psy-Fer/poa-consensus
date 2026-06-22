# Adaptive Two-Pass Mode

`consensus_adaptive` runs a first pass to build the graph and compute `GraphStats`, then
inspects those stats to decide whether a second pass is needed and what form it should take.

## Decision table

| Condition (checked in order) | Action | `AdaptiveAction` |
|---|---|---|
| 1-3 bubbles, minority arm ≥ `min_allele_freq × n` | Multi-allele split on pass-1 graph | `MultiAllele` |
| Consensus < 60% of median read length (banded only) | Retry unbanded | `TruncationRetry { recovered }` |
| `single_support_fraction > 0.30` | Rebuild with `min_coverage_fraction` raised to ≥ 0.6 | `NoisyTighten` |
| Coverage CV > 1.5 and mode is `Global` | Rebuild with `SemiGlobal` | `SemiGlobalFallback` |
| Otherwise | Return pass-1 result | `PassThrough` |

Conditions are checked in order; the first that matches fires.

## Return type

```rust
pub struct AdaptiveResult {
    pub consensuses: Vec<Consensus>,  // one or two elements
    pub action:      AdaptiveAction,
}
```

The `action` field records which branch fired. This is useful for logging filter reasons:
a truncation that was detected and corrected (`TruncationRetry { recovered: true }`) should
be reported differently from a clean pass. Inspect `action` rather than re-running
`diagnose` after the fact.

## Truncation detection

Banded DP on pure-repetitive sequence can silently converge to a shorter-than-correct
diagonal without the alignment ever approaching the band edge. The consensus is shorter than
it should be but no error is returned. The truncation detector catches this by comparing:

```
consensus_len / median_input_read_len < 0.60
```

If this fires and alignment was banded, a second graph is built with `band_width = 0`
(unbanded) and the result replaces the pass-1 consensus. `recovered` in the action is `true`
if the unbanded consensus clears the 0.60 threshold.

This is the workaround for the RFC1 5-mer repeat bug: using `pad = 20` bp unique flanking
context causes the banded-DP-truncated consensus to trigger the ratio check; the unbanded
retry recovers ~106× AAAAG (vs truth ~115×, within tolerance).

## Multi-allele detection

The bubble-count guard (`1-3 bubbles`) prevents accidental multi-allele splits on noisy
single-allele data. At ONT substitution rates (4-8%), substitution errors at the same
position across multiple reads can create a bubble above the 25% threshold. The upper bound
of 3 bubbles limits the split to plausible diploid structure; a graph with 15 bubbles is
almost certainly a noisy single allele, not 16 haplotypes.

For ONT data, raising `min_allele_freq` to 0.40 suppresses the most common false-positive
multi-allele calls from error-rate-driven bubbles.

## When to use `consensus_adaptive` vs `consensus_multi`

Use `consensus_adaptive` when you want the library to decide. Use `consensus_multi` directly
when you have already determined from prior analysis (e.g. from read length bimodality) that
the locus is heterozygous and you want the split unconditionally.

`consensus_adaptive` calls `consensus_multi` internally when it detects the right conditions,
so the two paths produce the same phasing result for clear heterozygous cases.

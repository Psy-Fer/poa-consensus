# Limitations

What poa-consensus cannot do, or does not do reliably, and what to do about it. Numbers
come from the benchmarks in `bench/` and simulated read sets (str-replicant ONT R10 and
HiFi presets); treat them as indications of scale, not guarantees.

## Repeat lengths at the edge of what the reads can resolve

A consensus caller can only be as right as its reads allow. For a repeat whose reads
disagree by a unit or so, there is a hard limit on how often *any* method can call the
exact length, and sometimes the reads point at the wrong length altogether.

- **The ceiling.** With `n` reads whose per-read repeat length spreads with standard
  deviation `σ` (in units), the best any unbiased method can do is call the exact integer
  with probability `2·Φ(√n / 2σ) − 1`. At 10 ONT reads with `σ ≈ 1.2` units that is about
  0.8: one locus in five is wrong whatever the caller. [`analysis::max_achievable_accuracy`]
  computes it, and [`analysis::count_credible_interval`] gives the interval for the true
  length from the per-read lengths.
- **Biased reads.** ONT reads in repeats are deletion-biased, so the reads' own estimate
  of a repeat length runs short. In simulation, the true length lies outside the 95%
  interval of the reads' lengths for roughly 4-10% of haplotypes. No method that follows
  the reads can get those right; a caller that does is lucky, or carries a built-in bias
  toward the longer length.
- **Where poa-consensus stands.** Its heaviest-path consensus leans toward the full-length
  path, which partly offsets the deletion bias: on simulated clean ONT expansions it calls
  the exact length 98.5% of the time against an average ceiling of 96.3%, and 91.7% for
  rounding the reads' mean length. It still misses some near-split loci, and the misses
  are not the same ones other callers miss. On a set of eight hard ONT `(CA)n` cases no
  caller tested (abPOA, SPOA, POASTA, MAFFT, Kalign, medoid, pileup, poa-consensus) got all
  eight, and a length estimator tuned to those cases with hindsight did no better than
  poa-consensus. Each fix that flipped one of these cases broke another.

**What to do.** Treat an exact repeat count as uncertain where
[support quality](library/diagnostics.md#per-base-support-quality) is below about Q10 inside
the repeat: there the reads split on the length, and a different caller giving a
different answer is expected, not a sign that either is broken. More reads raise the
ceiling (it grows with `√n`). Platform-specific length bias is best modelled by the repeat
caller downstream (a stutter / error model, as in dedicated STR genotypers), not by the
consensus.

## Remaining sensitivity to read order

Reads are added to the graph one at a time, so the order can matter. A read whose extra
length belongs in one place can align more cheaply by reusing a node another read created
elsewhere, carrying a wrong extra repeat unit onto the heaviest path. Tandem gap costs fix
the common form of this (a slip next to a repeat-like flank), but some cases remain,
notably a repeat directly next to a homopolymer whose length varies between reads. These
calls show low support quality at the affected bases. If a call matters, running with the
reads in a different order and comparing is a direct check.

## Support quality measures agreement, not truth

Support quality says how decisively the reads back each base over its strongest
alternative. An error every read shares (a systematic homopolymer undercall, or the
deletion bias above) still scores high. Inside a repeat, where an indel sits is not
determined by the sequence, so a disputed length shows as low support over part of the
repeat rather than at one exact base.

## Homopolymer compression and two-base runs

Homopolymer handling (on by default) occasionally calls one unit short on clean ONT
expansions of repeats whose units contain a two-base run, such as `GGC` or `CCG`: about 2
in 100 in simulation. The same mechanism is what fixes `(AAAAT)n`-type repeats, and no
guard tried separated the two without giving up far more correct calls. Run with
`--no-homopolymer` to compare when a GC-rich expansion matters. See
[When to turn these off](library/single-allele.md#when-to-turn-these-off).

## Very low depth over multi-kb loci

With only 2-4 reads over several kilobases, the adaptive band can occasionally align a
read along the wrong diagonal; in one simulated 5 kb locus with 2 reads the banded call
had several hundred more errors than the unbanded one. The truncation warning catches
gross cases. At such low depth, cross-check with an unbanded run (`band_width = 0`,
`adaptive_band = false`; CLI `--band-width 0 --no-adaptive-band`), memory permitting.

## Repeats without flanking sequence

- **Rotation ambiguity.** Reads entering a repeat at different phases (`GAA`, `AAG`,
  `AGA`) align ambiguously, and boundary trim can drop the start of the repeat when most
  reads are phase-shifted relative to the backbone. Anchor every read on unique flanking
  sequence ([`consensus_flanked`](library/flanking-anchor.md)).
- **Long expansions.** When most reads are shorter than the expansion, use
  [`consensus_flanked`](library/flanking-anchor.md), which builds from spanning reads, or
  [`bridged_consensus`](library/bridged-consensus.md) when no read spans the locus.

## Multi-allele calling

- Each allele is the heaviest path of its own graph. Best-fit selection, homopolymer
  handling and tandem gap costs apply to single-allele calls only.
- Phasing by shared differences is confused by the phase ambiguity of periodic repeats,
  so SNP haplotypes inside repeats, alleles a unit or two apart, and skewed minority
  alleles remain hard. See [Multi-Allele Phasing](algorithm/multi-allele.md#known-limitation).
- At ONT error rates a single-base bubble can pass a low `min_allele_freq` on a clean
  locus. Raise `min_allele_freq` to about 0.4 on ONT and use at least 10 reads per allele.

## Not a repeat genotyper

poa-consensus builds a consensus sequence and reports how well the reads support it. It
does not model repeat units, stutter, or platform length bias, and it does not report
repeat counts or genotype likelihoods. For those, pass its output to a repeat caller
(for example bladerunner) that does.

[`analysis::max_achievable_accuracy`]: https://docs.rs/poa-consensus/latest/poa_consensus/analysis/fn.max_achievable_accuracy.html
[`analysis::count_credible_interval`]: https://docs.rs/poa-consensus/latest/poa_consensus/analysis/fn.count_credible_interval.html

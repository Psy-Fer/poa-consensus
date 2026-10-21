# Single-Allele Consensus

## Basic usage

```rust
use poa_consensus::{consensus, select_seed, SeedSelection, PoaConfig};

let reads: Vec<&[u8]> = /* your reads */;

// Let the library pick the seed automatically.
let seed_idx = select_seed(&reads, &SeedSelection::Auto)?;
let result   = consensus(&reads, seed_idx, &PoaConfig::default())?;

println!("{}", String::from_utf8_lossy(&result.sequence));
println!("length: {} bp, depth: {}", result.sequence.len(), result.n_reads);
```

## Alignment mode

`AlignmentMode::SemiGlobal` is the default and should be used for extracted STR reads.

**Why not Global?** Global alignment has a subtle failure at homopolymer-flanked loci. If
the seed has more leading identical bases than most reads (e.g. a 17-A seed when the majority
have 15-16 A's), shorter reads must place their deletions somewhere inside the homopolymer
run. The aligner puts them at the end of the run rather than the start, because both
positions are equally penalised. This leaves the extra 1-2 A-nodes at the boundary of the
run on the only spine path, with no shortcut alternative for those reads. The boundary trim
approaches from the left and stops as soon as the first node has adequate coverage; the
extra nodes are unreachable by the trim. The result is a consensus 1-2 bases longer than
the true majority length.

Semi-global eliminates this: shorter reads take free terminal gaps and never traverse the
extra nodes, giving them genuine low coverage that the boundary trim removes.

## Coverage output

`result.coverage` is the per-position read depth along the consensus path:

```rust
let min_cov = result.coverage.iter().copied().min().unwrap_or(0);
let mean_cov = result.coverage.iter().sum::<u32>() as f64 / result.coverage.len() as f64;
println!("min coverage: {min_cov}, mean coverage: {mean_cov:.1}");
```

Coverage can drop below the mean at coverage gaps (see
[Bridged Consensus](bridged-consensus.md)) and at boundaries where partial reads don't reach.

## How the consensus is chosen

With the default `ConsensusMode::BestFit`, `consensus()` builds the graph on the
median-length read (the `seed_idx` you pass is only bounds-checked) and extracts two
candidates:

- the **heaviest path** (abPOA-style: follow the heaviest out-edge from each node, tie-broken
  by downstream weight), with low-coverage leading/trailing nodes trimmed;
- the **majority-frequency** consensus (the plurality base of each MSA column, counting read
  deletions explicitly).

If they differ, each candidate is realigned against every read and both indel measures are
totalled: indel bases, and the affine indel penalty (every insert or delete event costs
`|gap_open|` and every indel base `|gap_extend|`). Mismatches are not counted (on noisy
reads they are mostly sequencing error). The majority-frequency candidate replaces the
heaviest path only if it is strictly better on **both**. Each measure alone fails on
repeats: counting bases alone rewards a frame-broken candidate halfway between two repeat
lengths (one base from every read), and the affine penalty alone breaks base-count ties
toward the shorter candidate, because a longer read can absorb several differences into
one insertion run.

`ConsensusMode::MajorityFrequency` skips the comparison and always returns the column
majority.

## Length variation

With length variation from **errors** (reads that are ±1 repeat unit due to error), both
candidates follow the plurality length: a minority of reads with an extra or missing unit
does not shift the consensus.

With two genuine alleles (e.g. 10 × CAG×20 and 10 × CAG×30), a single-allele consensus
reports one length and the other appears only as a bubble site. Use `consensus_multi` for
diploid loci (see [Multi-Allele Consensus](multi-allele.md)), or check
`diagnose(...).structural_competing` to detect the case.

## Homopolymer handling

Inside long, length-variable homopolymer runs, where a single interrupting base sits is not
fixed by the alignment. Take reads of `(A)12 G (A)11` and `(A)11 G (A)12`, or an
`(AAAAT)n` repeat whose A runs vary in length. Under the 2/-4/-4/-3 scoring the aligner
prefers two G/A mismatches over a deletion plus an insertion, so the one biological G is
split across two MSA columns and can be lost or duplicated in the consensus. Two
`PoaConfig` options address this, and both are on by default:

- `homopolymer_min_run: usize` (default `2`; `0` = off). Runs of at least this many
  identical bases are aligned as a single "poly" symbol, and each run's length is
  re-estimated afterwards as the median over the reads that support it. `1` compresses
  every run; `2` uses the `A, C, G, T, polyA..polyT` alphabet (single bases stay plain).
- `homopolymer_gap: Option<(i32, i32)>` (default `Some((-2, -1))`; `None` = off).
  `Some((open, extend))` replaces `gap_open` / `gap_extend` for indels in a homopolymer
  context, so a base shifted within a run aligns as cheap indels rather than a mismatch
  pair.

The plain engine consensus is always built too and is the fallback. Each enabled option
builds an extra candidate, which is **spliced** into the plain consensus only inside
homopolymer runs (so calls elsewhere, such as a CAG unit count, never change). The splice
is kept only if it fits the reads **strictly better** than the plain consensus (judged with
mismatches counted). Otherwise you get the plain result.

Scope: applies to `consensus`, and through it to `consensus_flanked` and
`bridged_consensus`. `consensus_multi` and `consensus_multi_flanked` ignore both options.

Cost is about 2.4x the plain engine's CPU time (the plain engine is itself 2.5-4x faster
than in v0.6.0, so the default is still faster than v0.6.0 was). It matters mainly on noisy
reads (ONT); on clean HiFi it rarely changes anything. Turn it off to save CPU on clean
HiFi data, or to get the plain engine's call (CLI: `--no-homopolymer`, or
`--homopolymer-gap off` to drop only the gap costs). See also
[When to turn these off](#when-to-turn-these-off):

```rust
use poa_consensus::{consensus, PoaConfig};

let config = PoaConfig {
    homopolymer_min_run: 0,
    homopolymer_gap: None,
    ..PoaConfig::default()
};
let result = consensus(&reads, 0, &config)?;
```

## Repeat-unit slips (tandem gap costs)

Inside a tandem repeat the commonest real difference between reads is one more or one
fewer copy of the unit: polymerase slippage, and sequencing stutter. Plain affine gap
costs price that like any other indel, so reads reach for alignments that only look
cheaper. A read whose extra length is really a `CT` unit in a `(CT)n` flank can slide
across the boundary into a neighbouring `(CA)n` repeat with two mismatches and pass
through an extra `CA` node another read created, because matching an existing node is
free. If enough reads do that, the heaviest path carries a wrong extra unit, and whether
it happens depends on read order.

`tandem_gap` (default `Some((-2, -1))`) makes the honest alignment cheap: an indel that
copies its neighbour in a locally periodic stretch of the read (a 2-6 base unit
immediately followed by a copy of itself) pays these costs instead of `gap_open` /
`gap_extend`. The test is local to each alignment step; no repeat tract or unit count is
identified, and nothing about repeats is reported.

Like homopolymer handling it only builds a candidate. The candidate replaces the plain
call only if it fits the reads strictly better, and it is built only when the plain call
has a base with [support quality](diagnostics.md#per-base-support-quality) below 10
(somewhere the reads disagree), so clean loci cost nothing extra. The guard matters:
applied unguarded, the same costs pull a genuine length split toward the shorter
length. Like the homopolymer options it applies to `consensus`, `consensus_flanked` and
`bridged_consensus`; `consensus_multi` ignores it.

## When to turn these off

The defaults give more correct calls overall, so leave them on unless one of these
applies:

- **Reproducing a plain POA call.** For a like-for-like comparison with another POA tool,
  or with earlier poa-consensus output, turn all three off (`homopolymer_min_run: 0`,
  `homopolymer_gap: None`, `tandem_gap: None`; CLI `--no-homopolymer --tandem-gap off`).
  Earlier versions also differ in the best-fit rule, so v0.6.0 output is not reproduced
  exactly.
- **Speed on clean data.** Homopolymer handling costs about 2.4x the plain engine's CPU on
  every call. On HiFi reads it rarely changes a call; `--no-homopolymer` saves that time.
  Tandem gap costs are only paid where the reads disagree, so leaving them on is cheap.
- **Checking a call that matters.** Run once with the options off. If the call changes,
  the reads are contested there (the support quality will usually say so too) and the
  locus deserves a closer look.
- **Known weak spots.** Homopolymer compression can call one unit short on clean ONT
  expansions of repeats whose units contain two-base runs (`GGC`, `CCG`): about 2 in 100
  in simulation. On heavily smeared (mosaic) expansions tandem gap costs recover calls
  that otherwise break, but where both versions give an intact call they are roughly a
  wash on length.

```rust
use poa_consensus::{consensus, PoaConfig};

// the plain engine: no homopolymer handling, no tandem gap costs
let plain = PoaConfig {
    homopolymer_min_run: 0,
    homopolymer_gap: None,
    tandem_gap: None,
    ..PoaConfig::default()
};
let result = consensus(&reads, 0, &plain)?;
```

## Checking result quality

```rust
use poa_consensus::{diagnose, DiagnoseConfig};

let warnings = diagnose(&result, &DiagnoseConfig::default());
if !warnings.is_clean() {
    for (is_warning, msg) in warnings.messages("my_locus") {
        eprintln!("{}: {msg}", if is_warning { "WARN" } else { "NOTE" });
    }
}

// Structured access
if let Some(trunc) = &warnings.truncation_suspected {
    eprintln!("truncation: consensus {}bp vs median read {}bp (ratio {:.2})",
        trunc.consensus_len, trunc.median_read_len, trunc.ratio);
}
```

See [Diagnostics and Warnings](diagnostics.md) for all available checks.

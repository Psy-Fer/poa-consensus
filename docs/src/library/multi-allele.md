# Multi-Allele Consensus

## Basic usage

```rust
use poa_consensus::{consensus_multi, select_seed, SeedSelection, PoaConfig};

let reads: Vec<&[u8]> = /* reads from a diploid STR locus */;
let seed_idx = select_seed(&reads, &SeedSelection::Auto)?;
let alleles  = consensus_multi(&reads, seed_idx, &PoaConfig::default())?;

for (i, allele) in alleles.iter().enumerate() {
    println!("allele {}: {} bp ({} reads)",
        i, allele.sequence.len(), allele.n_reads);
}
```

`consensus_multi` returns one `Consensus` per detected allele. If no second allele is
found the result is a single-element `Vec`. This is close to, but not identical to,
`consensus`: each multi-allele output (including the single-allele fallback) is the
heaviest-path consensus of its read group, so `consensus_mode` (best-fit) and the
homopolymer options do not apply, and `read_indices` is populated. Like `consensus`, the
engine seeds on the median-length read; `seed_idx` is only bounds-checked.

Phasing enforces the `min_allele_freq` frequency floor: an allele group whose read
support falls below `min_allele_freq` of the locus is not called as a separate allele.
This rejects stutter and valley over-splits, where sequencing noise on a homogeneous
repeat would otherwise fragment one true allele into several near-duplicates.

## Read indices

Each returned `Consensus` carries `read_indices: Vec<usize>`, which holds the indices into
the original `reads` slice that contributed to that allele. This allows you to assign per-read
metadata to the correct haplotype without re-running alignment:

```rust
for (allele_id, allele) in alleles.iter().enumerate() {
    for &read_idx in &allele.read_indices {
        println!("read {} → allele {}", read_idx, allele_id);
    }
}
```

`read_indices` is always populated by `consensus_multi`, including when it returns a
single allele (then it lists every non-empty input read). It is empty only for outputs of
`consensus()`, `consensus_flanked()` and `bridged_consensus()`. `consensus_multi_flanked`
reports indices into the original `reads` slice even when it phases only the spanning
subset.

## Depth requirements

`min_reads` (default 3) applies twice. The whole input needs at least `min_reads` non-empty
reads, or `consensus_multi` returns `Err(PoaError::InsufficientDepth)`. During phasing, a
candidate group smaller than `min_reads` is normally rejected or merged into its
neighbour, producing fewer alleles. If a confirmed group still ends up below `min_reads`,
the call returns `Err(InsufficientDepth)` rather than a low-depth allele. For reliable
diploid calls, plan for at least 5-10 reads per allele.

## Multi-allele with explicit config

For a diploid STR locus with known depth:

```rust
let mut cfg = PoaConfig::default();
cfg.min_reads       = 5;    // require 5 reads per allele minimum
cfg.min_allele_freq = 0.30; // ONT: raise from default 0.2 to suppress error-rate bubbles

let alleles = consensus_multi(&reads, seed_idx, &cfg)?;
```

## Checking allele fractions

```rust
use poa_consensus::analysis::allele_fractions;

for allele in &alleles {
    for site in &allele.bubble_sites {
        let fracs = allele_fractions(site);
        println!("bubble at {}: arms = {:?}", site.consensus_pos, fracs);
    }
}
```

## When two alleles become one

`consensus_multi` falls back to a single consensus in three cases:

1. No bubble above `min_allele_freq` after building the graph (truly homozygous locus or
   insufficient depth to see the second allele).
2. Only one phasing group has `>= min_reads` reads after partitioning, or the confirmation
   step finds that the groups do not yield genuinely different alleles (an over-split).
3. No structural (length-changing) bubble appears on the banded graph. Before giving up on
   a length split, `consensus_multi` rebuilds once fully unbanded and looks again.

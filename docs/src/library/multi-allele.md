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

`consensus_multi` returns one `Consensus` per detected allele. If no heterozygous bubble
is found the result is a single-element `Vec` equivalent to `consensus`.

## Read indices

Each returned `Consensus` carries `read_indices: Vec<usize>` -- the indices into the
original `reads` slice that contributed to that allele. This allows you to assign per-read
metadata to the correct haplotype without re-running alignment:

```rust
for (allele_id, allele) in alleles.iter().enumerate() {
    for &read_idx in &allele.read_indices {
        println!("read {} → allele {}", read_idx, allele_id);
    }
}
```

For single-allele outputs (where `consensus_multi` found no bubble), `read_indices` is
empty, meaning all reads contributed.

## Depth requirements

`min_reads` (default 3) is enforced **per allele**. If the minority allele has fewer than
`min_reads` reads after phasing it is merged into the dominant group, producing a
single-allele result. For reliable diploid calls, plan for at least 5--10 reads per allele.

## Multi-allele with explicit config

For a diploid STR locus with known depth:

```rust
let mut cfg = PoaConfig::default();
cfg.min_reads       = 5;    // require 5 reads per allele minimum
cfg.min_allele_freq = 0.30; // ONT: raise from default 0.25 to suppress error-rate bubbles

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

## Stateful multi-allele

```rust
let mut graph = PoaGraph::new(reads[seed_idx], PoaConfig::default())?;
for (i, read) in reads.iter().enumerate() {
    if i == seed_idx { continue; }
    graph.add_read(read)?;
}

// Inspect graph before splitting
let stats = graph.stats();
println!("structural bubbles: {}", stats.bubble_count);

let alleles = graph.consensus_multi()?;
```

## When two alleles become one

`consensus_multi` falls back to a single consensus in three cases:

1. No bubble above `min_allele_freq` after building the graph (truly homozygous locus or
   insufficient depth to see the second allele).
2. Only one phasing group has `>= min_reads` reads after partitioning.
3. The graph is built with a narrow band that cannot represent a large Insert-arm bubble
   for the allele-length difference -- in this case `consensus_multi` automatically rebuilds
   with unbanded alignment and tries again before falling back.

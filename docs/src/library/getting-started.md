# Getting Started

## Add to your project

```toml
[dependencies]
poa-consensus = "0.1"
```

For visualisation support:

```toml
[dependencies]
poa-consensus = { version = "0.1", features = ["plot"] }
```

## Minimal example

```rust
use poa_consensus::{consensus, PoaConfig};

fn main() -> Result<(), poa_consensus::PoaError> {
    let reads: Vec<&[u8]> = vec![
        b"CATCATCATCAT",
        b"CATCATCATCAT",
        b"CATCATCATCAT",
        b"CATCATCATCAT",
        b"CATCATCATCATCAT",  // one longer read
    ];

    let result = consensus(&reads, 0, &PoaConfig::default())?;
    println!("{}", String::from_utf8_lossy(&result.sequence));
    println!("depth: {}", result.n_reads);
    println!("bubbles: {}", result.graph_stats.bubble_count);
    Ok(())
}
```

## Choosing an API

There are two API styles: **functional** (convenience wrappers) and **stateful** (direct
graph access).

### Functional

```rust
use poa_consensus::{consensus, consensus_multi, consensus_adaptive, PoaConfig};

// Single-allele
let c = consensus(&reads, seed_idx, &config)?;

// Multi-allele (returns one Consensus per detected allele)
let alleles = consensus_multi(&reads, seed_idx, &config)?;

// Adaptive: decide multi vs single based on graph statistics
let result = consensus_adaptive(&reads, seed_idx, &config)?;
let alleles = result.consensuses;
```

### Stateful

```rust
use poa_consensus::{PoaGraph, PoaConfig};

let mut graph = PoaGraph::new(reads[seed_idx], PoaConfig::default())?;
for (i, read) in reads.iter().enumerate() {
    if i == seed_idx { continue; }
    graph.add_read(read)?;
}

// Inspect state mid-build
let stats = graph.stats();
println!("bubble count after {} reads: {}", graph.n_reads(), stats.bubble_count);

// Extract consensus
let consensus = graph.consensus()?;

// Or multi-allele
let alleles = graph.consensus_multi()?;
```

The stateful API is preferable for high-throughput pipelines: the graph object pre-allocates
internal buffers and reuses them across calls, avoiding repeated heap allocation. It also
allows inspecting `GraphStats` mid-build to decide whether to continue adding reads or bail
early.

## `PoaConfig` defaults

```rust
PoaConfig {
    band_width: 50,                         // minimum band; 0 = unbanded
    adaptive_band: true,                    // abPOA formula: w = 10 + 0.01 * max(read_len, graph_nodes)
    adaptive_band_b: 10,
    adaptive_band_f: 0.01,
    match_score: 1,
    mismatch_score: -1,
    gap_open: -2,
    gap_extend: -1,
    min_coverage_fraction: 0.5,             // strict majority
    min_allele_freq: 0.25,
    min_reads: 3,
    alignment_mode: AlignmentMode::SemiGlobal,
    consensus_mode: ConsensusMode::HeaviestPath,
    warn_on_long_unbanded: true,
}
```

The defaults are tuned for STR loci with HiFi or ONT reads. For ONT data specifically,
raise `min_allele_freq` to 0.40 to reduce false-positive multi-allele calls caused by
the higher substitution rate.

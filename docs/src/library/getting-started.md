# Getting Started

## Add to your project

```toml
[dependencies]
poa-consensus = "0.6"
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

## The functional API

The public API is purely functional: you pass a slice of reads and get back a
`Consensus` (or several, for multi-allele loci). There is no graph object to
manage.

```rust
use poa_consensus::{consensus, consensus_multi, PoaConfig};

// Single-allele
let c = consensus(&reads, seed_idx, &config)?;

// Multi-allele (returns one Consensus per detected allele)
let alleles = consensus_multi(&reads, seed_idx, &config)?;
```

`consensus()` self-seeds on the median-length read internally, so the graph is
built from a representative read regardless of which index you pass. The
`seed_idx` argument is still validated for bounds (it must point at a real read),
but it does not bias the result.

The returned consensus is a **best-fit** consensus: `consensus()` computes both a
heaviest-path consensus and a majority-frequency (MSA-column) consensus and keeps
whichever the reads better support. You can force the majority-frequency result
with `ConsensusMode::MajorityFrequency` (see below).

## `PoaConfig` defaults

```rust
PoaConfig {
    band_width: 50,                         // minimum band; 0 = unbanded
    adaptive_band: true,                    // abPOA formula: w = 10 + 0.01 * max(read_len, graph_nodes)
    adaptive_band_b: 10,
    adaptive_band_f: 0.01,
    match_score: 2,
    mismatch_score: -4,
    gap_open: -4,
    gap_extend: -3,
    min_coverage_fraction: 0.0,             // 0.0 = use the (n/2 + 1).max(2) majority floor
    min_allele_freq: 0.2,                   // minimum arm frequency to call a second allele
    min_reads: 3,                           // Err(InsufficientDepth) below this
    alignment_mode: AlignmentMode::SemiGlobal,
    consensus_mode: ConsensusMode::BestFit,  // "best-fit" (see above)
    warn_on_long_unbanded: true,            // warn on stderr for unbanded multi-kb reads
    phasing_bubble_min_span: 10,            // min arm span (bp) for a structural (length) bubble
}
```

The defaults are tuned for STR loci with HiFi or ONT reads. For ONT data specifically,
raise `min_allele_freq` to 0.40 to reduce false-positive multi-allele calls caused by
the higher substitution rate.

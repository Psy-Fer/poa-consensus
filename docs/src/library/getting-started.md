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

`consensus()` and `consensus_multi()` self-seed on the median-length read internally,
so the graph is built from a representative read regardless of which index you pass.
The `seed_idx` argument is still validated for bounds (it must point at a real read),
and orientation helpers such as `auto_orient` use it as the reference strand, but it does
not bias the result.

The returned consensus is a **best-fit** consensus: `consensus()` computes both a
heaviest-path consensus and a majority-frequency (MSA-column) consensus and keeps
whichever the reads better support: each candidate is realigned to every read, and the
majority-frequency one replaces the heaviest path only if it needs strictly fewer indel
bases and a lower affine indel penalty. You can force the majority-frequency result
with `ConsensusMode::MajorityFrequency` (see below).

## `PoaConfig` defaults

```rust
PoaConfig {
    band_width: 50,                          // band half-width floor (see below)
    adaptive_band: true,                     // abPOA formula: w = 10 + 0.01 * read_len
    adaptive_band_b: 10,
    adaptive_band_f: 0.01,
    match_score: 2,
    mismatch_score: -4,
    gap_open: -4,                            // a gap of k bases costs gap_open + k * gap_extend
    gap_extend: -3,
    min_coverage_fraction: 0.0,              // 0.0 = use the (n/2 + 1).max(2) majority floor
    min_allele_freq: 0.2,                    // minimum allele frequency for multi-allele calls
    min_reads: 3,                            // Err(InsufficientDepth) below this
    alignment_mode: AlignmentMode::SemiGlobal,
    consensus_mode: ConsensusMode::BestFit,  // "best-fit" (see above)
    warn_on_long_unbanded: true,             // warn on stderr for unbanded reads > ~1 kb
    phasing_bubble_min_span: 10,             // min arm span (nodes) for a structural bubble
    homopolymer_min_run: 2,                  // align runs >= 2 bases as one symbol (0 = off)
    homopolymer_gap: Some((-2, -1)),         // gap open/extend inside homopolymers (None = off)
    tandem_gap: Some((-2, -1)),             // repeat-unit slip gap costs (None = off)
}
```

Build configs with struct-update syntax so new fields don't break your code:

```rust
let config = PoaConfig { min_reads: 5, ..PoaConfig::default() };
```

**Banding.** The band half-width is `max(band_width, adaptive_band_b + adaptive_band_f ×
read_len, 4)`, recomputed per read. Fully unbanded alignment needs both `band_width: 0` and
`adaptive_band: false`; `band_width: 0` alone still uses the adaptive band. Reads much
shorter than the graph (partial reads) are aligned unbanded automatically.

The defaults are tuned for STR loci with HiFi or ONT reads. For ONT data, raise
`min_allele_freq` to about 0.4 to reduce false-positive multi-allele calls caused by the
higher substitution rate. Homopolymer handling is on by default (single-allele only); see
[Homopolymer handling](single-allele.md#homopolymer-handling) for what it does and when
to turn it off.

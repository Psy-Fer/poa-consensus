# Adaptive Mode

`consensus_adaptive` is a two-pass wrapper that uses `GraphStats` from the first pass to
decide automatically whether a second pass is needed and what form it should take.

## Basic usage

```rust
use poa_consensus::{consensus_adaptive, select_seed, SeedSelection, PoaConfig};

let reads: Vec<&[u8]> = /* your reads */;
let seed_idx = select_seed(&reads, &SeedSelection::Auto)?;
let result   = consensus_adaptive(&reads, seed_idx, &PoaConfig::default())?;

// result.action tells you what happened
println!("action: {:?}", result.action);

for (i, allele) in result.consensuses.iter().enumerate() {
    println!("allele {}: {} bp", i, allele.sequence.len());
}
```

## Using the action for logging

```rust
use poa_consensus::AdaptiveAction;

match &result.action {
    AdaptiveAction::PassThrough => {
        // clean single-allele result, no second pass needed
    }
    AdaptiveAction::MultiAllele => {
        println!("{} alleles detected", result.consensuses.len());
    }
    AdaptiveAction::TruncationRetry { recovered } => {
        if *recovered {
            println!("truncation detected and corrected by unbanded retry");
        } else {
            eprintln!("truncation detected; unbanded retry did not recover full length");
        }
    }
    AdaptiveAction::NoisyTighten => {
        println!("high noise: coverage threshold tightened");
    }
    AdaptiveAction::SemiGlobalFallback => {
        println!("partial reads detected: switched to semi-global alignment");
    }
}
```

## When to use adaptive vs direct multi-allele

Use `consensus_adaptive` when:
- You don't know in advance whether the locus is homozygous or heterozygous
- You want truncation detection as a safety net for repetitive sequence
- You want a single call that handles all common failure modes

Use `consensus_multi` directly when:
- You have already determined the locus is heterozygous (e.g. from bimodal read lengths)
- You want the phasing unconditionally without the bubble-count guard

## Tuning for ONT data

ONT substitution rates (4--8%) cause error-driven bubbles that can look like alleles:

```rust
let mut cfg = PoaConfig::default();
cfg.min_allele_freq = 0.40;  // raise from 0.25 to suppress error-rate bubbles
let result = consensus_adaptive(&reads, seed_idx, &cfg)?;
```

With `min_allele_freq = 0.40` and depth ≥ 10, at least 4 reads must support the minority
arm before the multi-allele split fires. An error that affects 1/10 reads (10%) is below
the threshold and treated as noise.

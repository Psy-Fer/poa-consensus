# Diagnostics and Warnings

`diagnose` runs four independent checks on a completed `Consensus` and returns a
`ConsensusWarnings` struct with structured output.

## Basic usage

```rust
use poa_consensus::{diagnose, DiagnoseConfig};

let warnings = diagnose(&result, &DiagnoseConfig::default());

// Single go/no-go check
if warnings.is_clean() {
    println!("consensus looks good");
} else {
    // Human-readable messages (label appears in the output text)
    for (is_warning, msg) in warnings.messages("HTT") {
        let level = if is_warning { "WARN" } else { "NOTE" };
        eprintln!("{level}: {msg}");
    }
}
```

## The four checks

### 1. Low depth

```rust
if let Some(w) = &warnings.low_depth {
    println!("depth {} below threshold {}", w.depth, w.threshold);
}
```

Fires when `n_reads < DiagnoseConfig::min_depth` (default 5). Low depth means the coverage
threshold (`min_cov`) is set too conservatively for the data, and the consensus may reflect
noise rather than the true sequence.

### 2. Coverage gaps

```rust
if warnings.has_coverage_gaps {
    for gap in &result.gaps {
        println!("gap at {:?}: {:?}", gap.start..gap.end, gap.kind);
    }
}
```

Fires when `result.gaps` is non-empty. A gap means reads did not fully span the region. The
consensus sequence at the gap position is either seed-only (Spanning gap) or not present at
all (Unknown gap from `bridged_consensus`). The gap region should not be used for unit
counting or variant calling.

### 3. Interior low support

```rust
if let Some(w) = &warnings.interior_low_support {
    println!("min interior weight fraction: {:.2}", w.min_fraction);
}
```

Fires when the minimum `path_weights` fraction (weight / n_reads) in the **middle 60%** of
the consensus is below `interior_support_threshold` (default 0.15). Very low path weight in
the interior suggests the heaviest path is stitching across a coverage gap or picking a
low-support route through a complex bubble structure.

### 4. Truncation

```rust
if let Some(w) = &warnings.truncation_suspected {
    println!("ratio {:.2}: consensus {}bp vs median read {}bp",
        w.ratio, w.consensus_len, w.median_read_len);
}
```

Fires when `consensus_len / median_input_read_len < truncation_ratio_threshold` (default
0.60). A ratio below 0.60 indicates the consensus is substantially shorter than the reads
that built it. The most common cause is banded DP converging to the wrong diagonal on a
highly repetitive locus.

When this fires in `consensus_adaptive`, an unbanded retry is attempted automatically. When
it fires on a standalone `consensus` call, the caller should retry with `band_width = 0`.

## `DiagnoseConfig`

```rust
pub struct DiagnoseConfig {
    pub min_depth: usize,                    // default 5
    pub interior_support_threshold: f64,     // default 0.15
    pub truncation_ratio_threshold: f64,     // default 0.60; set 0.0 to disable
    pub structural_allele_freq: f64,         // default 0.25
}
```

## Analysis helpers

For more targeted analysis without the full `diagnose` call:

```rust
use poa_consensus::analysis::{
    min_coverage, low_coverage_regions, has_competing_allele,
    should_call_multiallele, consensus_confidence,
};

// Minimum per-position coverage
let floor = min_coverage(&result);

// Positions below half-depth
let low = low_coverage_regions(&result, result.n_reads as u32 / 2);

// Is there a competing allele above 25%?
if let Some(site) = has_competing_allele(&result, 0.25) {
    println!("bubble at pos {}: {:?}", site.consensus_pos, site.arm_read_counts);
}

// Quick boolean for pipeline branching
if should_call_multiallele(&result, 0.25) {
    // re-run with consensus_multi
}

// Summary struct for logging
let conf = consensus_confidence(&result, 0.25);
println!("min_cov={} mean_cov={:.1} competing_allele={}",
    conf.min_cov, conf.mean_cov, conf.competing_allele);
```

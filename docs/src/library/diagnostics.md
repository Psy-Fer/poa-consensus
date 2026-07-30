# Diagnostics and Warnings

`diagnose` runs five independent checks on a completed `Consensus` and returns a
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

## The five checks

### 1. Low depth

```rust
if let Some(w) = &warnings.low_depth {
    println!("depth {} below recommended {} (critical={})",
        w.n_reads, w.recommended_min, w.is_critical);
}
```

Fires when `n_reads` is below `DiagnoseConfig::depth_warn_threshold` (default 10), and is
marked critical below `depth_critical_threshold` (default 5). Low depth means the coverage
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
    println!("interior support {:.2} at position {}", w.fraction, w.position);
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

Diagnostics only report this signal; they do not retry. When it fires, the caller should
rebuild with `band_width = 0` (unbanded) and compare.

### 5. Structural competing allele

```rust
if let Some(w) = &warnings.structural_competing {
    println!("{} competing structural site(s); weakest minority arm {} reads",
        w.site_count, w.min_arm_reads);
}
```

Fires (single-allele mode only) when `Consensus::bubble_sites` contains a structural
(length-changing) bubble whose minority arm clears the allele-frequency floor — a hint the
locus may be heterozygous and worth re-running with `consensus_multi`. Suppressed when
`DiagnoseConfig::is_allele_partition = true` (the caller is already in multi-allele mode).

## `DiagnoseConfig`

```rust
pub struct DiagnoseConfig {
    pub depth_warn_threshold: usize,         // default 10
    pub depth_critical_threshold: usize,     // default 5
    pub interior_support_threshold: f32,     // default 0.15
    pub boundary_margin: f32,                // default 0.20 (check the middle 60%)
    pub is_allele_partition: bool,           // default false; true = per-allele partition
    pub depth_allele_threshold: usize,       // default 15 (used when is_allele_partition)
    pub truncation_ratio_threshold: f32,     // default 0.60; set 0.0 to disable
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

// Is there a competing allele above 20%?
if let Some(site) = has_competing_allele(&result, 0.2) {
    println!("bubble at pos {}: {:?}", site.consensus_pos, site.arm_read_counts);
}

// Quick boolean for pipeline branching
if should_call_multiallele(&result, 0.2) {
    // re-run with consensus_multi
}

// Summary struct for logging
let conf = consensus_confidence(&result, 0.2);
println!("min_cov={} mean_cov={:.1} competing_allele={}",
    conf.min_cov, conf.mean_cov, conf.competing_allele);
```

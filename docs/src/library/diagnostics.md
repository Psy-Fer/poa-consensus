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
rebuild fully unbanded and compare. That needs both `band_width: 0` and
`adaptive_band: false`: with the default `adaptive_band: true`, `band_width: 0` alone still
uses the adaptive band.

The flank-anchored functions slice the output to the repeat, so it is naturally much
shorter than the reads. Ignore this check for them (the CLI skips it when flanks are
given), or set `truncation_ratio_threshold: 0.0`.

### 5. Structural competing allele

```rust
if let Some(w) = &warnings.structural_competing {
    println!("{} competing structural site(s); weakest minority arm {} reads",
        w.site_count, w.min_arm_reads);
}
```

Fires when `Consensus::bubble_sites` contains any structural (length-changing) site, a
hint that the locus may be heterozygous and worth re-running with `consensus_multi`.
`diagnose` applies no frequency threshold of its own; it relies on `bubble_sites` holding
only sites whose arms met `min_allele_freq` when the consensus was built. `site_count` is the
number of structural sites; `min_arm_reads` is the smallest non-zero arm read count across
all of them (an arm, not necessarily a minority arm). Suppressed when
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

## Per-base support quality

`diagnose` looks at the consensus as a whole. To see where along it the reads disagree,
every `Consensus` carries `support_quality`: one Phred-scaled value per base (capped at
60) saying how decisively the reads back that base over its strongest alternative.

```rust
use poa_consensus::analysis::low_support_regions;

// stretches where the reads do not decisively back the call
for r in low_support_regions(&result, 10) {
    println!("weak support at {}..{}", r.start, r.end);
}
```

How it is computed: each read votes from its route through the graph. Following the
consensus into a base agrees with it; any other route (a substitution, a skipped base,
an extra insertion) counts against the bases it bypasses. Reads taking the same
alternative route are pooled, and each base is scored only against its strongest
alternative. Sequencing errors scatter over many one-read routes and do not add up, but
a competing call, such as a group of reads one repeat unit shorter, concentrates on one
route and shows. With `a` reads for the base and `m` for its strongest alternative, the
value is `-10 log10 P(Binomial(a + m + 1, 1/2) > a)`:

| Reads for / against | Support quality |
|---|---|
| 10 / 0 | 33 |
| 20 / 0 | 60 (cap) |
| 7 / 3 | 9 |
| 6 / 4 | 6 |
| 5 / 5, or no spanning reads | 3 |

**Choosing a threshold.** Q10 is a good cut for "not solid". On simulated ONT and HiFi
GC-rich repeat ladders at 10 reads per haplotype, a minimum below Q10 inside the repeat
marked 98-99% of wrong repeat calls and 3-6% of right ones. Q20 marked no more wrong
calls but 18-34% of right ones. The score reflects how certain the reads' majority is,
so the same 2:1 split scores higher at higher depth.

**What it does not tell you.** It measures agreement, not truth: an error every read
shares (for example a systematic homopolymer undercall) still scores high. Inside a
repeat, where an indel sits is not determined by the sequence, so a disputed repeat
length shows as low support over part of the repeat rather than at one exact base.

The CLI writes it as the FASTQ quality string with `--fastq`.

## Analysis helpers

For more targeted analysis without the full `diagnose` call:

```rust
use poa_consensus::analysis::{
    min_coverage, low_coverage_regions, low_support_regions, has_competing_allele,
    should_call_multiallele, consensus_confidence,
};

// Minimum per-position coverage
let floor = min_coverage(&result);

// Positions below half-depth
let low = low_coverage_regions(&result, result.n_reads as u32 / 2);

// Positions the reads do not decisively back (support quality below Q10)
let weak = low_support_regions(&result, 10);

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

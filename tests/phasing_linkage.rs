//! Validates the linkage signal (step 1 of the multi-allele redesign): a real
//! second allele produces strongly-linked het sites; a single noisy allele does
//! not. This is the discriminator the folded-graph arm-span method lacked.

mod synth;

use poa_consensus::PoaConfig;
use poa_consensus::poa2::Poa;
use synth::{AlleleSpec, Read, Spec};

fn cfg() -> PoaConfig {
    PoaConfig {
        min_reads: 3,
        min_allele_freq: 0.2,
        ..PoaConfig::default()
    }
}

/// Build a poa2 graph from reads, seeding on the median-length read.
fn build(reads: &[&[u8]], config: &PoaConfig) -> Poa {
    let mut order: Vec<usize> = (0..reads.len()).collect();
    order.sort_by_key(|&i| reads[i].len());
    let med = order[order.len() / 2];
    let mut g = Poa::new(reads[med], config.clone());
    for (i, r) in reads.iter().enumerate() {
        if i != med {
            g.add_read(r);
        }
    }
    g
}

fn floor(n: usize) -> u32 {
    ((n as f64 * 0.2).ceil() as u32).max(3)
}

fn max_linkage_of(spec: &Spec) -> (usize, f64) {
    let reads = synth::generate(spec);
    let seqs = Read::seqs(&reads);
    let g = build(&seqs, &cfg());
    let m = g.phasing_matrix(floor(reads.len()));
    (m.n_sites(), m.max_linkage())
}

fn diploid(motif: &[u8], ua: usize, ub: usize, depth: usize, err: f64, seed: u64) -> Spec {
    Spec {
        left: synth::LEFT.to_vec(),
        right: synth::RIGHT.to_vec(),
        alleles: vec![
            AlleleSpec {
                mid: synth::repeat(motif, ua),
                n_reads: depth,
            },
            AlleleSpec {
                mid: synth::repeat(motif, ub),
                n_reads: depth,
            },
        ],
        sub: err,
        ins: err,
        del: err,
        partial_frac: 0.0,
        unit_jitter: true,
        unit_len: motif.len(),
        seed,
    }
}

fn single(motif: &[u8], units: usize, depth: usize, err: f64, seed: u64) -> Spec {
    Spec {
        left: synth::LEFT.to_vec(),
        right: synth::RIGHT.to_vec(),
        alleles: vec![AlleleSpec {
            mid: synth::repeat(motif, units),
            n_reads: depth,
        }],
        sub: err,
        ins: err,
        del: err,
        partial_frac: 0.0,
        unit_jitter: true,
        unit_len: motif.len(),
        seed,
    }
}

fn consistency_of(spec: &Spec) -> (f64, f64, usize) {
    let reads = synth::generate(spec);
    let seqs = Read::seqs(&reads);
    let g = build(&seqs, &cfg());
    let b = g.phasing_matrix(floor(reads.len())).dominant_bipartition();
    (b.consistency, b.minority_frac, b.n_informative_sites)
}

#[test]
fn real_diploid_produces_het_sites_and_a_balanced_split() {
    // Step 1 must at minimum surface the signal: a clean diploid yields het sites
    // and a roughly balanced dominant bipartition. (Confirmation that the split
    // is a REAL allele is the splitter's consensus-difference test, not here.)
    for (name, motif) in [
        ("CAG", &b"CAG"[..]),
        ("GAA", &b"GAA"[..]),
        ("CATCAT", &b"CATCAT"[..]),
    ] {
        for seed in [1u64, 2, 3] {
            let (sites, _) = max_linkage_of(&diploid(motif, 15, 30, 12, 0.03, seed));
            let (_, minority, informative) =
                consistency_of(&diploid(motif, 15, 30, 12, 0.03, seed));
            assert!(
                sites >= 2,
                "{name} seed{seed}: expected het sites, got {sites}"
            );
            assert!(informative >= 1, "{name} seed{seed}: no informative sites");
            assert!(
                minority >= 0.30,
                "{name} seed{seed}: split too skewed ({minority:.2}) for balanced het"
            );
        }
    }
}

#[test]
fn bipartition_consistency_beats_noise_on_average() {
    // Robust distributional claim (not a fragile per-instance threshold): over a
    // grid, the mean dominant-bipartition consistency of a real diploid clearly
    // exceeds that of a single noisy allele. Per-instance tails overlap — which is
    // exactly why the splitter needs the consensus-difference confirmation.
    let motifs: [&[u8]; 4] = [b"CAG", b"GAA", b"AAAAG", b"CATCAT"];
    let mut real = Vec::new();
    let mut noise = Vec::new();
    for m in motifs {
        for seed in [1u64, 2, 3] {
            real.push(consistency_of(&diploid(m, 15, 30, 12, 0.03, seed)).0);
            noise.push(consistency_of(&single(m, 30, 24, 0.06, seed)).0);
        }
    }
    let mean = |v: &Vec<f64>| v.iter().sum::<f64>() / v.len() as f64;
    let (rm, nm) = (mean(&real), mean(&noise));
    assert!(
        rm > nm + 0.08,
        "consistency does not separate on average: real {rm:.2} vs noise {nm:.2}"
    );
}

#[test]
#[ignore = "diagnostic; run with --ignored --nocapture"]
fn dbg_matrix_detail() {
    for (tag, spec) in [
        (
            "REAL diploid CAG15v30 d12 e3%",
            diploid(b"CAG", 15, 30, 12, 0.03, 1),
        ),
        (
            "NOISE single CAGx30 d24 e6%",
            single(b"CAG", 30, 24, 0.06, 1),
        ),
    ] {
        let reads = synth::generate(&spec);
        let seqs = Read::seqs(&reads);
        let g = build(&seqs, &cfg());
        let m = g.phasing_matrix(floor(reads.len()));
        eprintln!(
            "\n== {tag} ==  n_reads={} n_sites={}",
            m.n_reads,
            m.n_sites()
        );
        for (s, site) in m.sites.iter().enumerate() {
            let cov = m.covered(s);
            eprintln!(
                "  site {s}: node {} branches {:?} support {:?} covered {cov}",
                site.node, site.branches, site.branch_support
            );
        }
        // pairwise linkage + co-coverage
        let ns = m.n_sites();
        for a in 0..ns {
            for b in (a + 1)..ns {
                let shared = (0..m.n_reads)
                    .filter(|&r| {
                        m.calls[r][a] != poa_consensus::phasing::MISSING
                            && m.calls[r][b] != poa_consensus::phasing::MISSING
                    })
                    .count();
                let lk = m.linkage(a, b);
                if lk > 0.5 {
                    eprintln!("  linkage({a},{b}) = {lk:.2}  shared={shared}");
                }
            }
        }
    }
}

#[test]
#[ignore = "characterization; run with --ignored --nocapture"]
fn dbg_consistency_distribution() {
    let cons = |spec: &Spec| -> (f64, f64, usize) {
        let reads = synth::generate(spec);
        let seqs = Read::seqs(&reads);
        let g = build(&seqs, &cfg());
        let b = g.phasing_matrix(floor(reads.len())).dominant_bipartition();
        (b.consistency, b.minority_frac, b.n_informative_sites)
    };
    let motifs: [&[u8]; 4] = [b"CAG", b"GAA", b"AAAAG", b"CATCAT"];
    let mut real = Vec::new();
    for m in motifs {
        for (ua, ub) in [(15, 30), (20, 40), (30, 45)] {
            for e in [0.0, 0.03, 0.06] {
                for s in [1u64, 2, 3] {
                    real.push(cons(&diploid(m, ua, ub, 12, e, s)));
                }
            }
        }
    }
    let mut noise = Vec::new();
    for m in motifs {
        for u in [20usize, 30, 40] {
            for e in [0.03, 0.06, 0.09] {
                for s in [1u64, 2, 3] {
                    noise.push(cons(&single(m, u, 24, e, s)));
                }
            }
        }
    }
    let summ = |v: &Vec<(f64, f64, usize)>, name: &str| {
        let mut c: Vec<f64> = v.iter().map(|x| x.0).collect();
        c.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let q = |p: f64| c[((c.len() as f64 - 1.0) * p) as usize];
        eprintln!(
            "{name} (n={}): consistency min={:.2} p10={:.2} p50={:.2} p90={:.2} max={:.2}",
            v.len(),
            c[0],
            q(0.1),
            q(0.5),
            q(0.9),
            c[c.len() - 1]
        );
    };
    summ(&real, "REAL diploid");
    summ(&noise, "NOISE single");
    // how well does a threshold separate?
    for th in [0.78, 0.80, 0.82, 0.84] {
        let real_ok = real.iter().filter(|x| x.0 >= th).count();
        let noise_bad = noise.iter().filter(|x| x.0 >= th).count();
        eprintln!(
            "  threshold {th:.2}: real>=th {real_ok}/{}  noise>=th(false-pos) {noise_bad}/{}",
            real.len(),
            noise.len()
        );
    }
}

#[test]
#[ignore = "diagnostic; run with --ignored --nocapture"]
fn dbg_linkage_split_counts() {
    let run = |spec: &Spec| -> usize {
        let reads = synth::generate(spec);
        let seqs = Read::seqs(&reads);
        match poa_consensus::poa2::linkage_consensus_multi(&seqs, &cfg()) {
            Ok(a) => a.len(),
            Err(_) => 0,
        }
    };
    eprintln!("diploid length CAG (expect 2):");
    for (ua, ub) in [(15, 30), (20, 40), (30, 45)] {
        for e in [0.0, 0.03, 0.06] {
            let gots: Vec<usize> = [1u64, 2, 3]
                .iter()
                .map(|&s| run(&diploid(b"CAG", ua, ub, 12, e, s)))
                .collect();
            eprintln!("  {ua}v{ub} e{:.0}%: got {gots:?}", e * 100.0);
        }
    }
    eprintln!("single CAGx30 (expect 1):");
    for e in [0.0, 0.03, 0.06] {
        let gots: Vec<usize> = [1u64, 2, 3]
            .iter()
            .map(|&s| run(&single(b"CAG", 30, 24, e, s)))
            .collect();
        eprintln!("  e{:.0}%: got {gots:?}", e * 100.0);
    }
}

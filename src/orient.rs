use crate::types::Strand;
use std::borrow::Cow;

pub fn reverse_complement(read: &[u8]) -> Vec<u8> {
    read.iter().rev().map(|&b| complement(b)).collect()
}

#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        other => other,
    }
}

/// Determine the strand of `read` relative to `seed` using k-mer matching.
///
/// Counts shared k-mers between `read` and `seed` in the forward orientation,
/// then between `reverse_complement(read)` and `seed`. Returns `Forward` if the
/// forward count is greater or equal, `Reverse` otherwise.
///
/// Complexity: O(n) in the length of both sequences using a hash map.
pub fn orient_to_seed(read: &[u8], seed: &[u8], k: usize) -> Strand {
    if read.len() < k || seed.len() < k {
        return Strand::Forward;
    }

    use std::collections::HashMap;

    // Compare case-insensitively. `reverse_complement` uppercases as it
    // complements, so a soft-masked (lowercase) forward read would otherwise
    // score zero forward k-mer hits (case mismatch) and full hits on its
    // uppercased reverse complement, spuriously flipping a correctly-oriented
    // read. Uppercase both sequences up front so case never drives the call.
    let seed_u: Vec<u8> = seed.iter().map(|b| b.to_ascii_uppercase()).collect();
    let read_u: Vec<u8> = read.iter().map(|b| b.to_ascii_uppercase()).collect();

    let mut seed_kmers: HashMap<&[u8], u32> = HashMap::new();
    for w in seed_u.windows(k) {
        *seed_kmers.entry(w).or_insert(0) += 1;
    }

    let fwd_count: u32 = read_u
        .windows(k)
        .map(|w| seed_kmers.get(w).copied().unwrap_or(0))
        .sum();

    let rc = reverse_complement(&read_u);
    let rev_count: u32 = rc
        .windows(k)
        .map(|w| seed_kmers.get(w).copied().unwrap_or(0))
        .sum();

    if fwd_count >= rev_count {
        Strand::Forward
    } else {
        Strand::Reverse
    }
}

/// Orient all reads to match the strand of `reads[seed_idx]`.
///
/// Returns a `Vec<Cow<[u8]>>` where each element borrows the original read
/// if it is already in the forward orientation, or owns the reverse complement
/// if it was flipped.
pub fn auto_orient<'a>(reads: &'a [Vec<u8>], seed_idx: usize) -> Vec<Cow<'a, [u8]>> {
    let seed = &reads[seed_idx];
    let k = 8.min(seed.len().saturating_sub(1).max(1));

    reads
        .iter()
        .map(|read| match orient_to_seed(read, seed, k) {
            Strand::Forward => Cow::Borrowed(read.as_slice()),
            Strand::Reverse => Cow::Owned(reverse_complement(read)),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_roundtrip() {
        let s = b"ACGTACGGTTA";
        assert_eq!(reverse_complement(&reverse_complement(s)), s.to_vec());
    }

    #[test]
    fn revcomp_known() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT".to_vec());
        assert_eq!(reverse_complement(b"AACG"), b"CGTT".to_vec());
        // N passes through; case is normalized to uppercase.
        assert_eq!(reverse_complement(b"acgtN"), b"NACGT".to_vec());
    }

    #[test]
    fn forward_read_stays_forward() {
        let seed = b"ACGTACGTTTGGCCAATTGGCCAAGTAC";
        assert_eq!(orient_to_seed(seed, seed, 8), Strand::Forward);
    }

    #[test]
    fn reverse_read_detected() {
        let seed = b"ACGTACGTTTGGCCAATTGGCCAAGTAC";
        let rc = reverse_complement(seed);
        assert_eq!(orient_to_seed(&rc, seed, 8), Strand::Reverse);
    }

    #[test]
    fn lowercase_forward_read_is_not_flipped() {
        // Regression for the soft-masked-read bug: a lowercase forward read must
        // orient Forward, not be spuriously reverse-complemented because the
        // uppercased RC happened to match the uppercase seed's k-mers.
        let seed = b"ACGTACGTTTGGCCAATTGGCCAAGTAC";
        let lower: Vec<u8> = seed.iter().map(|b| b.to_ascii_lowercase()).collect();
        assert_eq!(orient_to_seed(&lower, seed, 8), Strand::Forward);
    }

    #[test]
    fn auto_orient_flips_only_reverse_reads() {
        let a = b"ACGTACGTTTGGCCAATTGGCCAAGTAC".to_vec();
        let rc = reverse_complement(&a);
        let reads = vec![a.clone(), rc.clone(), a.clone()];
        let oriented = auto_orient(&reads, 0);
        assert!(matches!(oriented[0], Cow::Borrowed(_)));
        assert!(matches!(oriented[1], Cow::Owned(_)));
        assert!(matches!(oriented[2], Cow::Borrowed(_)));
        // Every read now matches the seed orientation.
        assert_eq!(oriented[1].as_ref(), a.as_slice());
    }
}

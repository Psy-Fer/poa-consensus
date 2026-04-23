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

    let mut seed_kmers: HashMap<&[u8], u32> = HashMap::new();
    for w in seed.windows(k) {
        *seed_kmers.entry(w).or_insert(0) += 1;
    }

    let fwd_count: u32 = read
        .windows(k)
        .map(|w| seed_kmers.get(w).copied().unwrap_or(0))
        .sum();

    let rc = reverse_complement(read);
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

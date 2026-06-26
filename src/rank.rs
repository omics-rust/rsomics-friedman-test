//! Within-block tie-averaged ranking, matching SciPy `_rankdata` with
//! `method='average'` applied along the treatments axis of each block.
//!
//! For one block of `k` treatment values, SciPy stable-argsorts the values and
//! gives every member of a tie group the average of the ordinal ranks it spans
//! (`first + (count-1)/2`, 1-based). Friedman consumes only the per-treatment
//! rank sums across blocks and the multiset of tie-group sizes, so a block is
//! reduced directly into those two aggregates after a single sort of `0..k`.

/// Scratch buffers reused across blocks so per-block ranking allocates nothing.
pub struct RankScratch {
    order: Vec<usize>,
}

impl RankScratch {
    #[must_use]
    pub fn with_capacity(k: usize) -> Self {
        Self {
            order: Vec::with_capacity(k),
        }
    }

    /// Rank one block's treatment values in place, adding each treatment's
    /// tie-averaged rank to `rank_sums[j]` and accumulating `Σ(t³−t)` over the
    /// block's tie groups into `tie_term`.
    ///
    /// `block[j]` is treatment `j`'s value in this block; `rank_sums` and
    /// `block` have the same length `k`.
    pub fn rank_block(&mut self, block: &[f64], rank_sums: &mut [f64], tie_term: &mut f64) {
        let k = block.len();
        self.order.clear();
        self.order.extend(0..k);
        // Stable sort by value matches SciPy's stable argsort; ties keep input order.
        self.order.sort_by(|&a, &b| block[a].total_cmp(&block[b]));

        let mut i = 0;
        while i < k {
            let v = block[self.order[i]];
            let mut j = i + 1;
            while j < k && block[self.order[j]] == v {
                j += 1;
            }
            let count = j - i;
            let avg_rank = (i + 1) as f64 + (count as f64 - 1.0) / 2.0;
            for &t in &self.order[i..j] {
                rank_sums[t] += avg_rank;
            }
            if count > 1 {
                let t = count as f64;
                *tie_term += t * t * t - t;
            }
            i = j;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::RankScratch;

    fn rank(block: &[f64]) -> (Vec<f64>, f64) {
        let mut scratch = RankScratch::with_capacity(block.len());
        let mut sums = vec![0.0; block.len()];
        let mut tie = 0.0;
        scratch.rank_block(block, &mut sums, &mut tie);
        (sums, tie)
    }

    #[test]
    fn distinct_values_get_ordinal_ranks() {
        // ascending sort: 1.0->1, 5.0->2, 9.0->3 ; positions preserved
        let (sums, tie) = rank([9.0, 1.0, 5.0].as_slice());
        assert_eq!(sums, vec![3.0, 1.0, 2.0]);
        assert_eq!(tie, 0.0);
    }

    #[test]
    fn full_tie_block_all_get_mean_rank() {
        // four equal values share ranks 1..4 → all 2.5; tie term 4^3-4 = 60
        let (sums, tie) = rank([7.0, 7.0, 7.0, 7.0].as_slice());
        assert_eq!(sums, vec![2.5, 2.5, 2.5, 2.5]);
        assert!((tie - 60.0).abs() < 1e-12);
    }

    #[test]
    fn partial_tie_block() {
        // values [2,1,2,3]: one 1 (rank 1), two 2s (ranks 2,3 → 2.5), one 3 (rank 4)
        let (sums, tie) = rank([2.0, 1.0, 2.0, 3.0].as_slice());
        assert_eq!(sums, vec![2.5, 1.0, 2.5, 4.0]);
        // tie term: 2^3-2 = 6
        assert!((tie - 6.0).abs() < 1e-12);
    }
}

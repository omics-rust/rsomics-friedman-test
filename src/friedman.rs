//! Friedman chi-square test, matching `scipy.stats.friedmanchisquare`.
//!
//! For `k` treatments measured over `n` blocks, values are ranked within each
//! block (tie-averaged). With `Rⱼ` the rank sum of treatment `j` and
//! `ssbn = Σⱼ Rⱼ²`:
//!
//! ```text
//! c = 1 − Σ(t³−t) / (k·(k²−1)·n)
//! Q = (12 / (k·n·(k+1)) · ssbn − 3·n·(k+1)) / c
//! ```
//!
//! summed over all within-block tie groups of size `t`, with `p = chi2.sf(Q, k−1)`.
//! The arithmetic mirrors SciPy's expression order so `Q` is bit-identical.

use serde::Serialize;

use crate::igamc::chi2_sf;
use crate::rank::RankScratch;
use rsomics_common::{Result, RsomicsError};

/// Result of a Friedman chi-square test.
#[derive(Debug, Clone, Serialize)]
pub struct FriedmanResult {
    /// Tie-corrected Q statistic (chi-squared distributed under H₀).
    #[serde(rename = "Q")]
    pub q: f64,
    /// Degrees of freedom, k − 1.
    pub df: usize,
    /// Survival-function p-value, chi2.sf(Q, df).
    pub p: f64,
}

/// A block-by-treatment matrix: `rows` blocks, each a row of `k` treatment values.
pub struct Matrix {
    pub rows: Vec<Vec<f64>>,
    pub k: usize,
}

/// Compute the Friedman test over a block × treatment matrix. At least three
/// treatments and one block are required; every block must have `k` values and
/// contain no NaN (SciPy's default `nan_policy='propagate'` would yield NaN — we
/// fail loud instead).
pub fn friedman(matrix: &Matrix) -> Result<FriedmanResult> {
    let k = matrix.k;
    if k < 3 {
        return Err(RsomicsError::InvalidInput(format!(
            "Friedman test needs at least 3 treatments (columns), got {k}"
        )));
    }
    let n = matrix.rows.len();
    if n == 0 {
        return Err(RsomicsError::InvalidInput(
            "no blocks (rows) in input".into(),
        ));
    }

    let mut rank_sums = vec![0.0_f64; k];
    let mut tie_term = 0.0_f64;
    let mut scratch = RankScratch::with_capacity(k);

    for (bi, block) in matrix.rows.iter().enumerate() {
        if block.len() != k {
            return Err(RsomicsError::InvalidInput(format!(
                "block {} has {} values, expected {k}",
                bi + 1,
                block.len()
            )));
        }
        if block.iter().any(|v| v.is_nan()) {
            return Err(RsomicsError::InvalidInput(format!(
                "block {} contains NaN",
                bi + 1
            )));
        }
        scratch.rank_block(block, &mut rank_sums, &mut tie_term);
    }

    let kf = k as f64;
    let nf = n as f64;

    let c = 1.0 - tie_term / (kf * (kf * kf - 1.0) * nf);
    let ssbn: f64 = rank_sums.iter().map(|&r| r * r).sum();
    let q = (12.0 / (kf * nf * (kf + 1.0)) * ssbn - 3.0 * nf * (kf + 1.0)) / c;

    let df = k - 1;
    let p = chi2_sf(df as f64, q);

    Ok(FriedmanResult { q, df, p })
}

#[cfg(test)]
mod tests {
    use super::{Matrix, friedman};

    fn close(got: f64, want: f64, rel: f64) {
        let d = (got - want).abs() / want.abs().max(f64::MIN_POSITIVE);
        assert!(d <= rel, "got {got:e} want {want:e} rel {d:e} > {rel:e}");
    }

    fn mat(rows: Vec<Vec<f64>>) -> Matrix {
        let k = rows[0].len();
        Matrix { rows, k }
    }

    #[test]
    fn hand_case_no_ties() {
        // 3 treatments, 4 blocks, no within-row ties.
        // ranks per block:
        //  [1,2,3] -> (1,2,3)
        //  [2,1,3] -> (2,1,3)
        //  [3,1,2] -> (3,1,2)
        //  [1,3,2] -> (1,3,2)
        // Rⱼ = (7, 7, 10); ssbn = 49+49+100 = 198
        // Q = 12/(3*4*4)*198 - 3*4*4 = (12/48)*198 - 48 = 49.5 - 48 = 1.5 ; c=1
        let m = mat(vec![
            vec![1.0, 2.0, 3.0],
            vec![2.0, 1.0, 3.0],
            vec![3.0, 1.0, 2.0],
            vec![1.0, 3.0, 2.0],
        ]);
        let r = friedman(&m).unwrap();
        close(r.q, 1.5, 1e-12);
        assert_eq!(r.df, 2);
    }

    #[test]
    fn rejects_two_treatments() {
        let m = mat(vec![vec![1.0, 2.0], vec![2.0, 1.0]]);
        assert!(friedman(&m).is_err());
    }

    #[test]
    fn rejects_ragged_block() {
        let m = Matrix {
            rows: vec![vec![1.0, 2.0, 3.0], vec![1.0, 2.0]],
            k: 3,
        };
        assert!(friedman(&m).is_err());
    }

    #[test]
    fn rejects_nan() {
        let m = mat(vec![vec![1.0, f64::NAN, 3.0], vec![1.0, 2.0, 3.0]]);
        assert!(friedman(&m).is_err());
    }
}

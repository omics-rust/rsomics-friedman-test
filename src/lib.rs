//! Friedman chi-square test — `scipy.stats.friedmanchisquare` equivalent.
//!
//! Input is a TSV matrix: each row is a block (subject), each column a
//! treatment. Values are ranked within each block; the test reports `Q`, `df`,
//! and the chi-squared survival-function p-value.

mod friedman;
mod igamc;
mod rank;

use std::io::{BufRead, Write};

use rsomics_common::{Result, RsomicsError};

pub use friedman::{FriedmanResult, Matrix, friedman};

/// Parse a whitespace/tab-delimited block × treatment matrix. Each non-empty
/// line is one block; columns are treatments. All rows must share a width of at
/// least three, validated by `friedman`.
pub fn parse_matrix<R: BufRead>(reader: R) -> Result<Matrix> {
    let mut rows: Vec<Vec<f64>> = Vec::new();
    let mut width: Option<usize> = None;

    for (lineno, line) in reader.lines().enumerate() {
        let line = line.map_err(RsomicsError::Io)?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let mut row = Vec::new();
        for field in line.split_whitespace() {
            let v: f64 = field.parse().map_err(|_| {
                RsomicsError::InvalidInput(format!(
                    "line {}: value '{field}' is not a number",
                    lineno + 1
                ))
            })?;
            row.push(v);
        }
        match width {
            None => width = Some(row.len()),
            Some(w) if w != row.len() => {
                return Err(RsomicsError::InvalidInput(format!(
                    "line {}: row has {} columns, expected {w}",
                    lineno + 1,
                    row.len()
                )));
            }
            _ => {}
        }
        rows.push(row);
    }

    let k = width.ok_or_else(|| RsomicsError::InvalidInput("no rows in input".into()))?;
    Ok(Matrix { rows, k })
}

/// Run the test on a reader and write `Q<TAB>df<TAB>p` (no JSON; the framework
/// emits the JSON envelope when `--json` is set).
pub fn run_friedman<R: BufRead, W: Write>(reader: R, out: &mut W) -> Result<FriedmanResult> {
    let matrix = parse_matrix(reader)?;
    let result = friedman(&matrix)?;
    writeln!(out, "{}\t{}\t{}", result.q, result.df, result.p).map_err(RsomicsError::Io)?;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::{parse_matrix, run_friedman};

    #[test]
    fn parse_tab_matrix() {
        let input = "1\t2\t3\n3\t1\t2\n";
        let m = parse_matrix(input.as_bytes()).unwrap();
        assert_eq!(m.k, 3);
        assert_eq!(m.rows, vec![vec![1.0, 2.0, 3.0], vec![3.0, 1.0, 2.0]]);
    }

    #[test]
    fn parse_whitespace_matrix() {
        let input = "1  2   3\n3 1 2\n";
        let m = parse_matrix(input.as_bytes()).unwrap();
        assert_eq!(m.k, 3);
    }

    #[test]
    fn run_emits_three_fields() {
        let input = "1\t2\t3\n2\t1\t3\n3\t1\t2\n1\t3\t2\n";
        let mut out = Vec::new();
        let r = run_friedman(input.as_bytes(), &mut out).unwrap();
        let s = String::from_utf8(out).unwrap();
        let parts: Vec<&str> = s.trim().split('\t').collect();
        assert_eq!(parts.len(), 3);
        assert_eq!(parts[1], "2");
        assert!((r.q - 1.5).abs() < 1e-12);
    }

    #[test]
    fn rejects_ragged_rows() {
        assert!(parse_matrix("1\t2\t3\n1\t2\n".as_bytes()).is_err());
    }

    #[test]
    fn rejects_non_numeric() {
        assert!(parse_matrix("1\t2\tfoo\n".as_bytes()).is_err());
    }
}

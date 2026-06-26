# rsomics-friedman-test

The Friedman chi-square test for repeated measures — a value-exact, faster
replacement for `scipy.stats.friedmanchisquare`.

The Friedman test is the non-parametric analogue of the repeated-measures
one-way ANOVA: it tests whether several treatments measured on the same blocks
(subjects) come from the same distribution, ranking values within each block
rather than using raw values, so it is robust to non-normal data.

## Install

```bash
cargo install rsomics-friedman-test
```

## Usage

Input is a TSV/whitespace matrix: each row is one block (subject), each column
one treatment. At least three treatment columns are required. Values are ranked
within each block. Output is a single line `Q<TAB>df<TAB>p`.

```bash
rsomics-friedman-test data.tsv
# reads stdin when no path (or `-`) is given:
cat data.tsv | rsomics-friedman-test
```

Example input (4 blocks × 3 treatments):

```
1	2	3
2	1	3
3	1	2
1	3	2
```

Flags: `-t/--threads`, `-q/--quiet`, `--json` (machine-readable envelope).

## Method

For `k` treatments over `n` blocks, every block's `k` values are assigned
tie-averaged ranks. With `Rⱼ` the rank sum of treatment `j` and `ssbn = Σⱼ Rⱼ²`:

```
c = 1 − Σ(t³ − t) / (k·(k² − 1)·n)
Q = (12 / (k·n·(k + 1)) · ssbn − 3·n·(k + 1)) / c
```

summed over within-block tie groups of size `t`. The p-value is the chi-squared
survival function with `k − 1` degrees of freedom, `p = chi2.sf(Q, k−1)`.

This matches `scipy.stats.friedmanchisquare` exactly, including SciPy's
`_rankdata` within-block tie-averaging and tie correction. The chi-squared
survival function reproduces SciPy's `special.chdtrc(df, Q) = igamc(df/2, Q/2)`
via a port of the Cephes `igam`/`igamc` incomplete-gamma routines, so both the Q
statistic and the p-value match SciPy to full double precision (bit-identical on
the test goldens).

## Origin

This crate is an independent Rust reimplementation of
`scipy.stats.friedmanchisquare`. SciPy is BSD-3-Clause licensed, so reading and
citing its source is permitted; the implementation follows:

- M. Friedman, "The Use of Ranks to Avoid the Assumption of Normality Implicit
  in the Analysis of Variance", *Journal of the American Statistical
  Association*, 32(200), 675–701, 1937.
  <https://doi.org/10.1080/01621459.1937.10503522>
- SciPy `scipy.stats._stats_py.friedmanchisquare`, `_rankdata`, and the Cephes
  `chdtrc`/`igamc` incomplete-gamma routines.

Reference values for the compatibility tests were generated with SciPy 1.17.1
and frozen in `tests/golden/`; SciPy is not invoked at test time.

License: MIT OR Apache-2.0.
Upstream credit: SciPy <https://scipy.org> (BSD-3-Clause); Cephes Math Library
(Stephen L. Moshier, BSD-style).

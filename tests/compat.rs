//! Value-exact compatibility against `scipy.stats.friedmanchisquare`.
//!
//! Expected values were computed once with scipy 1.17.1 (`tests/golden/expected.json`)
//! and frozen here; no scipy is invoked at test time. The Q statistic is rational
//! rank arithmetic and must match to ~1e-12; the p-value goes through the Cephes
//! `igamc` port and must match to ~1e-12.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use rsomics_friedman_test::{friedman, parse_matrix};

struct Expected {
    name: &'static str,
    q: f64,
    df: usize,
    p: f64,
}

const EXPECTED: &[Expected] = &[
    Expected {
        name: "hand_ties",
        q: 0.960_000_000_000_013_6,
        df: 3,
        p: 0.810_929_469_086_764_9,
    },
    Expected {
        name: "no_ties",
        q: 0.75,
        df: 2,
        p: 0.687_289_278_790_972_1,
    },
    Expected {
        name: "doc_example",
        q: 0.228_571_428_571_427_76,
        df: 4,
        p: 0.993_946_268_476_305_2,
    },
    Expected {
        name: "heavy_ties",
        q: 2.036_809_815_950_916,
        df: 4,
        p: 0.728_988_473_840_503_5,
    },
];

fn rel(got: f64, want: f64) -> f64 {
    (got - want).abs() / want.abs().max(f64::MIN_POSITIVE)
}

#[test]
fn golden_matches_scipy() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    for e in EXPECTED {
        let path = dir.join(format!("{}.tsv", e.name));
        let f = File::open(&path).unwrap_or_else(|_| panic!("open {}", path.display()));
        let matrix = parse_matrix(BufReader::new(f)).expect("parse");
        let r = friedman(&matrix).expect("friedman");

        assert_eq!(r.df, e.df, "{}: df", e.name);
        let qr = rel(r.q, e.q);
        assert!(
            qr <= 1e-12,
            "{}: Q {} vs scipy {} (rel {qr:e})",
            e.name,
            r.q,
            e.q
        );
        let pr = rel(r.p, e.p);
        assert!(
            pr <= 1e-12,
            "{}: p {} vs scipy {} (rel {pr:e})",
            e.name,
            r.p,
            e.p
        );
    }
}

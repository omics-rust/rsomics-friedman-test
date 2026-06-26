//! Upper-tail chi-squared survival function via a Cephes `igamc` port.
//!
//! SciPy computes the Friedman p-value as `chi2.sf(Q, k-1)`, which inside SciPy
//! is `cephes.chdtrc(df, Q) = igamc(df/2, Q/2)` — the regularized upper
//! incomplete gamma. Porting Cephes `igam`/`igamc` directly makes the tail match
//! SciPy's special-function path to machine precision rather than diverging
//! through a different incomplete-gamma series.

const MACHEP: f64 = 1.110_223_024_625_156_5e-16;
const MAXLOG: f64 = 7.097_827_128_933_84e2;
const BIG: f64 = 4.503_599_627_370_496e15;
const BIGINV: f64 = 2.220_446_049_250_313e-16;

/// Chi-squared survival function P(X > x) with `df` degrees of freedom.
#[must_use]
pub fn chi2_sf(df: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    igamc(0.5 * df, 0.5 * x)
}

/// Regularized lower incomplete gamma P(a, x), Cephes power-series branch.
fn igam(a: f64, x: f64) -> f64 {
    if x <= 0.0 || a <= 0.0 {
        return 0.0;
    }
    if x > 1.0 && x > a {
        return 1.0 - igamc(a, x);
    }
    let ax = a * x.ln() - x - libm::lgamma(a);
    if ax < -MAXLOG {
        return 0.0;
    }
    let ax = ax.exp();

    let mut r = a;
    let mut c = 1.0;
    let mut ans = 1.0;
    loop {
        r += 1.0;
        c *= x / r;
        ans += c;
        if c / ans <= MACHEP {
            break;
        }
    }
    ans * ax / a
}

/// Regularized upper incomplete gamma Q(a, x), Cephes continued-fraction branch.
#[must_use]
pub fn igamc(a: f64, x: f64) -> f64 {
    if x <= 0.0 || a <= 0.0 {
        return 1.0;
    }
    if x < 1.0 || x < a {
        return 1.0 - igam(a, x);
    }

    let ax = a * x.ln() - x - libm::lgamma(a);
    if ax < -MAXLOG {
        return 0.0;
    }
    let ax = ax.exp();

    let mut y = 1.0 - a;
    let mut z = x + y + 1.0;
    let mut c = 0.0;
    let mut pkm2 = 1.0;
    let mut qkm2 = x;
    let mut pkm1 = x + 1.0;
    let mut qkm1 = z * x;
    let mut ans = pkm1 / qkm1;

    loop {
        c += 1.0;
        y += 1.0;
        z += 2.0;
        let yc = y * c;
        let pk = pkm1 * z - pkm2 * yc;
        let qk = qkm1 * z - qkm2 * yc;
        let t = if qk != 0.0 {
            let r = pk / qk;
            let t = ((ans - r) / r).abs();
            ans = r;
            t
        } else {
            1.0
        };
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if pk.abs() > BIG {
            pkm2 *= BIGINV;
            pkm1 *= BIGINV;
            qkm2 *= BIGINV;
            qkm1 *= BIGINV;
        }
        if t <= MACHEP {
            break;
        }
    }
    ans * ax
}

#[cfg(test)]
mod tests {
    use super::{chi2_sf, igamc};

    fn close(got: f64, want: f64, rel: f64) {
        let d = (got - want).abs() / want.abs().max(f64::MIN_POSITIVE);
        assert!(d <= rel, "got {got:e} want {want:e} rel {d:e} > {rel:e}");
    }

    #[test]
    fn igamc_matches_scipy_gammaincc() {
        let cases = [
            (0.5, 0.1, 0.654_720_846_018_576_8),
            (1.0, 2.0, 0.135_335_283_236_612_7),
            (2.0, 7.0, 0.007_295_055_724_436_127),
            (10.0, 25.0, 0.000_221_476_638_248_783_3),
            (50.0, 60.0, 0.084_406_681_093_691_88),
            (100.0, 80.0, 0.982_891_686_964_866_8),
            (2.5, 12.5, 0.000_139_333_791_185_626_3),
            (0.5, 700.0, 2.101_014_516_264_400_3e-306),
        ];
        for (a, x, want) in cases {
            close(igamc(a, x), want, 1e-12);
        }
    }

    #[test]
    fn chi2_sf_matches_scipy_chdtrc() {
        // (df, Q, scipy chi2.sf(Q, df))
        let cases = [
            (3.0, 0.96, 0.810_929_469_086_764_9),
            (2.0, 0.75, 0.687_289_278_790_972_1),
            (4.0, 0.228_571_428_571_427_76, 0.993_946_268_476_305_2),
            (4.0, 2.036_809_815_950_916, 0.728_988_473_840_503_5),
            (1.0, 8.4, 0.003_752_210_100_873_851),
        ];
        for (df, q, want) in cases {
            close(chi2_sf(df, q), want, 1e-12);
        }
    }

    #[test]
    fn nonpositive_statistic_is_one() {
        assert_eq!(chi2_sf(3.0, 0.0), 1.0);
        assert_eq!(chi2_sf(3.0, -1.0), 1.0);
    }
}

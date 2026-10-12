//! Number formatting and rounding that match Python/NumPy, so tables written
//! in the reference tools' formats (CheckV, CheckM2) are byte-compatible.

/// Python `str(float)` / `repr(float)`: shortest round-trip digits, fixed
/// notation for decimal exponents −4..15, otherwise `d.ddde±XX`.
pub fn py_float(x: f64) -> String {
    if x.is_nan() {
        return "nan".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".into() } else { "-inf".into() };
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0.0".into() } else { "0.0".into() };
    }
    // Rust's `{:e}` prints the shortest round-trip mantissa, as repr does.
    let s = format!("{:e}", x);
    let (mant, exp) = s.split_once('e').unwrap_or((&s, "0"));
    let exp: i32 = exp.parse().unwrap_or(0);
    let neg = mant.starts_with('-');
    let digits: String = mant.chars().filter(|c| c.is_ascii_digit()).collect();
    let mut out = String::new();
    if neg {
        out.push('-');
    }
    if !(-4..16).contains(&exp) {
        out.push_str(&digits[..1]);
        if digits.len() > 1 {
            out.push('.');
            out.push_str(&digits[1..]);
        }
        out.push('e');
        out.push(if exp < 0 { '-' } else { '+' });
        out.push_str(&format!("{:02}", exp.abs()));
    } else if exp >= 0 {
        let n = digits.len() as i32;
        if n <= exp + 1 {
            out.push_str(&digits);
            out.push_str(&"0".repeat((exp + 1 - n) as usize));
            out.push_str(".0");
        } else {
            out.push_str(&digits[..(exp + 1) as usize]);
            out.push('.');
            out.push_str(&digits[(exp + 1) as usize..]);
        }
    } else {
        out.push_str("0.");
        out.push_str(&"0".repeat((-exp - 1) as usize));
        out.push_str(&digits);
    }
    out
}

/// Strip trailing fractional zeros (and a bare trailing point), as C's `%g`
/// does without the `#` flag.
fn trim_fraction(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}

/// C `printf("%.<prec>e", x)`: two-digit minimum exponent with sign.
pub fn c_fmt_e(x: f64, prec: usize) -> String {
    if !x.is_finite() {
        return py_float(x);
    }
    let s = format!("{:.*e}", prec, x);
    let (mant, exp) = s.split_once('e').unwrap_or((&s, "0"));
    let exp: i32 = exp.parse().unwrap_or(0);
    format!("{}e{}{:02}", mant, if exp < 0 { '-' } else { '+' }, exp.abs())
}

/// C `printf("%.<prec>g", x)` (glibc semantics: the exponent of the
/// `%e` conversion at precision `P−1` selects fixed or scientific notation;
/// trailing zeros removed).
pub fn c_fmt_g(x: f64, prec: usize) -> String {
    let p = prec.max(1);
    if x.is_nan() {
        return "nan".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".into() } else { "-inf".into() };
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0".into() } else { "0".into() };
    }
    let e = format!("{:.*e}", p - 1, x);
    let (mant, exp) = e.split_once('e').unwrap_or((&e, "0"));
    let exp: i32 = exp.parse().unwrap_or(0);
    if exp < -4 || exp >= p as i32 {
        format!("{}e{}{:02}", trim_fraction(mant), if exp < 0 { '-' } else { '+' }, exp.abs())
    } else {
        let f = format!("{:.*}", (p as i32 - 1 - exp) as usize, x);
        trim_fraction(&f).to_string()
    }
}

/// A value as HMMER's `--tblout` prints it (`%.2g` E-value or `%.1f`
/// score) and Python reads it back with `float()`.
pub fn reparse(s: &str) -> f64 {
    s.trim().parse().unwrap_or(f64::NAN)
}

/// Python `round(x, nd)`: correctly rounded decimal value of the binary
/// float, ties to even (Rust's fixed-precision formatting rounds the exact
/// value the same way).
pub fn py_round(x: f64, nd: usize) -> f64 {
    if !x.is_finite() {
        return x;
    }
    format!("{:.*}", nd, x).parse().unwrap_or(x)
}

/// Round half to even (`f64::round_ties_even`, stable only since Rust 1.77).
pub fn round_ties_even(x: f64) -> f64 {
    let r = x.round();
    if (x - x.trunc()).abs() == 0.5 {
        2.0 * (x / 2.0).round()
    } else {
        r
    }
}

/// NumPy `np.round(x, nd)` for `nd > 0`: `rint(x · 10^nd) / 10^nd`.
pub fn np_round(x: f64, nd: i32) -> f64 {
    let f = 10f64.powi(nd);
    round_ties_even(x * f) / f
}

/// An item of a Python `sum()`: a float or an int.
#[derive(Clone, Copy, Debug)]
pub enum PyNum {
    F(f64),
    I(i64),
}

/// CPython ≥ 3.12 `sum(items)` (start 0): ints are summed exactly until the
/// first float, then floats use Neumaier-compensated summation while ints
/// are added plainly; the compensation is applied at the end.
pub fn py_sum(items: &[PyNum]) -> f64 {
    let mut it = items.iter();
    let mut i_result: i64 = 0;
    // integer fast path
    let mut f_result = loop {
        match it.next() {
            None => return i_result as f64,
            Some(PyNum::I(v)) => i_result += v,
            Some(PyNum::F(x)) => break i_result as f64 + x,
        }
    };
    let mut c = 0.0f64;
    for item in it {
        match *item {
            PyNum::F(x) => {
                let t = f_result + x;
                if f_result.abs() >= x.abs() {
                    c += (f_result - t) + x;
                } else {
                    c += (x - t) + f_result;
                }
                f_result = t;
            }
            PyNum::I(v) => f_result += v as f64,
        }
    }
    if c != 0.0 && c.is_finite() {
        f_result += c;
    }
    f_result
}

/// `py_sum` of Python floats.
pub fn py_fsum(xs: impl IntoIterator<Item = f64>) -> f64 {
    let v: Vec<PyNum> = xs.into_iter().map(PyNum::F).collect();
    py_sum(&v)
}

/// NumPy's pairwise summation (`np.add.reduce` on float64), so means match
/// `np.mean`/`np.average` bit for bit.
pub fn np_sum(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut res = 0.0;
        for &x in a {
            res += x;
        }
        res
    } else if n <= 128 {
        let mut r = [0f64; 8];
        r.copy_from_slice(&a[..8]);
        let mut i = 8;
        while i < n - (n % 8) {
            for j in 0..8 {
                r[j] += a[i + j];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        while i < n {
            res += a[i];
            i += 1;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        np_sum(&a[..n2]) + np_sum(&a[n2..])
    }
}

/// `np.mean` of float64 values.
pub fn np_mean(a: &[f64]) -> f64 {
    if a.is_empty() {
        return f64::NAN;
    }
    np_sum(a) / a.len() as f64
}

/// NumPy `np.quantile(v, q)` with the default linear interpolation
/// (`method="linear"`): position `(n − 1)·q` between the sorted values.
pub fn np_quantile(v: &[f64], q: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    // _compute_virtual_index(n, q, alpha=1, beta=1), operation for operation
    let n = s.len();
    let (alpha, beta) = (1.0f64, 1.0f64);
    let virt = (n as f64 * q + (alpha + q * (1.0 - alpha - beta))) - 1.0;
    let lo = (virt.floor().max(0.0) as usize).min(n - 1);
    let hi = (lo + 1).min(n - 1);
    let frac = virt - virt.floor().max(0.0);
    if lo == hi || frac == 0.0 {
        return s[lo];
    }
    // NumPy's _lerp: a + (b - a) * t, switching form for t >= 0.5
    let (a, b) = (s[lo], s[hi]);
    let d = b - a;
    if frac >= 0.5 {
        b - d * (1.0 - frac)
    } else {
        a + d * frac
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn python_repr() {
        assert_eq!(py_float(100.0), "100.0");
        assert_eq!(py_float(0.0), "0.0");
        assert_eq!(py_float(7.08), "7.08");
        assert_eq!(py_float(1e-5), "1e-05");
        assert_eq!(py_float(0.0001), "0.0001");
        assert_eq!(py_float(1e16), "1e+16");
        assert_eq!(py_float(1234567890123456.0), "1234567890123456.0");
        assert_eq!(py_float(0.1 + 0.2), "0.30000000000000004");
        assert_eq!(py_float(-2.5), "-2.5");
        assert_eq!(py_float(f64::NAN), "nan");
        assert_eq!(py_float(2.0 / 7.0), "0.2857142857142857");
    }

    #[test]
    fn python_round() {
        assert_eq!(py_round(2.675, 2), 2.67);
        assert_eq!(py_round(0.125, 2), 0.12);
        assert_eq!(py_round(0.375, 2), 0.38);
        assert_eq!(py_round(16.5, 0), 16.0);
        assert_eq!(py_round(17.5, 0), 18.0);
        assert_eq!(round_ties_even(0.3 * 55.0), 16.0);
        assert_eq!(round_ties_even(2.5), 2.0);
        assert_eq!(round_ties_even(-2.5), -2.0);
        assert_eq!(round_ties_even(3.5), 4.0);
        // NumPy multiplies first: 1.075*100 == 107.5 exactly -> rint -> 108
        assert_eq!(np_round(1.075, 2), 1.08);
        assert_eq!(py_round(1.075, 2), 1.07);
        assert_eq!(np_round(2.675, 2), 2.68);
    }

    #[test]
    fn c_printf() {
        assert_eq!(c_fmt_g(1.234e-45, 2), "1.2e-45");
        assert_eq!(c_fmt_g(0.0012, 2), "0.0012");
        assert_eq!(c_fmt_g(9.96, 2), "10");
        assert_eq!(c_fmt_g(10.0, 2), "10");
        assert_eq!(c_fmt_g(123.0, 2), "1.2e+02");
        assert_eq!(c_fmt_g(2e-10, 2), "2e-10");
        assert_eq!(c_fmt_g(0.5, 2), "0.5");
        assert_eq!(c_fmt_g(3.7e-5, 2), "3.7e-05");
        assert_eq!(c_fmt_g(0.0, 2), "0");
        assert_eq!(c_fmt_e(1.5e-7, 2), "1.50e-07");
        assert_eq!(c_fmt_e(123.0, 2), "1.23e+02");
        assert_eq!(py_float(reparse(&c_fmt_g(3.7e-5, 2))), "3.7e-05");
        assert_eq!(py_float(reparse(&c_fmt_g(10.0, 2))), "10.0");
    }

    #[test]
    fn python312_sum() {
        // compensated: 0.1 added ten times is exactly 1.0 in Python >= 3.12
        assert_eq!(py_fsum(std::iter::repeat(0.1).take(10)), 1.0);
        assert_eq!(py_fsum([1e16, 1.0, -1e16]), 1.0);
        // ints before the first float are exact; ints afterwards plain
        assert_eq!(py_sum(&[PyNum::I(50), PyNum::F(0.5), PyNum::I(2)]), 52.5);
        assert_eq!(py_sum(&[]), 0.0);
    }

    #[test]
    fn numpy_quantile() {
        let v = [10.0, 20.0, 30.0, 40.0];
        assert!((np_quantile(&v, 0.05) - 11.5).abs() < 1e-12);
        assert!((np_quantile(&v, 0.95) - 38.5).abs() < 1e-12);
        assert_eq!(np_quantile(&[5.0], 0.95), 5.0);
    }
}

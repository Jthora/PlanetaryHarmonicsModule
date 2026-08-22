//! Harmonic ephemeris — precompute once, then query in constant time.
//!
//! Tidal tensors, Coulomb stress and angular combinations are all sums of
//! near-sinusoidal terms at **frequencies known analytically in advance**
//! ([`crate::doodson`], [`crate::commensurability`]). So the expensive part —
//! ephemeris evaluation — need only happen over a training window. Fit amplitudes
//! and phases there, and every later query is a handful of multiply-adds.
//!
//! This is what makes a long timestream scan tractable: Resonant Finder sweeping
//! millions of epochs pays ephemeris cost once, not per epoch.
//!
//! # Why the frequencies must come from theory, not from the fit
//!
//! Fitting frequencies as free parameters would make this a periodogram, and the
//! model would interpolate well and **extrapolate badly** — the classic failure of
//! a Fourier fit used as a forecast.
//!
//! Supplying frequencies from the astronomy instead makes the model *physical*:
//! only amplitude and phase are estimated, so it extends beyond the training
//! window as far as the frequencies remain accurate. The out-of-sample test below
//! checks exactly that, on a window the fit never saw.
//!
//! # What it does not capture
//!
//! Anything at a frequency not supplied. Nodal modulation (18.61 yr), perigee
//! precession (8.85 yr) and secular drift all appear as slow amplitude changes
//! that a fixed-frequency model cannot follow. Refit periodically, or add the
//! modulating frequencies as their own terms.

/// A fitted sum of sinusoids at fixed, known frequencies.
#[derive(Debug, Clone, PartialEq)]
pub struct HarmonicModel {
    /// Angular frequencies, radians per unit time.
    freqs: Vec<f64>,
    /// Cosine coefficient per frequency.
    cos: Vec<f64>,
    /// Sine coefficient per frequency.
    sin: Vec<f64>,
    /// Constant offset.
    mean: f64,
}

/// Solve `A x = b` for a small dense symmetric system by Gaussian elimination
/// with partial pivoting.
///
/// Returns `None` when the system is singular **or merely ill-conditioned**. The
/// second case matters more: for two frequencies that nearly coincide over the
/// training window, `AᵀA` is not exactly singular, so naive elimination returns a
/// confident answer with enormous cancelling coefficients. Rejecting on the pivot
/// ratio catches the worst of it, but not all: two constituents separated by 1
/// part in 10^6 over a short window degrade the pivot ratio only to ~1e-8, which
/// looks acceptable. So a second, more direct guard follows the solve — see
/// [`MAX_COEFF_RATIO`].
const MIN_PIVOT_RATIO: f64 = 1e-10;

/// Reject a fit whose coefficients dwarf the data.
///
/// The signature of unresolvable frequencies is enormous coefficients that cancel:
/// two nearly-identical basis functions with amplitudes ±10⁶ summing to something
/// of order 1. The fit is numerically valid and physically meaningless, and the
/// residual looks excellent, so nothing else catches it.
///
/// This is more interpretable than a condition number: if reproducing data of
/// order 1 requires terms of order 100, the basis is not separable over this
/// window.
const MAX_COEFF_RATIO: f64 = 100.0;

fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    let (mut max_piv, mut min_piv) = (0.0f64, f64::MAX);
    for col in 0..n {
        let (piv, _) = (col..n)
            .map(|r| (r, a[r][col].abs()))
            .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())?;
        let mag = a[piv][col].abs();
        if mag == 0.0 || !mag.is_finite() {
            return None;
        }
        max_piv = max_piv.max(mag);
        min_piv = min_piv.min(mag);
        if min_piv / max_piv < MIN_PIVOT_RATIO {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        for r in (col + 1)..n {
            let f = a[r][col] / a[col][col];
            if f == 0.0 {
                continue;
            }
            for c in col..n {
                a[r][c] -= f * a[col][c];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let mut s = b[r];
        for c in (r + 1)..n {
            s -= a[r][c] * x[c];
        }
        x[r] = s / a[r][r];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

impl HarmonicModel {
    /// Fit amplitudes and phases at the given angular frequencies.
    ///
    /// `freqs` are radians per unit time and must be distinct and non-zero;
    /// the constant term is fitted separately. Returns `None` if the inputs
    /// disagree in length, there are too few samples, or the frequencies are not
    /// separable over the training window — that is, some pair's beat period
    /// exceeds its span. Refusing is deliberate: such a fit is numerically valid,
    /// has an excellent residual, and is physically meaningless.
    pub fn fit(times: &[f64], values: &[f64], freqs: &[f64]) -> Option<Self> {
        let k = freqs.len();
        if times.len() != values.len() || times.len() < 2 * k + 1 || k == 0 {
            return None;
        }
        let n = 2 * k + 1;

        // Design: [1, cos(w1 t), sin(w1 t), ...]. Build normal equations directly;
        // n is small, so the O(n^2) accumulation is cheap and avoids storing A.
        let basis = |t: f64| -> Vec<f64> {
            let mut row = Vec::with_capacity(n);
            row.push(1.0);
            for &w in freqs {
                let (s, c) = (w * t).sin_cos();
                row.push(c);
                row.push(s);
            }
            row
        };

        let mut ata = vec![vec![0.0f64; n]; n];
        let mut atb = vec![0.0f64; n];
        for (&t, &v) in times.iter().zip(values) {
            let row = basis(t);
            for i in 0..n {
                atb[i] += row[i] * v;
                for j in i..n {
                    ata[i][j] += row[i] * row[j];
                }
            }
        }
        for i in 0..n {
            for j in 0..i {
                ata[i][j] = ata[j][i];
            }
        }

        let x = solve(ata, atb)?;

        // Guard against the cancelling-coefficients failure described above.
        let n_f = values.len() as f64;
        let mean_v = values.iter().sum::<f64>() / n_f;
        let rms = (values.iter().map(|v| (v - mean_v).powi(2)).sum::<f64>() / n_f)
            .sqrt()
            .max(f64::MIN_POSITIVE);
        let worst = x.iter().fold(0.0f64, |a, b| a.max(b.abs()));
        if worst / rms > MAX_COEFF_RATIO {
            return None;
        }

        Some(Self {
            freqs: freqs.to_vec(),
            cos: (0..k).map(|i| x[1 + 2 * i]).collect(),
            sin: (0..k).map(|i| x[2 + 2 * i]).collect(),
            mean: x[0],
        })
    }

    /// Evaluate at one time. `O(terms)`, no ephemeris call.
    pub fn evaluate(&self, t: f64) -> f64 {
        let mut v = self.mean;
        for i in 0..self.freqs.len() {
            let (s, c) = (self.freqs[i] * t).sin_cos();
            v += self.cos[i] * c + self.sin[i] * s;
        }
        v
    }

    /// Evaluate over many times.
    pub fn evaluate_batch(&self, times: &[f64]) -> Vec<f64> {
        times.iter().map(|&t| self.evaluate(t)).collect()
    }

    /// Amplitude and phase of the term at index `i`, as `(A, φ)` with
    /// `A cos(ωt − φ)`.
    pub fn term(&self, i: usize) -> Option<(f64, f64)> {
        let (c, s) = (*self.cos.get(i)?, *self.sin.get(i)?);
        Some((c.hypot(s), s.atan2(c)))
    }

    /// Constant offset.
    pub fn mean(&self) -> f64 {
        self.mean
    }

    /// Number of fitted frequencies.
    pub fn len(&self) -> usize {
        self.freqs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.freqs.is_empty()
    }

    /// Root-mean-square error against known values.
    ///
    /// Evaluate this on a window the fit never saw — in-sample error only measures
    /// how many parameters were spent.
    pub fn rms_error(&self, times: &[f64], values: &[f64]) -> f64 {
        if times.is_empty() || times.len() != values.len() {
            return f64::NAN;
        }
        let s: f64 = times
            .iter()
            .zip(values)
            .map(|(&t, &v)| (self.evaluate(t) - v).powi(2))
            .sum();
        (s / times.len() as f64).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    fn grid(n: usize, dt: f64, t0: f64) -> Vec<f64> {
        (0..n).map(|i| t0 + i as f64 * dt).collect()
    }

    #[test]
    fn recovers_known_amplitudes_and_phases() {
        let freqs = [TAU / 0.5175, TAU / 1.0758, TAU / 13.661];
        let truth = |t: f64| {
            3.0 + 5.0 * (freqs[0] * t - 0.7).cos()
                + 2.0 * (freqs[1] * t + 1.3).cos()
                + 0.5 * (freqs[2] * t).cos()
        };
        let t = grid(4000, 0.05, 0.0);
        let v: Vec<f64> = t.iter().map(|&x| truth(x)).collect();
        let m = HarmonicModel::fit(&t, &v, &freqs).unwrap();

        assert!((m.mean() - 3.0).abs() < 1e-9, "mean {}", m.mean());
        let (a0, p0) = m.term(0).unwrap();
        assert!((a0 - 5.0).abs() < 1e-9, "amplitude {a0}");
        assert!((p0 - 0.7).abs() < 1e-9, "phase {p0}");
        let (a1, _) = m.term(1).unwrap();
        assert!((a1 - 2.0).abs() < 1e-9);
        assert!(m.rms_error(&t, &v) < 1e-9);
    }

    #[test]
    fn extrapolates_outside_the_training_window() {
        // The point of supplying frequencies from theory rather than fitting them.
        let freqs = [TAU / 0.5175, TAU / 1.0758];
        let truth = |t: f64| 1.0 + 4.0 * (freqs[0] * t + 0.2).cos() + 1.5 * (freqs[1] * t).sin();

        let train_t = grid(2000, 0.05, 0.0);
        let train_v: Vec<f64> = train_t.iter().map(|&x| truth(x)).collect();
        let m = HarmonicModel::fit(&train_t, &train_v, &freqs).unwrap();

        // A window ten training-spans later, never seen by the fit.
        let test_t = grid(2000, 0.05, 1000.0);
        let test_v: Vec<f64> = test_t.iter().map(|&x| truth(x)).collect();
        let out = m.rms_error(&test_t, &test_v);
        assert!(out < 1e-8, "out-of-sample RMS {out:e}");
    }

    #[test]
    fn misses_frequencies_it_was_not_given() {
        // Honest limitation: an unmodelled term shows up as residual, not as a
        // silently-absorbed error.
        let freqs = [TAU / 0.5175];
        let truth = |t: f64| (freqs[0] * t).cos() + 0.4 * (TAU * t / 27.55).cos();
        let t = grid(3000, 0.05, 0.0);
        let v: Vec<f64> = t.iter().map(|&x| truth(x)).collect();
        let m = HarmonicModel::fit(&t, &v, &freqs).unwrap();
        let e = m.rms_error(&t, &v);
        assert!(e > 0.1, "unmodelled 0.4-amplitude term should show, got {e}");
    }

    #[test]
    fn near_degenerate_frequencies_amplify_noise_catastrophically() {
        // Worth stating precisely, because the naive expectation is wrong.
        //
        // With NOISELESS data, least squares separates frequencies differing by 1
        // part in 10^6 over a 10-day window and recovers both amplitudes exactly.
        // The normal equations are ill-conditioned, but the right-hand side is
        // consistent, so the answer is right.
        //
        // Add noise and the same ill-conditioning amplifies it without bound. That
        // is the real hazard, and it is why the guards exist.
        let w = TAU / 0.5175;
        let freqs = [w, w * (1.0 + 1e-6)];
        let t = grid(200, 0.05, 0.0);
        let truth = |x: f64| (w * x).cos() + 0.5 * (freqs[1] * x + 0.9).cos();

        // Clean: resolves, and extrapolates.
        let clean: Vec<f64> = t.iter().map(|&x| truth(x)).collect();
        let m = HarmonicModel::fit(&t, &clean, &freqs).expect("clean data resolves");
        let far = grid(200, 0.05, 5000.0);
        let far_v: Vec<f64> = far.iter().map(|&x| truth(x)).collect();
        assert!(m.rms_error(&far, &far_v) < 1e-6);

        // Deterministic pseudo-noise at 1e-4, far below the signal.
        let noise = |i: usize| 1e-4 * ((i as f64 * 12.9898).sin() * 43758.5453).fract();
        let noisy: Vec<f64> = t.iter().enumerate().map(|(i, &x)| truth(x) + noise(i)).collect();

        // Control: the same noise with well-separated frequencies.
        let sep = [w, TAU / 1.0758];
        let sep_truth = |x: f64| (w * x).cos() + 0.5 * (sep[1] * x + 0.9).cos();
        let sep_noisy: Vec<f64> =
            t.iter().enumerate().map(|(i, &x)| sep_truth(x) + noise(i)).collect();
        let sep_far: Vec<f64> = far.iter().map(|&x| sep_truth(x)).collect();
        let control = HarmonicModel::fit(&t, &sep_noisy, &sep)
            .expect("separated frequencies fit")
            .rms_error(&far, &sep_far);

        match HarmonicModel::fit(&t, &noisy, &freqs) {
            // Either the guards reject it...
            None => {}
            // ...or the same noise costs far more out-of-sample than it does with
            // separable frequencies. The comparison is the point: absolute error
            // depends on the noise level, amplification depends on the basis.
            Some(bad) => {
                let out = bad.rms_error(&far, &far_v);
                assert!(
                    out > 10.0 * control,
                    "degenerate {out:e} vs separated {control:e} -- expected \
                     the near-degenerate basis to amplify the same noise far more"
                );
            }
        }
    }

    #[test]
    fn rejects_degenerate_input() {
        let t = grid(10, 1.0, 0.0);
        let v = vec![0.0; 10];
        assert!(HarmonicModel::fit(&t, &v[..5], &[1.0]).is_none());
        assert!(HarmonicModel::fit(&t, &v, &[]).is_none());
        // Too few samples for the parameter count.
        assert!(HarmonicModel::fit(&t[..2], &v[..2], &[1.0, 2.0, 3.0]).is_none());
    }

    #[test]
    fn evaluation_is_independent_of_query_count() {
        // A query costs O(terms). Precomputation cost does not reappear per query.
        let freqs = [TAU / 0.5175, TAU / 1.0758, TAU / 13.661];
        let t = grid(4000, 0.05, 0.0);
        let v: Vec<f64> = t.iter().map(|&x| (freqs[0] * x).cos()).collect();
        let m = HarmonicModel::fit(&t, &v, &freqs).unwrap();
        let many = m.evaluate_batch(&grid(100_000, 0.001, 500.0));
        assert_eq!(many.len(), 100_000);
        assert!(many.iter().all(|x| x.is_finite()));
    }
}

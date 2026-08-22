//! Complex demodulation — isolating one constituent from a quasi-periodic series.
//!
//! To measure a **transfer function** we need response at several frequencies
//! separately. Folding on a trial period cannot do that: it mixes every
//! constituent whose period is near the fold, and on a detection-limited
//! catalogue it measures the detector (see `docs/07`, C1/D1).
//!
//! Complex demodulation separates cleanly. For a target period `P`, write
//! `θ = 2πt/P` and form
//!
//! ```text
//! z(t) = v(t) · e^{−iθ}
//! ```
//!
//! For `v = A cos(θ + φ₀)` this gives `z = (A/2)(e^{iφ₀} + e^{−i(2θ+φ₀)})`. The
//! second term oscillates at `2θ` and averages away, so a moving average over at
//! least one period leaves `(A/2)e^{iφ₀}`. Hence
//!
//! ```text
//! amplitude(t) = 2|z̄(t)|        band phase(t) = θ + arg z̄(t)
//! ```
//!
//! Both are *instantaneous*, so amplitude modulation (spring–neap, perigee) and
//! slow phase drift are tracked rather than smeared.
//!
//! # Choosing the window
//!
//! The window sets resolution. To separate two constituents the window must
//! exceed their **beat period**, `1/|1/P₁ − 1/P₂|`. M2 (0.5175 d) and S2 (0.5 d)
//! beat at **14.77 d**, so separating them needs a window well above that —
//! 60 d is a reasonable default for the semidiurnal band.
//!
//! Too long a window and genuine amplitude modulation is smoothed away; too short
//! and neighbouring constituents leak in. [`Band::leakage_note`] reports the beat
//! period against a companion constituent so the trade-off stays explicit.

/// One demodulated constituent band.
#[derive(Debug, Clone)]
pub struct Band {
    period: f64,
    window: f64,
    times: Vec<f64>,
    /// Slowly varying phase offset `arg z̄`, unwrapped so it interpolates safely.
    offset: Vec<f64>,
    /// Instantaneous amplitude `2|z̄|`.
    amplitude: Vec<f64>,
    /// Samples trimmed from each end, where the moving average is incomplete.
    trim: usize,
}

/// Demodulate `values` sampled at `times` at the given `period`.
///
/// `window` is the moving-average length, in the same units as `times`. Returns
/// `None` if the inputs disagree in length, or the window covers fewer than two
/// samples, or trimming would leave nothing.
pub fn demodulate(times: &[f64], values: &[f64], period: f64, window: f64) -> Option<Band> {
    if times.len() != values.len() || times.len() < 8 || period <= 0.0 || window <= 0.0 {
        return None;
    }
    let dt = (times[times.len() - 1] - times[0]) / (times.len() - 1) as f64;
    let half = ((window / dt) / 2.0).round() as usize;
    if half < 1 || 2 * half + 1 >= times.len() {
        return None;
    }

    let tau = std::f64::consts::TAU;
    let n = times.len();

    // Prefix sums give an O(n) moving average.
    let (mut pre_re, mut pre_im) = (vec![0.0f64; n + 1], vec![0.0f64; n + 1]);
    for k in 0..n {
        let (s, c) = (tau * times[k] / period).sin_cos();
        // z = v * e^{-i theta}
        pre_re[k + 1] = pre_re[k] + values[k] * c;
        pre_im[k + 1] = pre_im[k] + values[k] * (-s);
    }

    let mut offset = Vec::with_capacity(n);
    let mut amplitude = Vec::with_capacity(n);
    for k in 0..n {
        let lo = k.saturating_sub(half);
        let hi = (k + half + 1).min(n);
        let m = (hi - lo) as f64;
        let re = (pre_re[hi] - pre_re[lo]) / m;
        let im = (pre_im[hi] - pre_im[lo]) / m;
        amplitude.push(2.0 * (re * re + im * im).sqrt());
        offset.push(im.atan2(re));
    }

    // Unwrap so linear interpolation of the offset is safe across ±π.
    let pi = std::f64::consts::PI;
    for k in 1..n {
        let mut d = offset[k] - offset[k - 1];
        while d > pi {
            offset[k] -= tau;
            d -= tau;
        }
        while d < -pi {
            offset[k] += tau;
            d += tau;
        }
    }

    Some(Band {
        period,
        window,
        times: times.to_vec(),
        offset,
        amplitude,
        trim: half,
    })
}

impl Band {
    pub fn period(&self) -> f64 {
        self.period
    }

    pub fn window(&self) -> f64 {
        self.window
    }

    /// First and last times where the moving average is complete.
    pub fn valid_span(&self) -> (f64, f64) {
        (
            self.times[self.trim],
            self.times[self.times.len() - 1 - self.trim],
        )
    }

    fn interp(&self, series: &[f64], t: f64) -> Option<f64> {
        let (lo, hi) = self.valid_span();
        if t < lo || t > hi {
            return None;
        }
        let i = self.times.partition_point(|&x| x <= t).max(1) - 1;
        let j = (i + 1).min(self.times.len() - 1);
        if j == i {
            return Some(series[i]);
        }
        let f = (t - self.times[i]) / (self.times[j] - self.times[i]);
        Some(series[i] + f * (series[j] - series[i]))
    }

    /// Phase within this band at `t`, in `[0, 2π)`, with 0 at a band maximum.
    pub fn phase_at(&self, t: f64) -> Option<f64> {
        let off = self.interp(&self.offset, t)?;
        let p = std::f64::consts::TAU * t / self.period + off;
        let tau = std::f64::consts::TAU;
        let r = p % tau;
        Some(if r < 0.0 { r + tau } else { r })
    }

    /// Instantaneous band amplitude at `t`.
    pub fn amplitude_at(&self, t: f64) -> Option<f64> {
        self.interp(&self.amplitude, t)
    }

    /// Mean amplitude over the valid span — the constituent's strength.
    pub fn mean_amplitude(&self) -> f64 {
        let a = &self.amplitude[self.trim..self.amplitude.len() - self.trim];
        a.iter().sum::<f64>() / a.len() as f64
    }

    /// Phases for a set of event times, with a count of those outside the span.
    pub fn phases(&self, events: &[f64]) -> (Vec<f64>, usize) {
        let mut out = Vec::with_capacity(events.len());
        let mut dropped = 0;
        for &t in events {
            match self.phase_at(t) {
                Some(p) => out.push(p),
                None => dropped += 1,
            }
        }
        (out, dropped)
    }

    /// Beat period against a companion constituent — the window must exceed it
    /// for the two to be separated.
    pub fn leakage_note(&self, other_period: f64) -> f64 {
        let d = (1.0 / self.period - 1.0 / other_period).abs();
        if d <= 0.0 {
            f64::INFINITY
        } else {
            1.0 / d
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    fn grid(n: usize, dt: f64) -> Vec<f64> {
        (0..n).map(|i| i as f64 * dt).collect()
    }

    #[test]
    fn recovers_amplitude_and_phase_of_a_pure_tone() {
        let t = grid(20_000, 0.01);
        let (amp, period, phi0) = (3.0, 0.5175, 0.9);
        let v: Vec<f64> = t.iter().map(|&x| amp * (TAU * x / period + phi0).cos()).collect();
        let b = demodulate(&t, &v, period, 20.0).unwrap();

        assert!((b.mean_amplitude() - amp).abs() < 0.02, "{}", b.mean_amplitude());
        // Maxima are at TAU*t/P + phi0 = 2*pi*m, i.e. t = (m - phi0/TAU)*P.
        // Pick one near t=100, comfortably inside the trimmed span (10 to 190).
        let t_max = (193.0 - phi0 / TAU) * period;
        let (lo, hi) = b.valid_span();
        assert!(t_max > lo && t_max < hi, "t_max {t_max} outside {lo}..{hi}");
        let p = b.phase_at(t_max).unwrap();
        assert!(p < 0.05 || TAU - p < 0.05, "phase at maximum was {p}");
    }

    #[test]
    fn separates_constituents_when_the_window_exceeds_the_beat() {
        // M2 and S2: beat period 14.77 d.
        let (m2, s2) = (0.5175, 0.5);
        let t = grid(200_000, 0.005);
        let v: Vec<f64> = t
            .iter()
            .map(|&x| 1.0 * (TAU * x / m2).cos() + 4.0 * (TAU * x / s2).cos())
            .collect();

        let b = demodulate(&t, &v, m2, 90.0).unwrap();
        assert!(
            (b.leakage_note(s2) - 14.77).abs() < 0.1,
            "beat {}",
            b.leakage_note(s2)
        );
        // Despite S2 being 4x larger, demodulating at M2 must recover ~1.
        assert!(
            (b.mean_amplitude() - 1.0).abs() < 0.15,
            "recovered {} from a 4x larger neighbour",
            b.mean_amplitude()
        );
    }

    #[test]
    fn short_window_leaks_the_neighbour() {
        let (m2, s2) = (0.5175, 0.5);
        let t = grid(200_000, 0.005);
        let v: Vec<f64> = t
            .iter()
            .map(|&x| 1.0 * (TAU * x / m2).cos() + 4.0 * (TAU * x / s2).cos())
            .collect();
        // 2 d is far below the 14.77 d beat, so S2 must contaminate.
        let b = demodulate(&t, &v, m2, 2.0).unwrap();
        assert!(
            b.mean_amplitude() > 1.5,
            "expected leakage, got {}",
            b.mean_amplitude()
        );
    }

    #[test]
    fn tracks_an_amplitude_envelope() {
        let period = 0.5;
        let t = grid(100_000, 0.01);
        let v: Vec<f64> = t
            .iter()
            .map(|&x| (1.0 + 0.8 * (TAU * x / 200.0).sin()) * (TAU * x / period).cos())
            .collect();
        let b = demodulate(&t, &v, period, 20.0).unwrap();
        let (lo, hi) = b.valid_span();
        let mut mn = f64::MAX;
        let mut mx = f64::MIN;
        let mut x = lo;
        while x < hi {
            let a = b.amplitude_at(x).unwrap();
            mn = mn.min(a);
            mx = mx.max(a);
            x += 1.0;
        }
        assert!(mn < 0.35, "envelope min {mn}");
        assert!(mx > 1.7, "envelope max {mx}");
    }

    #[test]
    fn rejects_bad_input() {
        let t = grid(100, 0.1);
        assert!(demodulate(&t, &t[..50], 1.0, 1.0).is_none());
        assert!(demodulate(&t, &t, 1.0, 1e6).is_none());
        assert!(demodulate(&t, &t, -1.0, 1.0).is_none());
    }
}

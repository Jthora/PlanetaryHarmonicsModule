//! Tidal phase from a sampled quasi-periodic forcing.
//!
//! # Why not fold on a trial period
//!
//! Folding event times modulo a fixed period assumes the forcing is a pure
//! sinusoid. Real tidal forcing is **quasi-periodic** — many constituents beating,
//! frequency-modulated by the nodal and perigee cycles — so a fixed-period fold
//! smears the phase and, worse, makes the time-shift null degenerate: a global
//! shift then rotates the phase cluster without changing its concentration, so
//! `D²ₙ` is exactly invariant and the null has zero power
//! (see [`crate::stats`]).
//!
//! # The method
//!
//! The standard construction in the tidal-triggering literature (Tanaka et al.
//! 2002 and successors): sample the forcing finely, locate its **successive
//! maxima**, and assign each event a phase by interpolating between the maxima
//! that bracket it.
//!
//! ```text
//! phase(t) = 2π · (t − t_prev_max) / (t_next_max − t_prev_max)
//! ```
//!
//! Phase 0 is a maximum of the forcing. Because the bracketing maxima are the
//! *actual* ones, cycle-to-cycle variation in period is absorbed exactly — which
//! is what restores power to the time-shift null.

/// A sampled scalar forcing series.
///
/// Times and values must be the same length, with times strictly increasing.
#[derive(Debug, Clone)]
pub struct Forcing {
    times: Vec<f64>,
    values: Vec<f64>,
    maxima: Vec<f64>,
}

impl Forcing {
    /// Build from evenly or unevenly sampled sequences and locate the maxima.
    ///
    /// Peak times are refined by parabolic interpolation through each local
    /// maximum and its neighbours, giving sub-sample accuracy. Returns `None` if
    /// the inputs disagree in length or hold fewer than three samples.
    pub fn new(times: Vec<f64>, values: Vec<f64>) -> Option<Self> {
        if times.len() != values.len() || times.len() < 3 {
            return None;
        }
        let mut maxima = Vec::new();
        for i in 1..values.len() - 1 {
            let (a, b, c) = (values[i - 1], values[i], values[i + 1]);
            if b > a && b > c {
                // Parabolic vertex offset, in samples, clamped to the cell.
                let denom = a - 2.0 * b + c;
                let off = if denom.abs() > f64::EPSILON {
                    (0.5 * (a - c) / denom).clamp(-0.5, 0.5)
                } else {
                    0.0
                };
                let dt = if off >= 0.0 {
                    times[i + 1] - times[i]
                } else {
                    times[i] - times[i - 1]
                };
                maxima.push(times[i] + off * dt);
            }
        }
        Some(Self {
            times,
            values,
            maxima,
        })
    }

    /// Interpolated times of successive forcing maxima.
    pub fn maxima(&self) -> &[f64] {
        &self.maxima
    }

    /// Mean interval between successive maxima — the forcing's dominant period.
    pub fn mean_period(&self) -> Option<f64> {
        if self.maxima.len() < 2 {
            return None;
        }
        let n = self.maxima.len();
        Some((self.maxima[n - 1] - self.maxima[0]) / (n - 1) as f64)
    }

    /// Phase of the forcing at `t`, in `[0, 2π)`, with 0 at a maximum.
    ///
    /// Returns `None` when `t` falls outside the bracketing maxima — events before
    /// the first maximum or after the last cannot be assigned a phase and must be
    /// dropped rather than guessed.
    pub fn phase_at(&self, t: f64) -> Option<f64> {
        let m = &self.maxima;
        if m.len() < 2 || t < m[0] || t > m[m.len() - 1] {
            return None;
        }
        // Index of the last maximum at or before t.
        let i = match m.binary_search_by(|p| p.partial_cmp(&t).unwrap()) {
            Ok(i) => i.min(m.len() - 2),
            Err(i) => i - 1,
        };
        let (lo, hi) = (m[i], m[i + 1]);
        if hi <= lo {
            return None;
        }
        Some(std::f64::consts::TAU * (t - lo) / (hi - lo))
    }

    /// Phases for a set of event times, dropping those outside the maxima span.
    ///
    /// Returns the phases and the number dropped.
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

    /// Peak-to-trough amplitude of the forcing cycle containing `t`.
    ///
    /// The cycle is the interval between the maxima bracketing `t`, and the
    /// amplitude is the largest minus the smallest sampled value within it.
    ///
    /// This is the handle for testing the **amplitude law** `R̃/r = Δτ/(aσ̄)`:
    /// a detection artifact does not care how strong the tide is, so any
    /// dependence of response on cycle amplitude is physical rather than
    /// instrumental.
    pub fn cycle_amplitude_at(&self, t: f64) -> Option<f64> {
        let m = &self.maxima;
        if m.len() < 2 || t < m[0] || t > m[m.len() - 1] {
            return None;
        }
        let i = match m.binary_search_by(|p| p.partial_cmp(&t).unwrap()) {
            Ok(i) => i.min(m.len() - 2),
            Err(i) => i - 1,
        };
        let (lo, hi) = (m[i], m[i + 1]);
        let a = self.times.partition_point(|&x| x < lo);
        let b = self.times.partition_point(|&x| x <= hi);
        if b <= a {
            return None;
        }
        let slice = &self.values[a..b];
        let mx = slice.iter().cloned().fold(f64::MIN, f64::max);
        let mn = slice.iter().cloned().fold(f64::MAX, f64::min);
        Some(mx - mn)
    }

    /// Sampled times.
    pub fn times(&self) -> &[f64] {
        &self.times
    }

    /// Sampled values.
    pub fn values(&self) -> &[f64] {
        &self.values
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    fn sample<F: Fn(f64) -> f64>(f: F, n: usize, dt: f64) -> Forcing {
        let times: Vec<f64> = (0..n).map(|i| i as f64 * dt).collect();
        let values: Vec<f64> = times.iter().map(|&t| f(t)).collect();
        Forcing::new(times, values).unwrap()
    }

    #[test]
    fn finds_maxima_of_a_sinusoid() {
        let period = 10.0;
        // span = n*dt = 500 units = 50 cycles of a period-10 forcing
        let f = sample(|t| (TAU * t / period).sin(), 10_000, 0.05);
        let m = f.maxima();
        assert!(m.len() > 40, "got {} maxima", m.len());
        // First maximum of sin is at a quarter period.
        assert!((m[0] - 2.5).abs() < 0.02, "first max at {}", m[0]);
        assert!((f.mean_period().unwrap() - period).abs() < 0.02);
    }

    #[test]
    fn phase_is_zero_at_a_maximum() {
        let f = sample(|t| (TAU * t / 10.0).sin(), 10_000, 0.05);
        let m = f.maxima()[3];
        let p = f.phase_at(m).unwrap();
        assert!(p < 1e-9 || (TAU - p) < 1e-9, "phase at max was {p}");
    }

    #[test]
    fn phase_is_half_way_between_maxima() {
        let f = sample(|t| (TAU * t / 10.0).sin(), 10_000, 0.05);
        let (a, b) = (f.maxima()[3], f.maxima()[4]);
        let p = f.phase_at(0.5 * (a + b)).unwrap();
        assert!((p - std::f64::consts::PI).abs() < 1e-6, "got {p}");
    }

    #[test]
    fn events_outside_the_maxima_span_are_dropped() {
        let f = sample(|t| (TAU * t / 10.0).sin(), 10_000, 0.05);
        let first = f.maxima()[0];
        let (phases, dropped) = f.phases(&[first - 1.0, first + 1.0, 1e9]);
        assert_eq!(dropped, 2);
        assert_eq!(phases.len(), 1);
    }

    #[test]
    fn absorbs_a_varying_period() {
        // Frequency-modulated forcing: successive periods genuinely differ.
        let f = sample(
            |t| (TAU * (t / 10.0 + 0.8 * (TAU * t / 200.0).sin())).sin(),
            40_000,
            0.05,
        );
        let m = f.maxima();
        let intervals: Vec<f64> = m.windows(2).map(|w| w[1] - w[0]).collect();
        let (lo, hi) = intervals
            .iter()
            .fold((f64::MAX, f64::MIN), |(l, h), &x| (l.min(x), h.max(x)));
        assert!(hi - lo > 1.0, "period should vary, got {lo}..{hi}");

        // Events placed exactly at maxima must all land at phase ~0 regardless.
        for &t in m.iter().skip(1).take(20) {
            let p = f.phase_at(t).unwrap();
            assert!(p < 1e-6 || (TAU - p) < 1e-6, "phase {p} at a maximum");
        }
    }

    #[test]
    fn cycle_amplitude_tracks_a_modulated_envelope() {
        // Carrier of period 10 with a slow envelope: amplitude must vary with it.
        let f = sample(
            |t| (1.0 + 0.9 * (TAU * t / 1000.0).sin()) * (TAU * t / 10.0).sin(),
            40_000,
            0.05,
        );
        let m = f.maxima();
        let amps: Vec<f64> = m
            .iter()
            .skip(1)
            .take(m.len() - 2)
            .filter_map(|&t| f.cycle_amplitude_at(t + 0.1))
            .collect();
        let (lo, hi) = amps
            .iter()
            .fold((f64::MAX, f64::MIN), |(l, h), &x| (l.min(x), h.max(x)));
        // Envelope runs 0.1 to 1.9, so peak-to-trough spans roughly 0.2 to 3.8.
        assert!(lo < 0.6, "min amplitude {lo}");
        assert!(hi > 3.0, "max amplitude {hi}");
    }

    #[test]
    fn cycle_amplitude_is_none_outside_the_span() {
        let f = sample(|t| (TAU * t / 10.0).sin(), 10_000, 0.05);
        assert!(f.cycle_amplitude_at(-100.0).is_none());
        assert!(f.cycle_amplitude_at(1e9).is_none());
    }

    #[test]
    fn rejects_mismatched_or_short_input() {
        assert!(Forcing::new(vec![0.0, 1.0], vec![0.0]).is_none());
        assert!(Forcing::new(vec![0.0, 1.0], vec![0.0, 1.0]).is_none());
    }
}

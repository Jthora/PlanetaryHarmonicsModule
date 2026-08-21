//! The generalised Schuster test and time-shifted null distributions.
//!
//! Schuster (1897) tests whether event phases cluster. For `N` events with phases
//! `θᵢ`, the order-`n` resultant is
//!
//! ```text
//! D²ₙ = (Σ cos nθᵢ)² + (Σ sin nθᵢ)²      p = exp(−D²ₙ / N)
//! ```
//!
//! The classical test used throughout the tidal-triggering literature is the
//! `n = 1` case; generalising to arbitrary `n` is how "does base n matter?" becomes
//! a rigorous question (`docs/02-angular-encoding.md`).
//!
//! ⚠ The analytic p-value assumes **independent, uniformly sampled** phases.
//! Earthquake catalogues violate both — aftershock clustering makes them strongly
//! red. Treat `p` as a screening statistic and establish significance with
//! [`time_shift_null`], which preserves the catalogue's own correlation structure.
//!
//! # ⚠ The time-shift null is degenerate at a single exact frequency
//!
//! Found while writing this module's tests, and it constrains how the null may be
//! used.
//!
//! If the forcing phase is a pure sinusoid of fixed period `P`, then shifting every
//! event time by `δ` shifts every phase by the same `2πδ/P`. That **rotates** the
//! phase cluster without changing its concentration, so `D²ₙ` is *exactly
//! invariant* and the null has **zero power**.
//!
//! The null draws its power from the forcing being **quasi-periodic**, not
//! periodic. Real tidal phase is frequency-modulated by the 18.61 yr nodal cycle
//! and by perigee precession, and the full signal is many beating constituents —
//! so a shift genuinely decorrelates forcing from catalogue.
//!
//! **Consequence:** compute phase from the *full* quasi-periodic forcing, never
//! from an idealised single constituent. For a genuinely single-frequency test,
//! the time-shift null must be replaced — by phase randomisation, or by an
//! analytic Schuster p-value with a red-noise correction.

/// Result of a Schuster test at one harmonic order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Schuster {
    pub order: usize,
    /// Squared resultant length D²ₙ.
    pub d_squared: f64,
    /// Analytic p-value exp(−D²ₙ/N). Screening only — see the module note.
    pub p_value: f64,
    /// Mean resultant phase, radians.
    pub phase: f64,
}

/// Schuster statistics for orders `1..=max_order`.
pub fn schuster(phases: &[f64], max_order: usize) -> Vec<Schuster> {
    assert!(max_order > 0, "max_order must be at least 1");
    let n = phases.len();
    (1..=max_order)
        .map(|order| {
            let (mut a, mut b) = (0.0f64, 0.0f64);
            for &p in phases {
                let (s, c) = (order as f64 * p).sin_cos();
                a += c;
                b += s;
            }
            let d_squared = a * a + b * b;
            let p_value = if n == 0 {
                1.0
            } else {
                (-d_squared / n as f64).exp()
            };
            Schuster {
                order,
                d_squared,
                p_value,
                phase: b.atan2(a),
            }
        })
        .collect()
}

/// A null distribution built by shifting event times against the forcing.
#[derive(Debug, Clone)]
pub struct NullDistribution {
    /// Statistic from each shifted realisation, ascending.
    pub samples: Vec<f64>,
}

impl NullDistribution {
    /// Fraction of null samples at or above `observed` — an empirical p-value.
    ///
    /// ⚠ **The floor is `1/(n+1)`** for `n` samples, from the add-one correction.
    /// Reaching p < 0.05 needs at least 20 shifts; p < 0.005 needs at least 200.
    /// Budget the shift count for the significance you intend to claim.
    pub fn p_value(&self, observed: f64) -> f64 {
        if self.samples.is_empty() {
            return 1.0;
        }
        let ge = self.samples.iter().filter(|&&s| s >= observed).count();
        // Add-one correction: never report exactly zero from finite sampling.
        (ge as f64 + 1.0) / (self.samples.len() as f64 + 1.0)
    }
}

/// Build a null distribution by offsetting event times.
///
/// Shifting the whole catalogue preserves its internal clustering *and* the
/// forcing's structure, breaking only their alignment — unlike random shuffling,
/// which destroys autocorrelation and yields falsely tight nulls
/// (`docs/04-ml-architecture.md` §6a).
///
/// `phase_at` maps an event time to its forcing phase. Avoid offsets near integer
/// multiples of dominant periods (1 yr, 18.61 yr), which partially preserve
/// alignment; [`shift_offsets`] does this for you.
pub fn time_shift_null<F>(
    times: &[f64],
    offsets: &[f64],
    order: usize,
    mut phase_at: F,
) -> NullDistribution
where
    F: FnMut(f64) -> f64,
{
    let mut samples: Vec<f64> = offsets
        .iter()
        .map(|&off| {
            let phases: Vec<f64> = times.iter().map(|&t| phase_at(t + off)).collect();
            schuster(&phases, order)[order - 1].d_squared
        })
        .collect();
    samples.sort_by(|a, b| a.partial_cmp(b).unwrap());
    NullDistribution { samples }
}

/// Evenly spaced ± offsets that avoid resonance with the given periods.
///
/// `min`, `max`, and `avoid` share units with the event times (seconds by
/// convention here). An offset is rejected when it falls within `tol` of an
/// integer multiple of any avoided period.
pub fn shift_offsets(count: usize, min: f64, max: f64, avoid: &[f64], tol: f64) -> Vec<f64> {
    let mut out = Vec::with_capacity(count * 2);
    if count == 0 || max <= min {
        return out;
    }
    let step = (max - min) / count as f64;
    for i in 0..count {
        let mag = min + step * i as f64;
        let resonant = avoid.iter().any(|&p| {
            if p <= 0.0 {
                return false;
            }
            let r = (mag / p).round() * p;
            (mag - r).abs() < tol
        });
        if !resonant {
            out.push(mag);
            out.push(-mag);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    #[test]
    fn uniform_phases_are_not_significant() {
        let n = 5000;
        let phases: Vec<f64> = (0..n).map(|i| TAU * i as f64 / n as f64).collect();
        for s in schuster(&phases, 4) {
            assert!(s.p_value > 0.5, "order {} p={}", s.order, s.p_value);
        }
    }

    #[test]
    fn clustered_phases_are_significant() {
        let phases = vec![1.0f64; 400];
        let s = schuster(&phases, 1);
        assert!(s[0].p_value < 1e-12, "p={}", s[0].p_value);
    }

    #[test]
    fn null_p_value_is_bounded_and_never_zero() {
        let d = NullDistribution {
            samples: vec![1.0, 2.0, 3.0, 4.0],
        };
        assert!(d.p_value(100.0) > 0.0);
        assert!((d.p_value(0.0) - 1.0).abs() < 1e-12);
        let mid = d.p_value(3.0);
        assert!(mid > 0.0 && mid < 1.0, "mid={mid}");
    }

    #[test]
    fn shift_offsets_skip_resonant_lags() {
        // Avoid multiples of 10 with tolerance 1: 10, 20, 30 must not appear.
        let offs = shift_offsets(40, 1.0, 41.0, &[10.0], 1.0);
        for o in offs {
            let a = o.abs();
            assert!((a - 10.0).abs() >= 1.0 && (a - 20.0).abs() >= 1.0);
        }
    }

    /// Quasi-periodic forcing: carrier of period `P`, frequency-modulated by a
    /// much slower cycle — a stand-in for nodal and perigee modulation of the
    /// real tidal phase.
    fn fm_phase(t: f64) -> f64 {
        let p = 100.0;
        let p_mod = 9_000.0;
        super::super::harmonics::wrap_2pi(TAU * (t / p + 0.45 * (TAU * t / p_mod).sin()))
    }

    #[test]
    fn pure_sinusoid_null_is_degenerate() {
        // Documented caveat: a global shift rotates the phase cluster without
        // changing its concentration, so D² is invariant and the null is powerless.
        let period = 100.0;
        let times: Vec<f64> = (0..200).map(|i| i as f64 * period + 12.0).collect();
        let phase_at = |t: f64| TAU * (t / period).fract();

        let observed = schuster(
            &times.iter().map(|&t| phase_at(t)).collect::<Vec<_>>(),
            1,
        )[0]
        .d_squared;
        let null = time_shift_null(&times, &shift_offsets(40, 7.0, 400.0, &[], 0.0), 1, phase_at);

        for s in &null.samples {
            assert!((s - observed).abs() / observed < 1e-9, "D² should be invariant");
        }
    }

    #[test]
    fn time_shift_null_flags_a_real_signal() {
        // Events selected where the quasi-periodic forcing is near phase zero.
        let times: Vec<f64> = (0..200_000)
            .map(|i| i as f64 * 0.5)
            .filter(|&t| {
                let p = fm_phase(t);
                p < 0.10 || p > TAU - 0.10
            })
            .collect();
        assert!(times.len() > 200, "got {} events", times.len());

        let observed = schuster(
            &times.iter().map(|&t| fm_phase(t)).collect::<Vec<_>>(),
            1,
        )[0]
        .d_squared;

        // Avoid the modulation period only. Filtering on the 100 s carrier too
        // would reject nearly every offset and drive the p-value floor above 0.05.
        let offsets = shift_offsets(80, 500.0, 40_000.0, &[9_000.0], 300.0);
        assert!(offsets.len() >= 20, "need >=20 shifts for p<0.05, got {}", offsets.len());

        let null = time_shift_null(&times, &offsets, 1, fm_phase);
        let p = null.p_value(observed);
        assert!(p < 0.05, "p={p}, null max={:?}", null.samples.last());
    }
}

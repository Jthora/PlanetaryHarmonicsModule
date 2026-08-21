//! Fourier angular encodings and time-domain harmonic decomposition.
//!
//! Angles are encoded as `[cos nθ, sin nθ]` pairs rather than binned into sectors.
//! Keeping the sine term is what lets a model learn a **phase offset**, and phase
//! offset is where the tidal signal lives — see `docs/02-angular-encoding.md`.
//!
//! The cos/sin ladder is built by angle-sum recurrence: O(N) multiply-adds rather
//! than N transcendental calls.

/// Fourier features `[cos θ, sin θ, cos 2θ, sin 2θ, …, cos Nθ, sin Nθ]`.
///
/// Output length is `2 * max_order`. Panics if `max_order` is zero.
pub fn fourier_features(theta: f64, max_order: usize) -> Vec<f64> {
    assert!(max_order > 0, "max_order must be at least 1");
    let (s1, c1) = theta.sin_cos();
    let mut out = Vec::with_capacity(2 * max_order);
    let (mut c, mut s) = (c1, s1);
    for n in 1..=max_order {
        out.push(c);
        out.push(s);
        if n < max_order {
            // (c, s) at order n+1 by angle-sum from order n.
            let cn = c * c1 - s * s1;
            let sn = s * c1 + c * s1;
            c = cn;
            s = sn;
        }
    }
    out
}

/// Amplitude and phase of a harmonic fitted to phase-angle samples.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Harmonic {
    /// Harmonic order n.
    pub order: usize,
    /// Amplitude √(a² + b²), normalised by sample count.
    pub amplitude: f64,
    /// Phase atan2(b, a), radians in (−π, π].
    pub phase: f64,
}

/// Fit amplitude and phase at each order from a set of phase angles.
///
/// This is the resultant-vector calculation underlying the Schuster test; see
/// [`crate::stats::schuster`] for the significance of each order.
pub fn decompose(phases: &[f64], max_order: usize) -> Vec<Harmonic> {
    assert!(max_order > 0, "max_order must be at least 1");
    let n = phases.len();
    (1..=max_order)
        .map(|order| {
            let (mut a, mut b) = (0.0f64, 0.0f64);
            for &p in phases {
                let x = order as f64 * p;
                let (s, c) = x.sin_cos();
                a += c;
                b += s;
            }
            let amplitude = if n == 0 {
                0.0
            } else {
                (a * a + b * b).sqrt() / n as f64
            };
            Harmonic {
                order,
                amplitude,
                phase: b.atan2(a),
            }
        })
        .collect()
}

/// Wrap an angle into `[0, 2π)`.
pub fn wrap_2pi(theta: f64) -> f64 {
    let two_pi = std::f64::consts::TAU;
    let r = theta % two_pi;
    if r < 0.0 {
        r + two_pi
    } else {
        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{PI, TAU};

    #[test]
    fn recurrence_matches_direct_evaluation() {
        let theta = 0.7213;
        let f = fourier_features(theta, 32);
        for n in 1..=32usize {
            let (s, c) = (n as f64 * theta).sin_cos();
            assert!((f[2 * (n - 1)] - c).abs() < 1e-12, "cos order {n}");
            assert!((f[2 * (n - 1) + 1] - s).abs() < 1e-12, "sin order {n}");
        }
    }

    #[test]
    fn uniform_phases_have_near_zero_amplitude() {
        let n = 2000;
        let phases: Vec<f64> = (0..n).map(|i| TAU * i as f64 / n as f64).collect();
        for h in decompose(&phases, 8) {
            assert!(h.amplitude < 1e-9, "order {} leaked {}", h.order, h.amplitude);
        }
    }

    #[test]
    fn concentrated_phases_recover_order_and_phase() {
        // All events at θ = π/3 → every order has amplitude 1, phase n·π/3.
        let phases = vec![PI / 3.0; 500];
        let hs = decompose(&phases, 3);
        for h in &hs {
            assert!((h.amplitude - 1.0).abs() < 1e-12);
            let expected = wrap_2pi(h.order as f64 * PI / 3.0);
            assert!((wrap_2pi(h.phase) - expected).abs() < 1e-9);
        }
    }

    #[test]
    fn wrap_handles_negatives() {
        assert!((wrap_2pi(-0.5) - (TAU - 0.5)).abs() < 1e-12);
        assert!((wrap_2pi(TAU + 0.25) - 0.25).abs() < 1e-12);
    }
}

//! Angle-domain event finding — solving for crossings instead of sampling for them.
//!
//! # The problem this replaces
//!
//! "When does this angle reach this value?" is usually answered by sampling the
//! angle densely and looking for sign changes. That costs `O(sample rate)` and its
//! precision is bounded by the step: resolving an event to one second over a
//! century needs 3.2 billion evaluations.
//!
//! But the angles that matter here — Doodson arguments, relative longitudes,
//! aspect angles — are **near-linear in time**. A body's mean motion dominates;
//! periodic corrections are small. So crossings can be *predicted* analytically
//! and refined by Newton iteration:
//!
//! ```text
//! 1. linear estimate from the mean rate       ->  t_k
//! 2. Newton refine on theta(t) - target       ->  t*
//! 3. converges in 2-3 iterations
//! ```
//!
//! Cost becomes `O(number of events)`, independent of the precision requested, and
//! the answer is exact to machine tolerance rather than to a sample interval.
//!
//! This is what makes millisecond-resolution event timing tractable over long
//! spans — the requirement `docs/06-engine-architecture.md` §2 identifies for
//! Star Seer.
//!
//! # Scope
//!
//! Works for any angle with a well-defined mean rate and no reversals over the
//! search interval. **Retrograde motion breaks the assumption**: near a station a
//! relative longitude can cross the same target three times, and the linear
//! predictor will find one. [`crossings_bracketed`] handles that case by scanning
//! coarsely for sign changes first, trading some speed for correctness.

use std::f64::consts::TAU;

/// One located crossing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Crossing {
    /// Time of the crossing, in the caller's units.
    pub time: f64,
    /// Newton iterations used. High counts indicate a poor rate estimate.
    pub iterations: u32,
    /// Achieved residual `|θ(t) − target|` in radians.
    ///
    /// Compare against the tolerance you asked for. A residual above it means the
    /// precision floor was hit, not that the crossing is spurious.
    pub residual: f64,
}

/// Wrap to `(−π, π]`.
fn wrap_pi(x: f64) -> f64 {
    let mut r = x % TAU;
    if r > std::f64::consts::PI {
        r -= TAU;
    }
    if r <= -std::f64::consts::PI {
        r += TAU;
    }
    r
}

/// Find all times in `[t0, t1]` where `theta` reaches `target` (mod 2π).
///
/// `rate` is the angle's mean rate in radians per unit time — the mean motion.
/// It seeds the linear prediction and serves as the Newton derivative, which is
/// why a near-linear angle converges in a couple of iterations.
///
/// `tol` is the residual in radians at which refinement stops.
///
/// ⚠ **There is a precision floor.** `θ(t)` is evaluated unwrapped, so its absolute
/// magnitude grows with the span, and the achievable residual is bounded below by
/// roughly `|θ| · ε`. At an angle of magnitude 700 that is ~1.5e-13 radians, and a
/// smaller `tol` cannot be met. Crossings are **never silently dropped** for
/// failing to reach `tol` — refinement stops when it stops improving, and the
/// achieved [`Crossing::residual`] is reported for the caller to judge.
///
/// Returns crossings in ascending time. Cost is proportional to the **number of
/// crossings**, not to the span divided by the tolerance.
pub fn crossings<F>(theta: F, rate: f64, target: f64, t0: f64, t1: f64, tol: f64) -> Vec<Crossing>
where
    F: Fn(f64) -> f64,
{
    let mut out = Vec::new();
    if !(rate.is_finite()) || rate == 0.0 || t1 <= t0 || !tol.is_finite() || tol <= 0.0 {
        return out;
    }
    let period = TAU / rate.abs();

    // How far past the target is theta at t0? Step forward one cycle at a time.
    let start_offset = wrap_pi(target - theta(t0));
    let first = t0 + start_offset / rate;
    // Step back one period if the linear estimate landed before the window.
    let mut k = if first < t0 { 1.0 } else { 0.0 };

    loop {
        let guess = first + k * period * rate.signum() * rate.signum();
        if guess > t1 + period {
            break;
        }
        k += 1.0;
        if guess < t0 - period {
            continue;
        }

        // Newton on wrap_pi(theta(t) - target), derivative approximated by `rate`.
        // Stop at `tol`, or earlier if iteration stops improving -- past the
        // precision floor the residual plateaus and further steps only add noise.
        let mut t = guess;
        let mut iters = 0u32;
        let mut resid = wrap_pi(theta(t) - target);
        while resid.abs() > tol && iters < 32 {
            let next_t = t - resid / rate;
            let next_resid = wrap_pi(theta(next_t) - target);
            if next_resid.abs() >= resid.abs() {
                break;
            }
            t = next_t;
            resid = next_resid;
            iters += 1;
        }
        if t >= t0 && t <= t1 {
            // Guard against the predictor and Newton landing on the same root twice.
            if out
                .last()
                .is_none_or(|c: &Crossing| (t - c.time).abs() > period * 0.5)
            {
                out.push(Crossing {
                    time: t,
                    iterations: iters,
                    residual: resid.abs(),
                });
            }
        }
    }
    out.sort_by(|a, b| a.time.partial_cmp(&b.time).unwrap());
    out
}

/// Crossings located by coarse bracketing then bisection-Newton.
///
/// Slower than [`crossings`] but tolerant of **reversals** — retrograde motion,
/// where an angle can cross the same target several times within one nominal
/// period and a linear predictor finds only one.
///
/// `coarse_step` must be small enough that no pair of crossings falls inside a
/// single step; a few samples per nominal period is the usual choice.
pub fn crossings_bracketed<F>(
    theta: F,
    target: f64,
    t0: f64,
    t1: f64,
    coarse_step: f64,
    tol: f64,
) -> Vec<Crossing>
where
    F: Fn(f64) -> f64,
{
    let mut out = Vec::new();
    if t1 <= t0 || coarse_step <= 0.0 || tol <= 0.0 {
        return out;
    }
    let mut a = t0;
    let mut fa = wrap_pi(theta(a) - target);
    while a < t1 {
        let b = (a + coarse_step).min(t1);
        let fb = wrap_pi(theta(b) - target);
        // A sign change with both values small brackets a root; a sign change with
        // large values is the +/-pi wrap, not a crossing.
        if fa.signum() != fb.signum() && fa.abs() < 1.0 && fb.abs() < 1.0 {
            let (mut lo, mut hi, mut flo) = (a, b, fa);
            let mut iters = 0u32;
            while hi - lo > tol.abs() && iters < 64 {
                let mid = 0.5 * (lo + hi);
                let fm = wrap_pi(theta(mid) - target);
                if fm.signum() == flo.signum() {
                    lo = mid;
                    flo = fm;
                } else {
                    hi = mid;
                }
                iters += 1;
            }
            let t = 0.5 * (lo + hi);
            out.push(Crossing {
                time: t,
                iterations: iters,
                residual: wrap_pi(theta(t) - target).abs(),
            });
        }
        a = b;
        fa = fb;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    #[test]
    fn locates_crossings_of_a_pure_linear_angle() {
        // theta = 0.7 t, so theta = 0 every 2*pi/0.7 = 8.976 units.
        let rate = 0.7;
        let c = crossings(|t| rate * t, rate, 0.0, 0.0, 100.0, 1e-12);
        let period = TAU / rate;
        assert!(c.len() >= 11, "got {} crossings", c.len());
        for (i, x) in c.iter().enumerate() {
            assert!(
                (x.time - i as f64 * period).abs() < 1e-9,
                "crossing {i} at {} expected {}",
                x.time,
                i as f64 * period
            );
            assert!(x.residual < 1e-12);
        }
    }

    #[test]
    fn converges_in_a_couple_of_iterations_for_near_linear_angles() {
        // Linear plus a small periodic wobble, as a real argument behaves.
        let rate = 0.5;
        let f = |t: f64| rate * t + 0.02 * (0.31 * t).sin();
        let c = crossings(f, rate, 1.0, 0.0, 200.0, 1e-10);
        assert!(!c.is_empty());
        // Newton uses the mean rate as the derivative, which is ~1% off here, so
        // convergence is linear rather than quadratic -- a handful of steps.
        let worst = c.iter().map(|x| x.iterations).max().unwrap();
        assert!(worst <= 10, "worst case took {worst} iterations");
        assert!(c.iter().all(|x| x.residual < 1e-9));
    }

    #[test]
    fn costs_the_number_of_events_not_the_precision() {
        // The point of the module: evaluations should not scale with tolerance.
        let rate = 0.7;
        let count = |tol: f64| -> u32 {
            let n = Cell::new(0u32);
            let f = |t: f64| {
                n.set(n.get() + 1);
                rate * t
            };
            let c = crossings(f, rate, 0.0, 0.0, 1000.0, tol);
            assert!(c.len() > 100);
            n.get()
        };
        let coarse = count(1e-6);
        let fine = count(1e-11);
        // A sampling approach would need 10^5 times more work for 10^5 more
        // precision; this needs at most a few extra iterations per event.
        assert!(
            fine < coarse * 2,
            "evaluations grew from {coarse} to {fine} for 5 more digits"
        );
    }

    #[test]
    fn respects_the_search_window() {
        let rate = 1.0;
        let c = crossings(|t| rate * t, rate, 0.0, 10.0, 30.0, 1e-12);
        assert!(c.iter().all(|x| x.time >= 10.0 && x.time <= 30.0));
        assert!(!c.is_empty());
    }

    #[test]
    fn rejects_degenerate_input() {
        assert!(crossings(|t| t, 0.0, 0.0, 0.0, 10.0, 1e-9).is_empty());
        assert!(crossings(|t| t, 1.0, 0.0, 10.0, 0.0, 1e-9).is_empty());
        assert!(crossings(|t| t, 1.0, 0.0, 0.0, 10.0, 0.0).is_empty());
    }

    #[test]
    fn bracketed_finder_matches_a_dense_scan_through_a_reversal() {
        // Retrograde-like: the angle advances, reverses, and advances again, so it
        // crosses the target more than once per nominal period. Ground truth comes
        // from a dense scan rather than hand-computed roots -- the property that
        // matters is agreement with brute force.
        let f = |t: f64| 0.2 * t + 1.5 * (t).sin();
        let target = 3.0;
        let (t0, t1) = (0.0, 25.0);

        let mut truth = Vec::new();
        let step = 1e-4;
        let mut a = t0;
        let mut fa = wrap_pi(f(a) - target);
        while a < t1 {
            let b = a + step;
            let fb = wrap_pi(f(b) - target);
            if fa.signum() != fb.signum() && fa.abs() < 1.0 && fb.abs() < 1.0 {
                truth.push(0.5 * (a + b));
            }
            a = b;
            fa = fb;
        }
        assert!(truth.len() >= 3, "test function should reverse: {truth:?}");

        let found = crossings_bracketed(f, target, t0, t1, 0.05, 1e-9);
        assert_eq!(
            found.len(),
            truth.len(),
            "found {:?} vs dense scan {truth:?}",
            found.iter().map(|c| c.time).collect::<Vec<_>>()
        );
        for (c, &t) in found.iter().zip(&truth) {
            assert!((c.time - t).abs() < 1e-3, "{} vs {t}", c.time);
            assert!(c.residual < 1e-6);
        }
    }

    #[test]
    fn both_finders_agree_on_a_monotonic_angle() {
        let rate = 0.4;
        let f = |t: f64| rate * t + 0.01 * (0.7 * t).cos();
        let fast = crossings(f, rate, 0.0, 0.0, 120.0, 1e-10);
        let safe = crossings_bracketed(f, 0.0, 0.0, 120.0, 1.0, 1e-10);
        assert_eq!(fast.len(), safe.len(), "{fast:?} vs {safe:?}");
        for (a, b) in fast.iter().zip(&safe) {
            assert!((a.time - b.time).abs() < 1e-6);
        }
    }
}

//! Doodson fundamental arguments and constituent phases.
//!
//! C4 failed because constituent phase was derived by *demodulating* the composite
//! ΔCFS series. That made each band a near-pure tone (killing the time-shift null)
//! and made tide-free frequencies return the leaking constituent's phase rather
//! than a neutral baseline. See `docs/07`, traps 5 and 6.
//!
//! The fix is to compute constituent phase **analytically**, which is what the
//! tidal literature does. Each constituent's phase is an integer combination of
//! six slowly-varying astronomical arguments:
//!
//! | Symbol | Argument | Period |
//! |---|---|---|
//! | `τ`  | lunar time | 24.84 h |
//! | `s`  | Moon's mean longitude | 27.32 d |
//! | `h`  | Sun's mean longitude | 365.24 d |
//! | `p`  | longitude of lunar perigee | 8.85 yr |
//! | `N′` | negative longitude of lunar ascending node | 18.61 yr |
//! | `pₛ` | longitude of solar perigee | 20 940 yr |
//!
//! ```text
//! argument = n₁τ + n₂s + n₃h + n₄p + n₅N′ + n₆pₛ
//! ```
//!
//! **Why this fixes the null.** Each argument advances essentially linearly, so
//! the combined phase is uniform over long spans *by construction*. Events uniform
//! in time are then uniform in phase, which is exactly the condition the Schuster
//! statistic assumes and the condition demodulated phase failed to meet.
//!
//! Only the *arguments* are needed for phase; constituent **amplitudes** still
//! require a catalogue (HW95/KSM03). Phase is what the statistics need, so this
//! module unblocks C4 on its own.
//!
//! # Time and precision
//!
//! Input is **days since 2000-01-01T00:00 UTC**, matching the catalogue modules.
//! Two approximations, both far below what matters here:
//!
//! - UTC is used where UT1 is strictly wanted. Leap seconds bound `|UT1−UTC|` to
//!   0.9 s, i.e. 0.004° of Earth rotation.
//! - UTC is used where TT is strictly wanted in the polynomials. Over 2000–2025
//!   that is ~69 s, about 0.55° of M2 phase.
//!
//! Both are systematic but tiny. Applying them properly needs a leap-second
//! kernel; see `docs/10` on time-scale correctness mattering more than ephemeris
//! precision.



/// The six Doodson fundamental arguments at a time, in degrees.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Arguments {
    /// Lunar time τ.
    pub tau: f64,
    /// Moon's mean longitude s.
    pub s: f64,
    /// Sun's mean longitude h.
    pub h: f64,
    /// Longitude of lunar perigee p.
    pub p: f64,
    /// **Negative** longitude of the lunar ascending node, N′.
    pub n_prime: f64,
    /// Longitude of solar perigee pₛ.
    pub p_s: f64,
}

/// Evaluate the fundamental arguments at `days` since 2000-01-01T00:00 UTC.
pub fn arguments(days: f64) -> Arguments {
    // Julian centuries from J2000.
    let t = days / 36525.0;
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;

    let s = 218.3164477 + 481_267.881_234_21 * t - 0.0015786 * t2 + t3 / 538_841.0
        - t4 / 65_194_000.0;
    let h = 280.46646 + 36_000.76983 * t + 0.0003032 * t2;
    let p = 83.3532465 + 4069.013_728_7 * t - 0.0103200 * t2 - t3 / 80_053.0
        + t4 / 18_999_000.0;
    let node = 125.0445479 - 1934.136_289_1 * t + 0.0020754 * t2 + t3 / 467_441.0
        - t4 / 60_616_000.0;
    let p_s = 282.93735 + 1.71946 * t + 0.00046 * t2;

    // Greenwich mean sidereal time, degrees. `days` is measured from J2000 at
    // 00:00 UTC while GMST is referenced to J2000 at 12:00 TT, hence the -0.5.
    let gmst = 280.460_618_37 + 360.985_647_366_29 * (days - 0.5) + 0.000_387_933 * t2
        - t3 / 38_710_000.0;

    // Lunar time: mean solar-style angle of the Moon, tau = GMST + 180 - s.
    let tau = gmst + 180.0 - s;

    Arguments {
        tau: norm360(tau),
        s: norm360(s),
        h: norm360(h),
        p: norm360(p),
        n_prime: norm360(-node),
        p_s: norm360(p_s),
    }
}

fn norm360(x: f64) -> f64 {
    let r = x % 360.0;
    if r < 0.0 {
        r + 360.0
    } else {
        r
    }
}

/// A tidal constituent, identified by its Doodson coefficients.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Constituent {
    pub name: &'static str,
    /// Coefficients `(n₁ … n₆)` of `(τ, s, h, p, N′, pₛ)`.
    pub n: [i32; 6],
}

impl Constituent {
    /// Phase of this constituent at `days`, in `[0, 2π)`.
    ///
    /// Zero phase is a maximum of the constituent's contribution to the tide-
    /// generating potential.
    pub fn phase_at(&self, days: f64) -> f64 {
        let a = arguments(days);
        let deg = self.n[0] as f64 * a.tau
            + self.n[1] as f64 * a.s
            + self.n[2] as f64 * a.h
            + self.n[3] as f64 * a.p
            + self.n[4] as f64 * a.n_prime
            + self.n[5] as f64 * a.p_s;
        norm360(deg).to_radians()
    }

    /// Phases for a set of times.
    pub fn phases(&self, days: &[f64]) -> Vec<f64> {
        days.iter().map(|&d| self.phase_at(d)).collect()
    }

    /// Phase at a given **east longitude**, in `[0, 2π)`.
    ///
    /// [`Self::phase_at`] uses Greenwich lunar time. Tidal phase is *local*: the
    /// lunar time argument τ advances with longitude, so
    ///
    /// ```text
    /// τ_local = τ_Greenwich + longitude
    /// ```
    ///
    /// and the constituent argument gains `n₁ · λ`.
    ///
    /// **Essential for a global catalogue.** For semidiurnal constituents `n₁ = 2`,
    /// so a 180° longitude error is a full 360° of phase — events would be
    /// scattered across every phase and any real signal erased. Long-period
    /// constituents have `n₁ = 0` and are unaffected.
    pub fn phase_at_longitude(&self, days: f64, lon_deg: f64) -> f64 {
        let base = self.phase_at(days).to_degrees();
        norm360(base + self.n[0] as f64 * lon_deg).to_radians()
    }

    /// Phases for times paired with east longitudes.
    pub fn phases_at_longitudes(&self, days: &[f64], lons: &[f64]) -> Vec<f64> {
        days.iter()
            .zip(lons)
            .map(|(&d, &l)| self.phase_at_longitude(d, l))
            .collect()
    }

    /// Period in days, from the argument's mean rate of advance.
    ///
    /// Measured by central difference over a century rather than assumed, so it
    /// doubles as a check that the coefficients are right.
    pub fn period_days(&self) -> f64 {
        let (a, b) = (arguments(-18262.5), arguments(18262.5));
        let rate = |x: &Arguments| {
            self.n[0] as f64 * x.tau
                + self.n[1] as f64 * x.s
                + self.n[2] as f64 * x.h
                + self.n[3] as f64 * x.p
                + self.n[4] as f64 * x.n_prime
                + self.n[5] as f64 * x.p_s
        };
        // Use exact polynomial rates rather than wrapped values.
        let span = 36525.0;
        let cycles = (rate_unwrapped(self, 18262.5) - rate_unwrapped(self, -18262.5)) / 360.0;
        let _ = (a, b, rate);
        span / cycles.abs()
    }
}

/// Argument in degrees without wrapping, for rate calculations.
fn rate_unwrapped(c: &Constituent, days: f64) -> f64 {
    let t = days / 36525.0;
    let t2 = t * t;
    let s = 218.3164477 + 481_267.881_234_21 * t - 0.0015786 * t2;
    let h = 280.46646 + 36_000.76983 * t + 0.0003032 * t2;
    let p = 83.3532465 + 4069.013_728_7 * t - 0.0103200 * t2;
    let node = 125.0445479 - 1934.136_289_1 * t + 0.0020754 * t2;
    let p_s = 282.93735 + 1.71946 * t + 0.00046 * t2;
    let gmst = 280.460_618_37 + 360.985_647_366_29 * (days - 0.5) + 0.000_387_933 * t2;
    let tau = gmst + 180.0 - s;
    c.n[0] as f64 * tau
        + c.n[1] as f64 * s
        + c.n[2] as f64 * h
        + c.n[3] as f64 * p
        + c.n[4] as f64 * (-node)
        + c.n[5] as f64 * p_s
}

/// Named constituents, with Doodson coefficients `(τ, s, h, p, N′, pₛ)`.
pub const CONSTITUENTS: &[Constituent] = &[
    // Semidiurnal
    Constituent { name: "M2", n: [2, 0, 0, 0, 0, 0] },
    Constituent { name: "S2", n: [2, 2, -2, 0, 0, 0] },
    Constituent { name: "N2", n: [2, -1, 0, 1, 0, 0] },
    Constituent { name: "K2", n: [2, 2, 0, 0, 0, 0] },
    // Diurnal
    Constituent { name: "K1", n: [1, 1, 0, 0, 0, 0] },
    Constituent { name: "O1", n: [1, -1, 0, 0, 0, 0] },
    Constituent { name: "P1", n: [1, 1, -2, 0, 0, 0] },
    Constituent { name: "Q1", n: [1, -2, 0, 1, 0, 0] },
    // Long period
    Constituent { name: "Mf", n: [0, 2, 0, 0, 0, 0] },
    Constituent { name: "Msf", n: [0, 2, -2, 0, 0, 0] },
    Constituent { name: "Mm", n: [0, 1, 0, -1, 0, 0] },
    Constituent { name: "Ssa", n: [0, 0, 2, 0, 0, 0] },
    Constituent { name: "Sa", n: [0, 0, 1, 0, 0, -1] },
];

/// Look up a constituent by name.
pub fn constituent(name: &str) -> Option<&'static Constituent> {
    CONSTITUENTS.iter().find(|c| c.name == name)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Published constituent periods, in days.
    const KNOWN: &[(&str, f64)] = &[
        ("M2", 0.5175),
        ("S2", 0.5000),
        ("N2", 0.5274),
        ("K2", 0.4986),
        ("K1", 0.9973),
        ("O1", 1.0758),
        ("P1", 1.0027),
        ("Q1", 1.1195),
        ("Mf", 13.6608),
        ("Msf", 14.7653),
        ("Mm", 27.5546),
        ("Ssa", 182.621),
        ("Sa", 365.256),
    ];

    #[test]
    fn every_constituent_period_matches_the_published_value() {
        for (name, expect) in KNOWN {
            let c = constituent(name).unwrap();
            let got = c.period_days();
            let err = (got - expect).abs() / expect;
            assert!(err < 1e-3, "{name}: got {got:.6} d, expected {expect} d");
        }
    }

    #[test]
    fn arguments_stay_in_range() {
        for d in [-20000.0, -1.0, 0.0, 1234.5, 9000.0] {
            let a = arguments(d);
            for x in [a.tau, a.s, a.h, a.p, a.n_prime, a.p_s] {
                assert!((0.0..360.0).contains(&x), "{x} out of range at day {d}");
            }
        }
    }

    #[test]
    fn phases_are_in_range_and_advance_uniformly() {
        let m2 = constituent("M2").unwrap();
        // Sample far more finely than the period, then check the phase histogram
        // is flat -- this is the property demodulated phase failed to have.
        let n = 200_000;
        let mut hist = [0usize; 12];
        for i in 0..n {
            let d = i as f64 * 0.037;
            let p = m2.phase_at(d);
            assert!((0.0..std::f64::consts::TAU).contains(&p));
            hist[(p / std::f64::consts::TAU * 12.0) as usize % 12] += 1;
        }
        let expect = n as f64 / 12.0;
        for (i, &c) in hist.iter().enumerate() {
            let dev = (c as f64 - expect).abs() / expect;
            assert!(dev < 0.05, "bin {i} deviates {:.1}% from uniform", dev * 100.0);
        }
    }

    #[test]
    fn long_period_constituents_are_also_uniform() {
        for name in ["Mf", "Mm", "Sa"] {
            let c = constituent(name).unwrap();
            let n = 100_000;
            let mut hist = [0usize; 12];
            for i in 0..n {
                // ~25 years, the Parkfield span.
                let d = i as f64 * 0.09;
                hist[(c.phase_at(d) / std::f64::consts::TAU * 12.0) as usize % 12] += 1;
            }
            let expect = n as f64 / 12.0;
            let worst = hist
                .iter()
                .map(|&c| (c as f64 - expect).abs() / expect)
                .fold(0.0, f64::max);
            assert!(worst < 0.10, "{name} worst bin deviates {:.1}%", worst * 100.0);
        }
    }

    #[test]
    fn m2_and_s2_are_distinguishable() {
        let (m2, s2) = (constituent("M2").unwrap(), constituent("S2").unwrap());
        // S2 is locked to solar time: its phase must repeat every half solar day.
        let a = s2.phase_at(1000.0);
        let b = s2.phase_at(1000.5);
        assert!((a - b).abs() < 1e-6 || (std::f64::consts::TAU - (a - b).abs()) < 1e-6, "{a} vs {b}");
        // M2 must NOT, since a lunar day is 24.84 h.
        let c = m2.phase_at(1000.0);
        let d = m2.phase_at(1000.5);
        assert!((c - d).abs() > 0.1, "M2 should drift against solar time");
    }

    #[test]
    fn longitude_shifts_semidiurnal_phase_by_twice_the_angle() {
        let m2 = constituent("M2").unwrap();
        let base = m2.phase_at_longitude(1000.0, 0.0);
        // n1 = 2 for M2, so +90 deg longitude is +180 deg of phase.
        let shifted = m2.phase_at_longitude(1000.0, 90.0);
        let d = (shifted - base).rem_euclid(std::f64::consts::TAU);
        assert!(
            (d - std::f64::consts::PI).abs() < 1e-9,
            "expected pi, got {d}"
        );
        // A full 180 deg of longitude is a whole cycle: back to the start.
        let full = m2.phase_at_longitude(1000.0, 180.0);
        assert!((full - base).abs() < 1e-9, "{full} vs {base}");
    }

    #[test]
    fn longitude_does_not_affect_long_period_constituents() {
        for name in ["Mf", "Mm", "Ssa", "Sa"] {
            let c = constituent(name).unwrap();
            let a = c.phase_at_longitude(1000.0, 0.0);
            let b = c.phase_at_longitude(1000.0, 137.0);
            assert!((a - b).abs() < 1e-12, "{name} moved with longitude");
        }
    }

    #[test]
    fn phase_at_longitude_zero_matches_phase_at() {
        let m2 = constituent("M2").unwrap();
        for d in [0.0, 123.4, 9000.0] {
            assert!((m2.phase_at_longitude(d, 0.0) - m2.phase_at(d)).abs() < 1e-12);
        }
    }

    #[test]
    fn unknown_constituent_is_none() {
        assert!(constituent("XYZ").is_none());
    }
}

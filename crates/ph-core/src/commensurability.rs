//! Multi-body angular commensurabilities — the harmonic basis beyond base 12.
//!
//! The project's founding idea was to stop privileging the twelvefold division and
//! test *every* base. [`crate::harmonics`] does that for a single angle. This
//! module does it for **relationships among several bodies**.
//!
//! For angles `θ₁ … θₙ` the general combination is
//!
//! ```text
//! Φ_k(t) = Σ kᵢ θᵢ(t)          kᵢ ∈ ℤ
//! ```
//!
//! # The d'Alembert constraint is the whole point
//!
//! Combinations are restricted to those with
//!
//! ```text
//! Σ kᵢ = 0
//! ```
//!
//! This is the d'Alembert rule from the classical disturbing-function expansion,
//! and it enforces **rotational invariance**: no physical quantity can depend on
//! where the origin of longitude was arbitrarily placed.
//!
//! Two consequences, and they are the same consequence:
//!
//! 1. It is required for physical correctness.
//! 2. It **automatically removes every feature that depends on absolute zodiacal
//!    position** while keeping every feature built on relative geometry.
//!
//! The constraint that makes the physics right is the constraint that discards the
//! indefensible features. It also collapses the combinatorial space substantially,
//! which matters for the multiple-testing budget.
//!
//! # Frequencies, and why they decide usefulness
//!
//! Given mean motions `ṅᵢ`, a combination advances at `Σ kᵢ ṅᵢ`. Its period follows,
//! and **period decides whether a feature is measurable at all** — a combination
//! with a 3,000-year period is not testable against a century of data, however
//! elegant. [`Commensurability::period`] exists so the basis can be filtered by
//! what the data can actually resolve, rather than by taste.

/// Integer coefficients on a set of angles, constrained to `Σ kᵢ = 0`.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Commensurability {
    pub k: Vec<i32>,
}

impl Commensurability {
    /// Construct, returning `None` if the d'Alembert constraint is violated or all
    /// coefficients are zero.
    pub fn new(k: Vec<i32>) -> Option<Self> {
        if k.iter().sum::<i32>() != 0 || k.iter().all(|&x| x == 0) {
            return None;
        }
        Some(Self { k })
    }

    /// Combined angle `Σ kᵢ θᵢ`, radians, unwrapped.
    ///
    /// Returns `None` if the angle count does not match the coefficient count.
    pub fn angle(&self, thetas: &[f64]) -> Option<f64> {
        if thetas.len() != self.k.len() {
            return None;
        }
        Some(
            self.k
                .iter()
                .zip(thetas)
                .map(|(&ki, &t)| ki as f64 * t)
                .sum(),
        )
    }

    /// Combined rate `Σ kᵢ ṅᵢ`, in the units of `rates`.
    pub fn frequency(&self, rates: &[f64]) -> Option<f64> {
        if rates.len() != self.k.len() {
            return None;
        }
        Some(
            self.k
                .iter()
                .zip(rates)
                .map(|(&ki, &r)| ki as f64 * r)
                .sum(),
        )
    }

    /// Period of the combination, in the reciprocal units of `rates`.
    ///
    /// `None` when the combined rate is zero — a degenerate combination with no
    /// time dependence.
    pub fn period(&self, rates: &[f64]) -> Option<f64> {
        let f = self.frequency(rates)?;
        (f != 0.0).then(|| (std::f64::consts::TAU / f).abs())
    }

    /// Order of the combination, `Σ |kᵢ| / 2`.
    ///
    /// Low order means a simple relationship — order 1 is a plain difference of two
    /// angles, the classical conjunction-opposition axis.
    pub fn order(&self) -> u32 {
        self.k.iter().map(|x| x.unsigned_abs()).sum::<u32>() / 2
    }

    /// Whether the coefficients share a common factor, making this a harmonic of a
    /// simpler combination.
    ///
    /// `(2, −2, 0)` is the second harmonic of `(1, −1, 0)`. Both are legitimate
    /// basis elements — they carry different information — but a caller enumerating
    /// *distinct relationships* usually wants only the primitive ones.
    pub fn is_primitive(&self) -> bool {
        fn gcd(a: u32, b: u32) -> u32 {
            if b == 0 {
                a
            } else {
                gcd(b, a % b)
            }
        }
        self.k
            .iter()
            .map(|x| x.unsigned_abs())
            .filter(|&x| x != 0)
            .fold(0u32, gcd)
            == 1
    }

    /// Canonical sign: the first non-zero coefficient is positive.
    ///
    /// `Φ` and `−Φ` carry identical information under cosine, so a basis should
    /// contain one of each pair, not both.
    pub fn canonical(mut self) -> Self {
        if let Some(&first) = self.k.iter().find(|&&x| x != 0) {
            if first < 0 {
                for x in self.k.iter_mut() {
                    *x = -*x;
                }
            }
        }
        self
    }
}

/// Enumerate all commensurabilities over `n` bodies with `|kᵢ| ≤ max_coeff`.
///
/// The d'Alembert constraint and sign canonicalisation are applied, so each
/// distinct relationship appears once.
pub fn enumerate(n: usize, max_coeff: i32) -> Vec<Commensurability> {
    let mut out = Vec::new();
    if n == 0 || max_coeff < 1 {
        return out;
    }
    let span = (2 * max_coeff + 1) as usize;
    let total = span.checked_pow(n as u32).unwrap_or(usize::MAX);
    let mut seen = std::collections::HashSet::new();

    for idx in 0..total {
        let mut rem = idx;
        let mut k = Vec::with_capacity(n);
        for _ in 0..n {
            k.push((rem % span) as i32 - max_coeff);
            rem /= span;
        }
        if let Some(c) = Commensurability::new(k) {
            let c = c.canonical();
            if seen.insert(c.k.clone()) {
                out.push(c);
            }
        }
    }
    out.sort_by_key(|c| (c.order(), c.k.clone()));
    out
}

/// Fourier features for a set of commensurabilities at one instant.
///
/// Emits `[cos Φ, sin Φ]` per combination, in enumeration order. Keeping the sine
/// term is what lets a downstream model learn a **phase offset**; dropping it
/// assumes the response peaks at exact alignment
/// (`docs/02-angular-encoding.md`).
pub fn features(combos: &[Commensurability], thetas: &[f64]) -> Vec<f64> {
    let mut out = Vec::with_capacity(combos.len() * 2);
    for c in combos {
        match c.angle(thetas) {
            Some(phi) => {
                let (s, co) = phi.sin_cos();
                out.push(co);
                out.push(s);
            }
            None => {
                out.push(f64::NAN);
                out.push(f64::NAN);
            }
        }
    }
    out
}

/// Mean motions of the planets, radians per day, for period filtering.
///
/// Sidereal, from standard orbital elements. Ordered Mercury through Neptune.
pub const PLANET_RATES: [f64; 8] = [
    0.071_42,   // Mercury,  87.969 d
    0.027_96,   // Venus,   224.701 d
    0.017_20,   // Earth,   365.256 d
    0.009_146,  // Mars,    686.980 d
    0.001_450,  // Jupiter, 4332.59 d
    0.000_583_9,// Saturn, 10759.22 d
    0.000_204_8,// Uranus, 30688.5 d
    0.000_104_4,// Neptune, 60182 d
];

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    #[test]
    fn rejects_combinations_violating_dalembert() {
        assert!(Commensurability::new(vec![1, 0, 0]).is_none());
        assert!(Commensurability::new(vec![2, -1, 0]).is_none());
        assert!(Commensurability::new(vec![0, 0, 0]).is_none());
        assert!(Commensurability::new(vec![1, -1, 0]).is_some());
        assert!(Commensurability::new(vec![2, -3, 1]).is_some());
    }

    #[test]
    fn combined_angle_is_rotation_invariant() {
        // The reason for the constraint: adding a constant to every input angle --
        // moving the origin of longitude -- must not change the feature.
        let c = Commensurability::new(vec![3, -5, 2]).unwrap();
        let base = [0.4, 1.9, 2.7];
        let a = c.angle(&base).unwrap();
        for shift in [0.1, 1.0, -2.5, TAU] {
            let moved: Vec<f64> = base.iter().map(|x| x + shift).collect();
            let b = c.angle(&moved).unwrap();
            assert!((a - b).abs() < 1e-12, "shift {shift} changed the angle");
        }
    }

    #[test]
    fn recovers_the_jupiter_saturn_synodic_period() {
        // k = (1, -1) on Jupiter and Saturn: the great conjunction cycle, 19.86 yr.
        let c = Commensurability::new(vec![1, -1]).unwrap();
        let rates = [PLANET_RATES[4], PLANET_RATES[5]];
        let p = c.period(&rates).unwrap() / 365.25;
        assert!((p - 19.86).abs() < 0.1, "got {p:.3} yr, expected 19.86");
    }

    #[test]
    fn recovers_the_venus_earth_synodic_period() {
        let c = Commensurability::new(vec![1, -1]).unwrap();
        let rates = [PLANET_RATES[1], PLANET_RATES[2]];
        let p = c.period(&rates).unwrap() / 365.25;
        assert!((p - 1.599).abs() < 0.01, "got {p:.4} yr, expected 1.599");
    }

    #[test]
    fn recovers_the_venus_earth_near_resonance() {
        // 8 Venus periods (8 x 224.701 = 1798 d... ) -- more precisely, 13 Earth
        // years are 4748 d and 8 Venus years 1798 d; the classical near-resonance
        // is 8 Earth years ~ 13 Venus years, so the combination whose frequency
        // nearly vanishes is 8*n_Venus - 13*n_Earth.
        //
        // The third slot carries a zero-rate reference direction (a slowly moving
        // perihelion longitude, idealised as fixed) purely to satisfy d'Alembert.
        let c = Commensurability::new(vec![8, -13, 5]).unwrap();
        let rates = [PLANET_RATES[1], PLANET_RATES[2], 0.0];
        let p = c.period(&rates).unwrap() / 365.25;
        assert!(
            (100.0..1000.0).contains(&p),
            "near-resonant beat should be centuries, got {p:.1} yr"
        );

        // The reversed coefficients are NOT resonant -- a good check that the
        // near-cancellation is real and not an artefact of large integers.
        let d = Commensurability::new(vec![13, -8, -5]).unwrap();
        let q = d.period(&rates).unwrap() / 365.25;
        assert!(q < 1.0, "reversed coefficients should be fast, got {q:.2} yr");
    }

    #[test]
    fn period_is_none_for_a_degenerate_combination() {
        let c = Commensurability::new(vec![1, -1]).unwrap();
        // Two bodies with identical rates never separate.
        assert!(c.period(&[0.01, 0.01]).is_none());
    }

    #[test]
    fn order_and_primitivity() {
        assert_eq!(Commensurability::new(vec![1, -1]).unwrap().order(), 1);
        assert_eq!(Commensurability::new(vec![2, -2]).unwrap().order(), 2);
        assert_eq!(Commensurability::new(vec![3, -1, -2]).unwrap().order(), 3);
        assert!(Commensurability::new(vec![1, -1]).unwrap().is_primitive());
        assert!(!Commensurability::new(vec![2, -2]).unwrap().is_primitive());
        assert!(Commensurability::new(vec![2, -3, 1]).unwrap().is_primitive());
    }

    #[test]
    fn canonical_form_removes_sign_duplicates() {
        let a = Commensurability::new(vec![-1, 1]).unwrap().canonical();
        let b = Commensurability::new(vec![1, -1]).unwrap().canonical();
        assert_eq!(a, b);
    }

    #[test]
    fn enumeration_is_deduplicated_and_constrained() {
        let all = enumerate(3, 2);
        assert!(!all.is_empty());
        // Every entry satisfies d'Alembert.
        assert!(all.iter().all(|c| c.k.iter().sum::<i32>() == 0));
        // No sign duplicates.
        let mut seen = std::collections::HashSet::new();
        for c in &all {
            assert!(seen.insert(c.k.clone()), "duplicate {:?}", c.k);
            let neg: Vec<i32> = c.k.iter().map(|x| -x).collect();
            assert!(!seen.contains(&neg), "sign duplicate of {:?}", c.k);
        }
        // Ordered by increasing complexity.
        assert!(all.windows(2).all(|w| w[0].order() <= w[1].order()));
    }

    #[test]
    fn constraint_collapses_the_search_space() {
        // Unconstrained, 4 bodies at |k| <= 3 gives 7^4 = 2401 combinations.
        // d'Alembert plus sign canonicalisation removes the great majority.
        let n = enumerate(4, 3).len();
        assert!(n < 2401 / 2, "expected a large reduction, got {n}");
        assert!(n > 20, "expected a usable basis, got {n}");
    }

    #[test]
    fn features_emit_cosine_and_sine_per_combination() {
        let combos = enumerate(2, 1);
        let f = features(&combos, &[0.3, 1.1]);
        assert_eq!(f.len(), combos.len() * 2);
        assert!(f.iter().all(|x| x.is_finite()));
        // Wrong arity yields NaN rather than a silent wrong answer.
        let bad = features(&combos, &[0.3]);
        assert!(bad.iter().all(|x| x.is_nan()));
    }
}

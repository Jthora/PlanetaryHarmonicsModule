//! Elastic response — converting tidal tensors to stress.
//!
//! Everything upstream of this module produces the **tide-generating potential's
//! second derivatives**, in s⁻². That is stress *shape*: correct in geometry and
//! timing, but carrying no information about how hard the body is squeezed. All
//! phase results are unaffected by the missing scale, which is why `docs/07` could
//! report C2–C4 without it. Anything involving *magnitude* needs this module.
//!
//! # ⚠ Scope: a calibrated scale factor, not an elastic solution
//!
//! A rigorous treatment integrates the radial eigenfunctions of an elastic,
//! self-gravitating, layered Earth (PREM) to get the stress tensor at depth.
//! **That is not what this does.**
//!
//! This applies the standard degree-2 surface relation used throughout the tidal
//! triggering literature: strain follows from the potential via the Love numbers
//! `h₂` and `l₂`, and stress follows from strain via Hooke's law with crustal
//! moduli. It is an order-unity-accurate scalar calibration, adequate for
//! estimating `Aσ₀` and locating `T_a`, and **not** adequate for a per-component
//! stress tensor at depth.
//!
//! Treat outputs as good to a factor of ~2. Where that matters, say so.
//!
//! # The chain
//!
//! For a degree-2 solid harmonic the potential goes as `r²`, so
//! `∂²V/∂r² = 2V/R²`, giving the potential back from the tensor's radial
//! component:
//!
//! ```text
//! V = R² · T_rr / 2
//! ```
//!
//! Areal strain at the surface from the degree-2 tide is
//!
//! ```text
//! ε_areal = (2h₂ − 6l₂) · V / (g R)
//! ```
//!
//! and stress follows as `σ ≈ 2μ ε` for a shear modulus `μ`. Combining,
//!
//! ```text
//! σ ≈ 2μ · (2h₂ − 6l₂) · R · T_rr / (4 g)
//! ```
//!
//! With IERS elastic values this predicts **~1 kPa** for the M2 solid Earth tide,
//! which is the published figure — the test below checks exactly that.

/// Elastic parameters for converting tidal tensors to stress.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Elastic {
    /// Degree-2 radial displacement Love number.
    pub h2: f64,
    /// Degree-2 horizontal displacement Love number (Shida number).
    pub l2: f64,
    /// Degree-2 potential Love number. Not used in the strain path; carried
    /// because ocean loading and gravimetric work need it.
    pub k2: f64,
    /// Body radius, m.
    pub radius_m: f64,
    /// Surface gravity, m/s².
    pub gravity: f64,
    /// Shear modulus, Pa.
    pub shear_modulus: f64,
}

impl Elastic {
    /// Earth, with IERS Conventions degree-2 elastic Love numbers and typical
    /// crustal shear modulus.
    pub const EARTH: Elastic = Elastic {
        h2: 0.6078,
        l2: 0.0847,
        k2: 0.2980,
        radius_m: 6_371_000.0,
        gravity: 9.80665,
        shear_modulus: 3.0e10,
    };

    /// Moon. Love numbers from lunar laser ranging; the Moon is far less
    /// deformable than Earth relative to its size.
    pub const MOON: Elastic = Elastic {
        h2: 0.0424,
        l2: 0.0107,
        k2: 0.02405,
        radius_m: 1_737_400.0,
        gravity: 1.62,
        shear_modulus: 6.0e10,
    };

    /// The dimensionless areal-strain combination `2h₂ − 6l₂`.
    pub fn strain_factor(&self) -> f64 {
        2.0 * self.h2 - 6.0 * self.l2
    }

    /// Conversion from tensor radial component (s⁻²) to stress (Pa).
    ///
    /// Multiply a tensor component by this to obtain a stress magnitude.
    pub fn stress_per_tensor(&self) -> f64 {
        2.0 * self.shear_modulus * self.strain_factor() * self.radius_m
            / (4.0 * self.gravity)
    }

    /// Strain corresponding to a tensor component, dimensionless.
    pub fn strain(&self, tensor_component: f64) -> f64 {
        self.strain_factor() * self.radius_m * tensor_component / (4.0 * self.gravity)
    }

    /// Stress in Pa corresponding to a tensor component in s⁻².
    pub fn stress(&self, tensor_component: f64) -> f64 {
        self.stress_per_tensor() * tensor_component
    }
}

/// Dieterich's characteristic time `t_a = Aσ₀ / τ̇`, in the units of `stressing_rate`.
///
/// With `a_sigma` in Pa and `stressing_rate` in Pa/year, the result is years.
pub fn characteristic_time(a_sigma: f64, stressing_rate: f64) -> f64 {
    a_sigma / stressing_rate
}

/// The critical period `T_a = 2π t_a` at which the seismicity response peaks
/// (Ader et al. 2014, eq. 8).
pub fn critical_period(a_sigma: f64, stressing_rate: f64) -> f64 {
    std::f64::consts::TAU * characteristic_time(a_sigma, stressing_rate)
}

/// Invert a measured critical period for the stressing rate it implies.
///
/// Useful in the direction the data actually runs: `T_a` is observable from a
/// response spectrum, `Aσ₀` is constrained by other work, and `τ̇` is what
/// forecasting wants.
pub fn stressing_rate_from(a_sigma: f64, critical_period: f64) -> f64 {
    std::f64::consts::TAU * a_sigma / critical_period
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Moon-on-Earth tensor scale: GM_moon / d³, in s⁻².
    const T_M2: f64 = 4.9028e12 / (3.844e8_f64 * 3.844e8 * 3.844e8);

    #[test]
    fn m2_solid_earth_tide_is_about_a_kilopascal() {
        let s = Elastic::EARTH.stress(T_M2).abs();
        assert!(
            (300.0..10_000.0).contains(&s),
            "M2 solid Earth tidal stress came out {s:.0} Pa; published values are \
             of order 1 kPa"
        );
    }

    #[test]
    fn m2_strain_is_of_order_1e_minus_8() {
        let e = Elastic::EARTH.strain(T_M2).abs();
        assert!(
            (1e-9..1e-7).contains(&e),
            "M2 tidal strain came out {e:e}; published values are ~1e-8"
        );
    }

    #[test]
    fn strain_factor_matches_the_iers_love_numbers() {
        // 2(0.6078) - 6(0.0847) = 0.7074
        assert!((Elastic::EARTH.strain_factor() - 0.7074).abs() < 1e-4);
    }

    #[test]
    fn moon_is_much_stiffer_relative_to_its_size() {
        // Lower Love numbers and smaller radius: far less stress per unit tensor.
        assert!(Elastic::MOON.stress_per_tensor() < Elastic::EARTH.stress_per_tensor());
    }

    #[test]
    fn stress_is_linear_in_the_tensor() {
        let e = Elastic::EARTH;
        assert!((e.stress(2.0 * T_M2) - 2.0 * e.stress(T_M2)).abs() < 1e-6);
        assert_eq!(e.stress(0.0), 0.0);
    }

    #[test]
    fn critical_period_round_trips() {
        // Parkfield-like: A*sigma0 = 6e-4 MPa = 600 Pa.
        let a_sigma = 600.0;
        let rate = 3000.0; // Pa/yr
        let t_a = critical_period(a_sigma, rate);
        assert!((stressing_rate_from(a_sigma, t_a) - rate).abs() < 1e-6);
    }

    #[test]
    fn ordinary_crust_gives_a_decadal_critical_period() {
        // A*sigma0 ~ 0.01-0.1 MPa, tectonic stressing ~3 kPa/yr.
        let lo = critical_period(1.0e4, 3000.0);
        let hi = critical_period(1.0e5, 3000.0);
        assert!(
            (15.0..40.0).contains(&lo) && (150.0..400.0).contains(&hi),
            "expected roughly 20-200 yr, got {lo:.0} to {hi:.0}"
        );
    }
}

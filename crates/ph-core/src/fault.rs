//! Resolving a tidal tensor onto a fault plane.
//!
//! # Conventions
//!
//! Fault geometry follows Aki & Richards, in a local **North-East-Down** frame:
//!
//! - **strike** φ — clockwise from north, 0–360°
//! - **dip** δ — from horizontal, 0–90°
//! - **rake** λ — slip direction within the plane, measured from strike, −180–180°
//!
//! ```text
//! normal  n̂ = (−sin δ sin φ,  sin δ cos φ, −cos δ)
//! slip    û = ( cos λ cos φ + cos δ sin λ sin φ,
//!               cos λ sin φ − cos δ sin λ cos φ,
//!              −sin λ sin δ )
//! ```
//!
//! # ⚠ Units: this is stress *shape*, not stress magnitude
//!
//! The tidal tensor here is the gradient of the tide-generating potential, in
//! s⁻². Converting to actual stress in Pa requires the body's elastic response —
//! Love numbers and a radial structure model (see IERS Conventions Ch. 7).
//!
//! **That scale factor is irrelevant to phase analysis**, which is what Phase 1
//! needs: the *timing* of Coulomb maxima is unchanged by a positive constant. It
//! is essential for any absolute ΔCFS in MPa, so do not report magnitudes from
//! this module without applying Love numbers first.

use crate::tidal::TidalTensor;

/// A fault plane in the Aki & Richards convention, degrees.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FaultPlane {
    pub strike_deg: f64,
    pub dip_deg: f64,
    pub rake_deg: f64,
}

impl FaultPlane {
    pub fn new(strike_deg: f64, dip_deg: f64, rake_deg: f64) -> Self {
        Self {
            strike_deg,
            dip_deg,
            rake_deg,
        }
    }

    /// Unit normal to the fault plane, in North-East-Down.
    pub fn normal(&self) -> [f64; 3] {
        let (sf, cf) = self.strike_deg.to_radians().sin_cos();
        let (sd, cd) = self.dip_deg.to_radians().sin_cos();
        [-sd * sf, sd * cf, -cd]
    }

    /// Unit slip direction of the hanging wall, in North-East-Down.
    pub fn slip(&self) -> [f64; 3] {
        let (sf, cf) = self.strike_deg.to_radians().sin_cos();
        let (sd, cd) = self.dip_deg.to_radians().sin_cos();
        let (sl, cl) = self.rake_deg.to_radians().sin_cos();
        [
            cl * cf + cd * sl * sf,
            cl * sf - cd * sl * cf,
            -sl * sd,
        ]
    }
}

/// Tractions resolved on a plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Resolved {
    /// Normal component, positive in tension (unclamping).
    pub normal: f64,
    /// Shear component along the slip direction.
    pub shear: f64,
    /// Magnitude of the full shear traction, regardless of direction.
    pub shear_magnitude: f64,
}

/// Resolve a tensor, already expressed in North-East-Down, onto a plane.
pub fn resolve(t: &TidalTensor, plane: &FaultPlane) -> Resolved {
    let n = plane.normal();
    let u = plane.slip();

    // Traction vector t_i = T_ij n_j.
    let mut trac = [0.0f64; 3];
    for i in 0..3 {
        for j in 0..3 {
            trac[i] += t.m[i][j] * n[j];
        }
    }

    let normal: f64 = (0..3).map(|i| trac[i] * n[i]).sum();
    let shear: f64 = (0..3).map(|i| trac[i] * u[i]).sum();
    let shear_vec: [f64; 3] = std::array::from_fn(|i| trac[i] - normal * n[i]);
    let shear_magnitude = shear_vec.iter().map(|x| x * x).sum::<f64>().sqrt();

    Resolved {
        normal,
        shear,
        shear_magnitude,
    }
}

/// Coulomb failure function `ΔCFS = τ + μ′ σₙ`, in the tensor's units.
///
/// `mu` is the effective friction coefficient. Cochran et al. (2004) found the
/// best fit for shallow thrust faults at **μ = 0.4**, with good correlation over
/// 0.2–0.6.
pub fn coulomb(t: &TidalTensor, plane: &FaultPlane, mu: f64) -> f64 {
    let r = resolve(t, plane);
    r.shear + mu * r.normal
}

/// Rotate a body-fixed tensor into the local North-East-Down frame at a point.
///
/// `lat_deg` and `lon_deg` are body-fixed (planetocentric) coordinates in the same
/// frame as the tensor — for the Moon, `MOON_PA`.
pub fn to_local_ned(t: &TidalTensor, lat_deg: f64, lon_deg: f64) -> TidalTensor {
    let (sp, cp) = lat_deg.to_radians().sin_cos();
    let (sl, cl) = lon_deg.to_radians().sin_cos();

    // Rows: North, East, Down. Down is the inward radial direction.
    let r = [
        [-sp * cl, -sp * sl, cp],
        [-sl, cl, 0.0],
        [-cp * cl, -cp * sl, -sp],
    ];

    // T' = R T Rᵀ
    let mut rt = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                rt[i][j] += r[i][k] * t.m[k][j];
            }
        }
    }
    let mut out = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                out[i][j] += rt[i][k] * r[j][k];
            }
        }
    }
    TidalTensor { m: out }
}

/// A regular grid over strike, dip and rake.
///
/// Weber, Bills & Johnson (2009) search exactly this space for deep moonquakes,
/// since Apollo's geometry cannot constrain fault parameters from first motions.
pub fn plane_grid(strike_step_deg: f64, dip_step_deg: f64, rake_step_deg: f64) -> Vec<FaultPlane> {
    let mut out = Vec::new();
    let mut strike = 0.0;
    while strike < 360.0 {
        let mut dip = dip_step_deg;
        while dip <= 90.0 {
            let mut rake = -180.0;
            while rake < 180.0 {
                out.push(FaultPlane::new(strike, dip, rake));
                rake += rake_step_deg;
            }
            dip += dip_step_deg;
        }
        strike += strike_step_deg;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
        (0..3).map(|i| a[i] * b[i]).sum()
    }

    fn norm(a: [f64; 3]) -> f64 {
        dot(a, a).sqrt()
    }

    #[test]
    fn normal_and_slip_are_unit_and_orthogonal() {
        for p in [
            FaultPlane::new(0.0, 90.0, 0.0),
            FaultPlane::new(37.0, 45.0, 90.0),
            FaultPlane::new(210.0, 15.0, -120.0),
        ] {
            assert!((norm(p.normal()) - 1.0).abs() < 1e-12);
            assert!((norm(p.slip()) - 1.0).abs() < 1e-12);
            // Slip lies in the plane, so it is perpendicular to the normal.
            assert!(dot(p.normal(), p.slip()).abs() < 1e-12);
        }
    }

    #[test]
    fn vertical_strike_slip_has_horizontal_normal() {
        // Strike north, dip 90 -> plane is the N-D plane, normal points east.
        let p = FaultPlane::new(0.0, 90.0, 0.0);
        let n = p.normal();
        assert!((n[1] - 1.0).abs() < 1e-12, "normal {n:?} should be +East");
        assert!(n[2].abs() < 1e-12, "normal should be horizontal");
    }

    #[test]
    fn traction_decomposition_is_consistent() {
        // |t|^2 == sigma_n^2 + |shear|^2 for any tensor and plane.
        let t = TidalTensor {
            m: [[2.0, 0.3, -0.4], [0.3, -1.0, 0.7], [-0.4, 0.7, -1.0]],
        };
        let p = FaultPlane::new(53.0, 37.0, 21.0);
        let r = resolve(&t, &p);

        let n = p.normal();
        let mut trac = [0.0f64; 3];
        for i in 0..3 {
            for j in 0..3 {
                trac[i] += t.m[i][j] * n[j];
            }
        }
        let total = dot(trac, trac);
        assert!(
            (total - (r.normal * r.normal + r.shear_magnitude * r.shear_magnitude)).abs() < 1e-12
        );
        // Shear along slip cannot exceed the full shear magnitude.
        assert!(r.shear.abs() <= r.shear_magnitude + 1e-12);
    }

    #[test]
    fn rake_rotation_projects_the_shear_traction() {
        // Varying rake only rotates the slip direction in the plane, so shear
        // along slip traces a cosine and peaks at the full shear magnitude.
        let t = TidalTensor {
            m: [[1.5, 0.2, 0.1], [0.2, -0.5, -0.3], [0.1, -0.3, -1.0]],
        };
        let mag = resolve(&t, &FaultPlane::new(20.0, 60.0, 0.0)).shear_magnitude;
        let best = (-180..180)
            .map(|r| resolve(&t, &FaultPlane::new(20.0, 60.0, r as f64)).shear)
            .fold(f64::MIN, f64::max);
        assert!((best - mag).abs() < 2e-3, "best {best} vs magnitude {mag}");
    }

    #[test]
    fn local_rotation_preserves_invariants() {
        let t = TidalTensor::from_body(4902.8, [1.0e5, -2.0e5, 3.0e5]);
        let l = to_local_ned(&t, 33.0, -117.0);
        // Trace and eigenvalues are rotation invariant.
        assert!((l.trace() - t.trace()).abs() < 1e-20);
        let (a, _) = t.eigen();
        let (b, _) = l.eigen();
        for i in 0..3 {
            assert!((a[i] - b[i]).abs() < 1e-20, "eigenvalue {i}: {a:?} vs {b:?}");
        }
    }

    #[test]
    fn down_axis_points_inward_at_the_equator() {
        // At lat 0, lon 0 the body-fixed +X is up, so Down is -X.
        let t = TidalTensor {
            m: [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        };
        let l = to_local_ned(&t, 0.0, 0.0);
        // The X-X component maps entirely onto Down-Down.
        assert!((l.m[2][2] - 1.0).abs() < 1e-12, "{:?}", l.m);
    }

    #[test]
    fn grid_covers_the_space() {
        let g = plane_grid(90.0, 30.0, 90.0);
        assert_eq!(g.len(), 4 * 3 * 4);
        assert!(g.iter().all(|p| p.dip_deg > 0.0 && p.dip_deg <= 90.0));
    }
}

/// The six independent components of a symmetric tensor, ordered
/// `(xx, yy, zz, xy, xz, yz)`.
pub fn components(t: &TidalTensor) -> [f64; 6] {
    [
        t.m[0][0], t.m[1][1], t.m[2][2],
        t.m[0][1], t.m[0][2], t.m[1][2],
    ]
}

/// Coefficients `c` such that `ΔCFS = c · components(T)`.
///
/// Coulomb stress is **linear in the tensor**:
///
/// ```text
/// ΔCFS = T_ij · C_ij      C_ij = ½(uᵢnⱼ + uⱼnᵢ) + μ nᵢnⱼ
/// ```
///
/// This is what makes a fault-orientation grid search cheap. Rather than
/// resolving every tensor against every candidate plane, precompute the
/// **covariance of the tensor components** over a time series once, then the
/// variance of ΔCFS for any plane is the quadratic form `cᵀ Σ c` — a few dozen
/// flops instead of a full pass over the series.
///
/// For a grid of thousands of planes over thousands of epochs that turns hours
/// into milliseconds, and it is exact rather than an approximation.
pub fn coulomb_coefficients(plane: &FaultPlane, mu: f64) -> [f64; 6] {
    let n = plane.normal();
    let u = plane.slip();
    let c = |i: usize, j: usize| 0.5 * (u[i] * n[j] + u[j] * n[i]) + mu * n[i] * n[j];
    // Off-diagonal terms appear twice in the symmetric contraction.
    [
        c(0, 0),
        c(1, 1),
        c(2, 2),
        2.0 * c(0, 1),
        2.0 * c(0, 2),
        2.0 * c(1, 2),
    ]
}

/// Covariance of tensor components over a set of tensors, as a symmetric 6×6.
///
/// Pair with [`coulomb_coefficients`]: `Var(ΔCFS) = cᵀ Σ c`.
pub fn component_covariance(tensors: &[TidalTensor]) -> [[f64; 6]; 6] {
    let n = tensors.len();
    let mut cov = [[0.0f64; 6]; 6];
    if n < 2 {
        return cov;
    }
    let vs: Vec<[f64; 6]> = tensors.iter().map(components).collect();
    let mut mean = [0.0f64; 6];
    for v in &vs {
        for k in 0..6 {
            mean[k] += v[k];
        }
    }
    for m in mean.iter_mut() {
        *m /= n as f64;
    }
    for v in &vs {
        for i in 0..6 {
            for j in 0..6 {
                cov[i][j] += (v[i] - mean[i]) * (v[j] - mean[j]);
            }
        }
    }
    let d = (n - 1) as f64;
    for row in cov.iter_mut() {
        for x in row.iter_mut() {
            *x /= d;
        }
    }
    cov
}

/// Standard deviation of ΔCFS implied by a component covariance, via `√(cᵀ Σ c)`.
pub fn coulomb_std(cov: &[[f64; 6]; 6], c: &[f64; 6]) -> f64 {
    let mut acc = 0.0;
    for i in 0..6 {
        for j in 0..6 {
            acc += c[i] * cov[i][j] * c[j];
        }
    }
    acc.max(0.0).sqrt()
}

#[cfg(test)]
mod linear_tests {
    use super::*;

    fn tensor(seed: f64) -> TidalTensor {
        TidalTensor::from_body(4902.8, [1.0e5 + seed, -2.0e5, 3.0e5 - 2.0 * seed])
    }

    #[test]
    fn coefficients_reproduce_the_direct_calculation() {
        for (s, d, r, mu) in [
            (0.0, 90.0, 0.0, 0.4),
            (37.0, 45.0, 90.0, 0.0),
            (210.0, 15.0, -120.0, 0.6),
        ] {
            let p = FaultPlane::new(s, d, r);
            let c = coulomb_coefficients(&p, mu);
            for k in 0..5 {
                let t = tensor(k as f64 * 1.0e4);
                let direct = coulomb(&t, &p, mu);
                let linear: f64 = c
                    .iter()
                    .zip(components(&t))
                    .map(|(a, b)| a * b)
                    .sum();
                assert!(
                    (direct - linear).abs() <= 1e-12 * direct.abs().max(1e-30),
                    "direct {direct:e} vs linear {linear:e}"
                );
            }
        }
    }

    #[test]
    fn quadratic_form_matches_the_sample_std() {
        let tensors: Vec<TidalTensor> =
            (0..200).map(|i| tensor(i as f64 * 500.0)).collect();
        let p = FaultPlane::new(53.0, 37.0, 21.0);
        let c = coulomb_coefficients(&p, 0.4);

        let cfs: Vec<f64> = tensors.iter().map(|t| coulomb(t, &p, 0.4)).collect();
        let mean = cfs.iter().sum::<f64>() / cfs.len() as f64;
        let direct = (cfs.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
            / (cfs.len() - 1) as f64)
            .sqrt();

        let via_cov = coulomb_std(&component_covariance(&tensors), &c);
        assert!(
            (direct - via_cov).abs() <= 1e-9 * direct,
            "direct {direct:e} vs covariance {via_cov:e}"
        );
    }

    #[test]
    fn covariance_of_too_few_tensors_is_zero() {
        assert_eq!(component_covariance(&[]), [[0.0; 6]; 6]);
    }
}

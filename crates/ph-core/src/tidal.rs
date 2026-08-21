//! Tide-generating potential, the tidal tensor, and concentration nodes.
//!
//! The degree-2 tide-generating potential from a body of mass `M` at geocentric
//! distance `d` produces, at Earth's centre, the trace-free tensor
//!
//! ```text
//! T_ij = (GM / d³) (3 n̂_i n̂_j − δ_ij)
//! ```
//!
//! with `n̂` the unit vector toward the body: `+2GM/d³` along the body axis,
//! `−GM/d³` transverse. Tensors from separate bodies **superpose linearly**, so
//! the complete global degree-2 field at an instant is five independent numbers.
//!
//! Units: positions in km, GM in km³/s², tensor components in s⁻².

/// A symmetric trace-free 3×3 tidal tensor, in s⁻².
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct TidalTensor {
    /// Row-major components. Symmetric: `m[i][j] == m[j][i]`.
    pub m: [[f64; 3]; 3],
}

impl TidalTensor {
    /// The degree-2 tensor from one body.
    ///
    /// `gm` is the body's gravitational parameter (km³/s²); `pos` is its
    /// geocentric position (km). Use **geometric** positions — aberration
    /// correction `None` — since tidal force acts on the instantaneous
    /// configuration (docs/10 §2).
    pub fn from_body(gm: f64, pos: [f64; 3]) -> Self {
        let d = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
        let n = [pos[0] / d, pos[1] / d, pos[2] / d];
        let k = gm / (d * d * d);

        let mut m = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                m[i][j] = k * (3.0 * n[i] * n[j] - delta);
            }
        }
        Self { m }
    }

    /// Sum of tensors — the combined field of several bodies.
    pub fn sum(tensors: impl IntoIterator<Item = TidalTensor>) -> Self {
        let mut acc = [[0.0f64; 3]; 3];
        for t in tensors {
            for i in 0..3 {
                for j in 0..3 {
                    acc[i][j] += t.m[i][j];
                }
            }
        }
        Self { m: acc }
    }

    /// Trace, which is analytically zero. Useful as a numerical check.
    pub fn trace(&self) -> f64 {
        self.m[0][0] + self.m[1][1] + self.m[2][2]
    }

    /// Eigenvalues and eigenvectors, ascending by eigenvalue.
    ///
    /// Cyclic Jacobi rotation. The tensor is symmetric, so this always converges.
    pub fn eigen(&self) -> ([f64; 3], [[f64; 3]; 3]) {
        let mut a = self.m;
        let mut v = [[0.0f64; 3]; 3];
        for i in 0..3 {
            v[i][i] = 1.0;
        }

        for _ in 0..64 {
            let off = a[0][1].abs() + a[0][2].abs() + a[1][2].abs();
            if off < 1e-300 {
                break;
            }
            for (p, q) in [(0usize, 1usize), (0, 2), (1, 2)] {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;

                for k in 0..3 {
                    let akp = a[k][p];
                    let akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..3 {
                    let apk = a[p][k];
                    let aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for k in 0..3 {
                    let vkp = v[k][p];
                    let vkq = v[k][q];
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }

        let mut vals = [a[0][0], a[1][1], a[2][2]];
        let mut idx = [0usize, 1, 2];
        idx.sort_by(|&i, &j| vals[i].partial_cmp(&vals[j]).unwrap());
        let sorted_vals = [vals[idx[0]], vals[idx[1]], vals[idx[2]]];
        let mut vecs = [[0.0f64; 3]; 3];
        for (c, &src) in idx.iter().enumerate() {
            for r in 0..3 {
                vecs[r][c] = v[r][src];
            }
        }
        vals = sorted_vals;
        (vals, vecs)
    }

    /// Unit vector along the principal (largest-eigenvalue) axis.
    ///
    /// This is the **tidal concentration node** of `docs/03-tidal-tensor.md` §2:
    /// where it pierces the surface gives two antipodal concentration points.
    /// Sign is arbitrary — the axis is a direction, not an orientation.
    pub fn principal_axis(&self) -> [f64; 3] {
        let (_, vecs) = self.eigen();
        [vecs[0][2], vecs[1][2], vecs[2][2]]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GM_MOON: f64 = 4902.800118;
    const D_MOON: f64 = 384_400.0;

    #[test]
    fn single_body_has_expected_eigenvalues() {
        // Body on +z: eigenvalues (−k, −k, +2k) with k = GM/d³.
        let t = TidalTensor::from_body(GM_MOON, [0.0, 0.0, D_MOON]);
        let k = GM_MOON / D_MOON.powi(3);
        let (vals, _) = t.eigen();
        assert!((vals[0] + k).abs() < 1e-18, "got {:?}", vals);
        assert!((vals[1] + k).abs() < 1e-18, "got {:?}", vals);
        assert!((vals[2] - 2.0 * k).abs() < 1e-18, "got {:?}", vals);
    }

    #[test]
    fn tensor_is_trace_free() {
        let t = TidalTensor::from_body(GM_MOON, [1.0e5, 2.0e5, 3.0e5]);
        assert!(t.trace().abs() < 1e-20);
    }

    #[test]
    fn principal_axis_points_along_the_body() {
        let t = TidalTensor::from_body(GM_MOON, [0.0, 0.0, D_MOON]);
        let ax = t.principal_axis();
        assert!(ax[2].abs() > 0.999_999, "axis {:?} should be ±z", ax);
    }

    #[test]
    fn superposition_is_linear() {
        let a = TidalTensor::from_body(GM_MOON, [D_MOON, 0.0, 0.0]);
        let b = TidalTensor::from_body(GM_MOON, [0.0, D_MOON, 0.0]);
        let s = TidalTensor::sum([a, b]);
        for i in 0..3 {
            for j in 0..3 {
                assert!((s.m[i][j] - (a.m[i][j] + b.m[i][j])).abs() < 1e-24);
            }
        }
        assert!(s.trace().abs() < 1e-20);
    }
}

//! Batched geometric states, wrapping `rustspice-core`.
//!
//! Two conventions matter here and both are easy to get wrong:
//!
//! - **Aberration correction is `None`.** Tidal force acts on the instantaneous
//!   geometric configuration; light-time correction would bias it with no physical
//!   justification (`docs/10-rustspice-requirements.md` §2).
//! - **Time-scale correctness dominates ephemeris precision.** For M2, 1° of tidal
//!   phase is 124 s, so a UTC/TAI mixup is a systematic 0.3° error — while DE440
//!   is four orders of magnitude more accurate than we need.

use rustspice_core::{Aberration, Et, Session, Vec3};

/// Bodies contributing to the tidal tensor, in decreasing order of effect.
///
/// Moon and Sun dominate; planetary terms are ~10⁻⁵ of lunar and are carried for
/// the generalised expansion, not because they are expected to matter directly
/// (`docs/03-tidal-tensor.md` §4).
pub const TIDAL_BODIES: &[&str] = &[
    "MOON",
    "SUN",
    "VENUS BARYCENTER",
    "JUPITER BARYCENTER",
    "MARS BARYCENTER",
    "MERCURY BARYCENTER",
    "SATURN BARYCENTER",
];

/// A body's gravitational parameter, km³/s², from the loaded PCK.
///
/// Requires a GM kernel such as `gm_de440.tpc` — the tensor is GM/d³, so this is
/// not optional.
pub fn gm(session: &mut Session, body: &str) -> rustspice_core::Result<f64> {
    Ok(session.constant(body, "GM")?[0])
}

/// Geometric positions of one body over a series of epochs, observer-centred.
pub fn positions_at(
    session: &mut Session,
    target: &str,
    epochs: &[Et],
    frame: &str,
    observer: &str,
) -> rustspice_core::Result<Vec<Vec3>> {
    session.positions(target, epochs, frame, observer, Aberration::None)
}

/// Evenly spaced epochs, `count` samples at `step_s` seconds from `start`.
pub fn epoch_series(start: Et, step_s: f64, count: usize) -> Vec<Et> {
    (0..count).map(|i| start.offset(step_s * i as f64)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn epoch_series_is_evenly_spaced() {
        let e = epoch_series(Et(0.0), 3600.0, 4);
        assert_eq!(e.len(), 4);
        assert!((e[3].0 - 10800.0).abs() < 1e-9);
    }

    #[test]
    fn tidal_bodies_lead_with_moon_and_sun() {
        assert_eq!(TIDAL_BODIES[0], "MOON");
        assert_eq!(TIDAL_BODIES[1], "SUN");
    }
}

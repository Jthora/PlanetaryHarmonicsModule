//! Tidal fields from real ephemeris geometry.
//!
//! Bridges [`crate::ephemeris`] and [`crate::tidal`]: given a body to be stressed
//! and the bodies stressing it, compute the combined tidal tensor at any epoch.
//!
//! The same type serves both validation phases, only the arguments change:
//!
//! | Stressed body | Frame | Stressing bodies |
//! |---|---|---|
//! | Moon (Phase 1) | `MOON_PA` | Earth, Sun |
//! | Earth (Phase 3) | `ITRF93` | Moon, Sun, planets |
//!
//! Positions are **geometric** — aberration correction `None`. Tidal force acts on
//! the instantaneous configuration, so light-time correction would bias it with no
//! physical justification.

use crate::tidal::TidalTensor;
use rustspice_core::{Aberration, Et, Result, Session, Vec3};

/// Bodies raising tides on the Moon. Earth dominates; the Sun contributes ~1/2%.
pub const LUNAR_TIDE_BODIES: &[&str] = &["EARTH", "SUN"];

/// Bodies raising tides on Earth, in decreasing order of effect.
pub const TERRESTRIAL_TIDE_BODIES: &[&str] = &[
    "MOON",
    "SUN",
    "VENUS BARYCENTER",
    "JUPITER BARYCENTER",
    "MARS BARYCENTER",
    "MERCURY BARYCENTER",
    "SATURN BARYCENTER",
];

/// A configured tidal field: which body is stressed, by what, in which frame.
pub struct TidalField {
    observer: String,
    frame: String,
    /// Body name paired with its GM, km³/s². Resolved once at construction.
    bodies: Vec<(String, f64)>,
}

impl TidalField {
    /// Resolve GM for each stressing body and configure the field.
    ///
    /// Requires a GM kernel (`gm_de440.tpc`) to be loaded — the tensor is GM/d³,
    /// so this is not optional.
    pub fn new(
        session: &mut Session,
        observer: &str,
        frame: &str,
        bodies: &[&str],
    ) -> Result<Self> {
        let mut resolved = Vec::with_capacity(bodies.len());
        for b in bodies {
            let gm = session.constant(b, "GM")?[0];
            resolved.push((b.to_string(), gm));
        }
        Ok(Self {
            observer: observer.to_string(),
            frame: frame.to_string(),
            bodies: resolved,
        })
    }

    /// The field on the Moon, in the lunar principal-axes frame.
    ///
    /// `MOON_PA` is libration-aware, so a fixed surface or interior location keeps
    /// fixed coordinates — necessary for per-nest analysis.
    pub fn on_moon(session: &mut Session) -> Result<Self> {
        Self::new(session, "MOON", "MOON_PA", LUNAR_TIDE_BODIES)
    }

    /// The field on Earth, in an Earth-fixed frame.
    pub fn on_earth(session: &mut Session, frame: &str) -> Result<Self> {
        Self::new(session, "EARTH", frame, TERRESTRIAL_TIDE_BODIES)
    }

    /// Combined tidal tensor at one epoch.
    pub fn tensor_at(&self, session: &mut Session, et: Et) -> Result<TidalTensor> {
        let mut parts = Vec::with_capacity(self.bodies.len());
        for (name, gm) in &self.bodies {
            let p = session.position(name, et, &self.frame, &self.observer, Aberration::None)?;
            parts.push(TidalTensor::from_body(*gm, [p.x, p.y, p.z]));
        }
        Ok(TidalTensor::sum(parts))
    }

    /// Combined tidal tensors over a series of epochs.
    ///
    /// Queries each body once across all epochs rather than once per epoch, which
    /// keeps the SPICE segment search warm.
    pub fn tensors(&self, session: &mut Session, epochs: &[Et]) -> Result<Vec<TidalTensor>> {
        let mut acc: Vec<TidalTensor> = vec![TidalTensor::default(); epochs.len()];
        for (name, gm) in &self.bodies {
            let ps: Vec<Vec3> =
                session.positions(name, epochs, &self.frame, &self.observer, Aberration::None)?;
            for (a, p) in acc.iter_mut().zip(ps) {
                let t = TidalTensor::from_body(*gm, [p.x, p.y, p.z]);
                for i in 0..3 {
                    for j in 0..3 {
                        a.m[i][j] += t.m[i][j];
                    }
                }
            }
        }
        Ok(acc)
    }

    /// Bodies and their gravitational parameters, as resolved.
    pub fn bodies(&self) -> &[(String, f64)] {
        &self.bodies
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn body_lists_are_ordered_by_effect() {
        assert_eq!(LUNAR_TIDE_BODIES, &["EARTH", "SUN"]);
        assert_eq!(TERRESTRIAL_TIDE_BODIES[0], "MOON");
        assert_eq!(TERRESTRIAL_TIDE_BODIES[1], "SUN");
    }
}

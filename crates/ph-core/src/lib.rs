//! Tidal tensors, harmonic decomposition, and phase statistics.
//!
//! This crate is the scientific computation layer of PlanetaryHarmonicsModule.
//! It produces numbers with documented units and provenance. It does not
//! interpret them — symbolic mappings live upstream in `AstrologyCore`.
//!
//! Ephemerides come from [`rustspice_core`], consumed as a Rust library. There is
//! no WASM boundary within this chain; see `docs/10-rustspice-requirements.md`.
//!
//! # Layout
//!
//! - [`tidal`] — tide-generating potential, the tidal tensor, concentration nodes
//! - [`harmonics`] — Fourier angular encodings and least-squares decomposition
//! - [`phase`] — tidal phase from a sampled quasi-periodic forcing
//! - [`demod`] — complex demodulation; isolating one constituent for R(ω)
//! - [`doodson`] — analytic constituent phases from the fundamental arguments
//! - [`stats`] — the generalised Schuster test and time-shifted null distributions
//! - [`ephemeris`] — batched geometric states, wrapping `rustspice-core`
//! - [`field`] — tidal fields from real ephemeris geometry
//! - [`fault`] — resolving a tensor onto a fault plane; Coulomb failure stress
//! - [`love`] — elastic response; tensors to stress in Pa, and `T_a`
//! - [`catalog`] — event catalogues
//! - [`apollo`] — Apollo PSE catalogue ingestion (Phase 1 validation)
//! - [`parkfield`] — Parkfield LFE catalogue ingestion (Phase 2 testbed)
//!
//! # A standing invariant
//!
//! Oscillatory stress does not change the mean event rate — only the timing
//! (Heimisson & Avouac 2020, eq. 6: `⟨R⟩ = r` exactly). Tides redistribute *when*
//! events occur; they do not create them. Any output implying otherwise is a bug.

pub mod apollo;
pub mod catalog;
pub mod demod;
pub mod doodson;
pub mod ephemeris;
pub mod fault;
pub mod field;
pub mod harmonics;
pub mod love;
pub mod parkfield;
pub mod phase;
pub mod stats;
pub mod tidal;

pub use rustspice_core;

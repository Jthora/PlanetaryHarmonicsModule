//! Tidal tensors, harmonic decomposition, and phase statistics.
//!
//! This crate is the scientific computation layer of PlanetaryHarmonicsModule.
//! It produces numbers with documented units and provenance. It does not
//! interpret them — symbolic mappings live upstream in `AstrologyCore`.
//!
//! Ephemerides come from [`rustspice_core`], consumed as a Rust library. There is
//! no WASM boundary within this chain; see `docs/10-rustspice-requirements.md`.
//!
//! # Scope
//!
//! This crate holds only what **more than one** downstream project needs. Catalogue
//! parsers, research examples and domain analyses live with the application that
//! needs them — the seismology programme is in
//! [EarthquakeForecastModule](https://github.com/Jthora/EarthquakeForecastModule).
//! A module only one consumer needs is application code, not library code.
//!
//! # Layout
//!
//! - [`tidal`] — tide-generating potential, the tidal tensor, concentration nodes
//! - [`harmonics`] — Fourier angular encodings and least-squares decomposition
//! - [`chart`] — body states in geocentric, heliocentric and barycentric frames
//! - [`chart_features`] — aspects, declinations, shape and resonance from a chart
//! - [`commensurability`] — multi-body angular relationships under d'Alembert
//! - [`harmonic_model`] — precompute a harmonic ephemeris; O(1) timestream queries
//! - [`phase`] — tidal phase from a sampled quasi-periodic forcing
//! - [`demod`] — complex demodulation; isolating one constituent for R(ω)
//! - [`doodson`] — analytic constituent phases from the fundamental arguments
//! - [`events`] — angle-domain event finding; solve for crossings, don't sample
//! - [`stats`] — the generalised Schuster test and time-shifted null distributions
//! - [`ephemeris`] — batched geometric states, wrapping `rustspice-core`
//! - [`field`] — tidal fields from real ephemeris geometry
//! - [`fault`] — resolving a tensor onto a fault plane; Coulomb failure stress
//! - [`love`] — elastic response; tensors to stress in Pa, and `T_a`
//! - [`catalog`] — event types; catalogue *ingestion* lives with its application
//!
//! # A standing invariant
//!
//! Oscillatory stress does not change the mean event rate — only the timing
//! (Heimisson & Avouac 2020, eq. 6: `⟨R⟩ = r` exactly). Tides redistribute *when*
//! events occur; they do not create them. Any output implying otherwise is a bug.

pub mod catalog;
pub mod chart;
pub mod chart_features;
pub mod chart_cycles;
pub mod chart_local;
pub mod commensurability;
pub mod demod;
pub mod doodson;
pub mod ephemeris;
pub mod events;
pub mod fault;
pub mod field;
pub mod harmonic_model;
pub mod harmonics;
pub mod love;
pub mod phase;
pub mod stats;
pub mod tidal;

pub use rustspice_core;

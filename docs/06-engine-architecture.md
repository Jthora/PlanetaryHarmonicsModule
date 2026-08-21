# 06 — Engine Architecture

Rust core, WASM + TypeScript surface. Design notes for the computation layer.

---

## 1. Spherical-harmonic precompute — the key to global bulk evaluation

The tidal field is a **low-degree spherical harmonic field**: degree 2 captures
nearly everything, plus a small lunar degree-3 term (doc 03).

Therefore:

```
per timestep:   compute time-varying SH coefficients      O(bodies)   ~5 numbers
per location:   synthesise from coefficients              O(coeffs)
```

Global evaluation becomes **O(coefficients), not O(locations × bodies)**. Cost
per location collapses to a handful of flops rather than a full ephemeris
evaluation.

This is the actual answer to "global scale, not just one point." Use **HEALPix**
for the mesh — equal-area, and natively compatible with SH synthesis.

## 2. Root-find in the angle domain, don't sample time

For Resonant Finder and Star Seer: **do not sample the timestream densely and hunt
for peaks.**

Angular combinations Φ_k(t) = Σ kᵢθᵢ(t) are near-linear in time plus small
periodic corrections. So crossing times can be **predicted analytically and
Newton-refined** to machine precision:

```
1. linear estimate from mean motions          →  t₀
2. Newton refine on Φ_k(t) − target            →  t*
3. converges in 2–3 iterations
```

Event detection goes from **O(sample rate)** to **O(number of events)**. Orders of
magnitude faster, and it yields exact millisecond event times — precisely Star
Seer's requirement, which dense sampling cannot deliver at any reasonable cost.

## 3. Harmonic ephemeris — O(1) timestream queries

Because these combinations have near-constant frequencies, precompute a table of
(frequency, phase, amplitude) triples once, then evaluate any epoch in closed
form. Effectively Fourier-transform the system once; queries become O(1).

This is exactly what the tidal community already does with harmonic constituents —
a validated pattern, not a novel gamble.

## 4. WASM boundary discipline

**Boundary crossings dominate cost.** Design accordingly:

- **Columnar batch API only.** Arrays of times in → typed arrays of results out.
  Never expose a scalar per-call interface.
- Return `Float64Array` views over WASM linear memory; avoid per-call allocation.
- Use **SIMD128** (2×f64 or 4×f32 lanes) for the inner evaluation loops.
- Evaluate Chebyshev ephemeris polynomials via **Clenshaw recurrence**, batched
  across time — this is what SPICE SPK segments store natively.
- Compute cos/sin harmonic ladders by **angle-sum recurrence** from n = 1 rather
  than N transcendental calls (doc 02).

## 5. RustSPICE integration

Vendored at `modules/RustSPICE` (submodule → `github.com/Jthora/RustSPICE`).

Currently available, per its own docs: time system, coordinate/frame system,
CSPICE conversion strategy, WASM build pipeline (`wasm-pack-build.sh`,
`build-cspice-wasm.sh`), TypeScript integration layer, benchmarks, kernels.

Branches to be aware of: `main` (= tag `archive/from-scratch-attempt`),
`pivot/wasm-ts-layer`, tag `v1.0.0`. **Confirm which branch is canonical before
pinning** — the submodule currently tracks the `main` commit.

Boundary: RustSPICE supplies **ephemeris, time transforms, and reference frames**.
PlanetaryHarmonics consumes those and produces **tidal tensors, harmonic
decompositions, and feature vectors**. Keep the seam clean — no astrology-specific
or feature-specific logic inside RustSPICE.

## 6. Output metadata — required, not optional

Every emitted feature vector carries provenance:

- Reference frame and epoch system (ecliptic of date vs. J2000, TT vs. TDB vs. UTC)
- Ephemeris kernel identifiers and versions
- Feature tier (A / B / C per doc 00)
- Units, explicitly

Feature provenance is what makes results reproducible and reviewable. A feature
whose frame is ambiguous is worthless in a statistical test.

## 7. Precision notes

- **Time scales matter at this precision.** TT/TDB/UTC confusion produces
  systematic phase errors that will masquerade as signal in harmonic analysis.
  Leap seconds are not optional.
- Long-period constituents (Ssa, Sa, 18.61 yr node) require **long, stable time
  baselines** — accumulated floating-point drift over a multi-decade timestream is
  a real risk. Prefer two-part Julian date representations (high/low) throughout.
- The 18.61-year nodal cycle needs ~40+ years of data to constrain. The USGS
  catalogue supports this; be careful that Mc changes over that span do not
  alias into it (doc 04).

## 8. Repository layout (proposed)

```
src/
  ephemeris/     RustSPICE bindings, batch evaluation
  tidal/         tensor assembly, CFS projection, ocean loading
  harmonics/     Doodson expansion, Fourier encodings, commensurabilities
  spatial/       HEALPix mesh, SH synthesis
  events/        angle-domain root finding (Star Seer, Resonant Finder)
  features/      feature vector assembly + provenance metadata
  wasm/          boundary layer, columnar API
ts/              TypeScript surface
docs/            this directory
modules/
  RustSPICE/     submodule
```

# 10 — Data and API Requirements for RustSPICE

What PlanetaryHarmonicsModule needs from `modules/RustSPICE`, derived from the
research in docs 01–09. Written for the RustSPICE maintainers.

**Summary in one line:** we need **batched geometric states in several frames, with
rigorous time-scale handling** — and the *time scales matter far more than the
ephemeris precision*.

---

## 1. The headline: precision is not the constraint, correctness is

Worth stating up front so effort goes to the right place.

**Ephemeris precision is massively over-specified for our use.** The tidal
potential scales as `GM·R²/d³`, so a relative distance error ε produces a 3ε
amplitude error. DE440 lunar positions are accurate to metres — relative error
~10⁻⁸, giving a tidal amplitude error ~3×10⁻⁸. We need perhaps 10⁻⁴. We have four
orders of magnitude of headroom.

**Time-scale handling is the real risk.** For the M2 constituent (12.42 h), 1° of
tidal phase is 124 seconds. So:

- Confusing UTC with TAI (currently 37 s offset) → **0.3° M2 phase error**, small
  but *systematic* and it will alias into every harmonic analysis
- Missing or wrong leap seconds → same class of error, and it varies over the
  catalogue's 130-year span, which is far worse than a constant offset
- Sub-second timing accuracy is **not** needed

**Please prioritise leap-second and time-scale correctness over ephemeris
accuracy.** A silently wrong TAI–UTC is the single most likely way this project
gets a fake result.

---

## 2. Geometric states, NOT light-time corrected ⚠

**We need `NONE` aberration correction — geometric positions.**

Tidal force depends on the instantaneous geometric configuration. The standard
`LT+S` default would be wrong here: the Sun's light-time is 8.3 minutes, which is
~0.03° of S2 phase — small, but a systematic bias with no physical justification
in a dynamical calculation.

If the API defaults to `LT+S`, we will need an explicit way to request `NONE`, and
ideally the default should be reconsidered for dynamics use.

---

## 3. Bodies required

| Body | NAIF ID | Purpose |
|---|---|---|
| Moon | 301 | Dominant tidal term |
| Sun | 10 | Second tidal term |
| Mercury–Neptune | 199, 299, 499, 599, 699, 799, 899 | Generalised Doodson expansion (doc 02) |
| Earth | 399 | Observer |
| Earth–Moon barycentre | 3 | Intermediate computations |
| Pluto | 999 | Optional, completeness only |

Use planet **barycentres** (1–8) where the ephemeris provides them more accurately;
for tidal purposes barycentre vs. body centre is irrelevant given point 1.

---

## 4. State vectors

**Positions and velocities.** Velocities are required — `dΔCFS/dt` is central to
the rate-and-state analysis (docs 03, 09) and finite-differencing positions across
the WASM boundary would be wasteful and less accurate.

Observer: **geocentric** primarily; barycentric for some intermediate work.

---

## 5. Reference frames

| Frame | Purpose | Notes |
|---|---|---|
| `J2000` / ICRF | Base frame | |
| **Ecliptic of date** | Angular/harmonic features (doc 02) | Needs precession; "of date" not J2000 |
| **ITRF93** or high-precision `IAU_EARTH` | Surface locations, sub-body point tracks (doc 03) | Requires precession, nutation, Earth rotation, polar motion |
| **`MOON_PA`** (principal axes) | **Deep moonquake work** (doc 05) | Libration-aware; needed for Phase 1 validation |

The Earth-fixed transform is the demanding one — it needs full Earth orientation,
not a simple sidereal-time approximation. Polar motion at the ~0.3 arcsec level is
below our threshold, but nutation is not.

**`MOON_PA` matters more than it might appear.** Deep moonquakes are the project's
Phase 1 validation dataset (doc 05 §1) — the known-answer test that proves the
pipeline works before any terrestrial claim. That work needs selenographic
coordinates with proper libration.

---

## 6. Time systems

Required conversions, both directions:

```
UTC  ↔  TAI  ↔  TT  ↔  TDB
ISO 8601 / calendar strings  →  internal epoch
```

- **Leap seconds** via LSK, covering **1900–2030** minimum
- **Two-part Julian date** (high/low) or equivalent internal representation —
  single-f64 JD loses precision over a 130-year baseline, and the 18.61-year nodal
  constituent needs a stable long baseline (doc 06 §7)
- Earthquake catalogues are **UTC**; ephemerides want **TDB**. This conversion is
  on the hot path for every event.

---

## 7. Kernels to bundle or document

| Kernel | Purpose |
|---|---|
| `de440.bsp` (or `de441` for extended range) | Planetary + lunar ephemeris |
| `naif0012.tls` | Leap seconds |
| `pck00011.tpc` | Body constants — radii, orientation |
| **`gm_de440.tpc`** | **GM values — required for tidal amplitudes** |
| `earth_latest_high_prec.bpc` | Earth orientation → ITRF |
| `moon_pa_de440_*.bpc` + `moon_de440_*.tf` | Lunar frames for moonquake work |

**`gm_de440.tpc` is easy to overlook** but the tidal tensor is `GM/d³` — without
correct GM values there is no calculation.

**Coverage:** DE440 spans 1550–2650, which covers USGS (1900–present) and Apollo
(1969–1977) comfortably. DE441 only if we extend to deep historical catalogues.

---

## 8. API shape — the most important section

### Columnar batch evaluation, mandatory

```rust
// what we need
fn states_batch(
    target: NaifId,
    epochs_tdb: &[f64],      // N epochs
    frame: Frame,
    observer: NaifId,
    aberration: Correction,  // we pass NONE
) -> Result<StateBatch, SpiceError>;   // 6N contiguous f64
```

**Not** a scalar `state(target, epoch, ...)` called N times from TypeScript.
WASM boundary crossings dominate cost — a per-call interface would make the
boundary, not the computation, the bottleneck.

Ideally also multi-target: `targets: &[NaifId]` × `epochs: &[f64]` → one call
returning a contiguous block, since we always need all bodies at each epoch.

### Zero-copy returns

Return a pointer + length into WASM linear memory that TypeScript wraps as a
`Float64Array` view. No per-call allocation, no serialisation.

### Kernel loading from bytes

Browser has no filesystem. We need:

```rust
fn load_kernel_bytes(name: &str, data: &[u8]) -> Result<(), SpiceError>;
```

so kernels can be fetched as `ArrayBuffer` and handed in directly.

### No panics across the boundary

Every fallible operation returns `Result`. A Rust panic in WASM aborts the module
and cannot be recovered by the host — in a long batch run over millions of epochs
that means losing the whole job to one bad input.

### Reentrancy

If we parallelise (web workers, or native multi-threaded runs), kernel state must
either be shareable read-only across threads or cleanly instanced per thread.
Global mutable kernel pools would force us to serialise.

---

## 9. Throughput targets — smaller than you might expect

Worth stating so the layer is not over-engineered.

The spherical-harmonic precompute (doc 06 §1) means **ephemeris is needed per
epoch, not per surface location** — the global field comes from ~5 coefficients
synthesised onto the mesh. So:

| Workload | Epochs | × bodies | Total states |
|---|---|---|---|
| Global hourly, 1900–2030 | 1.14 × 10⁶ | 11 | **1.3 × 10⁷** |
| Per-event, full USGS catalogue | ~3 × 10⁶ | 11 | **3.3 × 10⁷** |
| Apollo moonquake catalogue | ~10⁴ | 3 | 3 × 10⁴ |

**Order 10⁷ state evaluations for a full production run.** At 1 µs per state that
is ~30 seconds. Ephemeris evaluation is comfortably *not* the bottleneck — the
ocean tidal loading convolution is (doc 08 §8).

So: **correctness and batch API shape matter; extreme micro-optimisation does
not.** SIMD and Clenshaw batching are nice-to-have, not critical path.

**One exception — Star Seer** (doc 06 §2) has a different profile: few epochs,
latency-sensitive, real-time. A low-latency single-epoch path is worth keeping
alongside the batch API, but the batch path is what production runs need.

---

## 10. What we do NOT need

Scope-limiting, to keep the RustSPICE surface small:

- ❌ Light-time / stellar aberration corrections (see §2 — we want them *off*)
- ❌ Spacecraft ephemerides, CK/instrument kernels, FOV geometry
- ❌ Occultation, eclipse, or ray-surface intercept geometry
- ❌ DSK / shape models
- ❌ Star catalogues
- ❌ Any astrology-specific or feature-specific logic — that belongs in
  PlanetaryHarmonics, not RustSPICE. Keep the seam clean (doc 06 §5).

---

## 11. Priority order

```
1  Time systems + leap seconds (UTC/TAI/TT/TDB), two-part JD    highest risk
2  Batched geometric states, geocentric, J2000                  core
3  Kernel loading from bytes; Result-based error handling        blocks browser use
4  Ecliptic-of-date frame                                       harmonic features
5  Earth-fixed frame with full orientation                      surface locations
6  MOON_PA frame                                                Phase 1 validation
7  Velocities                                                   dCFS/dt
8  Zero-copy returns, multi-target batching                      performance
9  SIMD / Clenshaw micro-optimisation                            nice to have
```

Items 1–3 unblock essentially all of our work. Item 6 is needed earlier than its
position suggests, because deep moonquakes are the first validation phase.

---

## 12. Open questions for the RustSPICE maintainers

1. Which branch is canonical — `main`, `pivot/wasm-ts-layer`, or `v1.0.0`? The
   submodule currently pins `main`; we will repoint on request.
2. Is high-precision Earth orientation (binary PCK / ITRF) in scope, or should we
   implement the ICRF→ITRF transform ourselves from IERS EOP series?
3. Are lunar frames (`MOON_PA`) supported or planned?
4. What is the current kernel-loading model — filesystem, embedded, or bytes?
5. Is the existing API scalar or batch? If scalar, is a columnar path feasible?

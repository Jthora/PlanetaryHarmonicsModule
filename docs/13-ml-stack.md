# 13 — ML Stack: Layers, Ladders, and What Actually Earns a Model

Response to the proposal of one ML layer per stage of the dependency chain —
ephemeris/RustSPICE → PlanetaryHarmonics → Earthquake Finder.

**Verdict: the structure is right, but the three layers are not equally deserving
of ML, and one of them barely deserves any.**

---

## The governing principle

ML earns its place where there is a **gap between what can be computed and what is
needed**. No gap, no model. Applied to the chain, the layers differ along one axis:

> **How much of the answer does physics already give you?**

| Layer | Physics gives | ML's real role | Validation quality |
|---|---|---|---|
| **Ephemeris** | ~everything | almost none — build algorithms | **exact** |
| **Harmonics** | the forward model | surrogates + inverse problems | **computable ground truth** |
| **Earthquake** | almost nothing | irreducible statistical inference | **noisy, contested** |

Two consequences fall straight out:

1. **ML's share of the work grows exactly as physics' share shrinks.**
2. **Validation quality degrades along the same axis.**

Which sets the build order. Each layer's validation is cleaner than the next's, so
build outward from the clean end. That is not a scheduling preference — it is the
only order in which a failure can be attributed to the right layer.

---

## Layer 1 — Ephemeris / RustSPICE: **do not build ML here**

This is the critique. The instinct to put a "next-gen ML" at the ephemeris layer is
the weakest part of the proposal, and it is worth being blunt about why.

**There is no gap.** Ephemerides are a solved problem:

- DE440 is accurate to **metres**; we need ~10⁻⁴ relative, four orders of magnitude
  looser (doc 10 §1)
- `rsspice` is **bit-identical** to ANISE — `0.0e0` km difference
- Evaluation is **Chebyshev polynomial** — already nanoseconds

A learned surrogate would be **less accurate and slower than the polynomial it
replaces.** There is nothing to reduce and nothing to accelerate.

### What "next-gen" actually means at this layer

Algorithms, not models. The genuine wins are already identified in doc 06:

| Win | Why it beats ML |
|---|---|
| **Angle-domain root finding** | O(events) instead of O(sample rate). Analytic + Newton, converges in 2–3 iterations, machine precision. No model could match it. |
| **Harmonic ephemeris precompute** | Frequencies are near-constant; precompute (freq, phase, amplitude) once → O(1) queries. |
| **Kernel subsetting** | A real unsolved problem RustSPICE lists as unbuilt. It is a data-engineering problem, not a learning problem. |

The one arguable ML case is **learned compression for edge deployment** — a compact
representation for Star Seer in a browser instead of shipping DE440. Even there,
refitting Chebyshev coefficients over a narrow time window is the classical answer
and probably wins on both size and accuracy.

**Recommendation: Layer 1 stays at L0 (no model), permanently.** Effort spent on ML
here is effort not spent where a gap exists.

---

## Layer 2 — PlanetaryHarmonics: **the intermediate ML, and it is real**

This is the layer the redirect was reaching for, and it does have genuine gaps.

Its defining advantage: **Layer 2 ML can be validated against ground truth we can
compute ourselves.** A surrogate is checked against direct computation. A
decomposition is checked against reconstruction error. **No seismological
uncertainty enters.** That is a luxury Layer 3 never has, and it is why this layer
should be built and validated first.

### 2a. Surrogates — speed where computation is genuinely expensive

**Target:** the ocean tidal loading convolution (Farrell Green's functions × global
ocean model), and later spherical-harmonic synthesis over a global HEALPix mesh.

**Not** the tidal tensor — that is ~10 flops and already free (see
`crates/ph-core/src/tidal.rs`). A surrogate there would be slower than the truth.

**Method ladder:** per-location constituent admittances (a lookup, not a model) →
polynomial/Chebyshev fits → small MLP if those fail.

**Validation:** max and RMS error against direct computation, with a published error
bound shipped alongside every surrogate.

This is the honest answer to the original brief's "hyperspeed bulk processing."
Speed comes from **precompute plus surrogates for the few expensive pieces**, not
from learning the physics we already have in closed form.

### 2b. Representation learning — taming 12,935 waves

**Target:** HW95 has 12,935 waves, 1,483 of them planetary. That is far too many
features for any downstream model, and a multiple-testing catastrophe.

**Method ladder:** amplitude-prior selection (physics, not ML) → PCA/SVD →
autoencoder only if the linear methods leave structure behind.

**Honest expectation:** physics-guided selection by tidal potential amplitude will
probably beat PCA, because the basis is already physically structured. Test both;
expect the boring answer to win.

### 2c. Transfer function estimation — the real Layer 2 science ML ★

**This is the intermediate model the project actually needs.**

**Target:** given known forcing amplitudes (HW95/KSM03) and observed event times,
estimate the response `R(ω)` and phase `φ(ω)`.

This is an **inverse problem, not forecasting.** Its output is a physical
measurement — `Aσ₀`, `τ̇`, `t_n` — not an opaque score.

**Method ladder:**

| Level | Method | Notes |
|---|---|---|
| L1 | Parametric fit of the rate-and-state form | Few parameters, directly interpretable |
| L2 | Gaussian process over log-frequency | Non-parametric, gives uncertainty, no functional-form commitment |
| L3 | **Fourier Neural Operator** | FNOs learn operators *in Fourier space* — an unusually exact match to learning a forcing→response operator |

The FNO fit is worth flagging: we are literally learning a mapping between
functions in the frequency domain, which is the problem class FNOs were designed
for. Worth trying **after** L1 and L2 establish a baseline, never before.

### 2d. Symbolic regression — for the gap in the theory

**This one is unusually well-motivated.** Heimisson & Avouac state their analytical
model is **invalid near T ≈ t_a** — which is exactly the region containing the
spectral peak we care about (doc 07, third pass).

So we have a theory that fails precisely where we need it. That is the textbook
case for **symbolic regression** (PySR, SINDy): recover the functional form from
data, then compare it against both asymptotic regimes the theory *does* provide.

Success would be a closed-form response law valid through the transition — a real
contribution to the rate-and-state literature, independent of any forecasting
result.

---

## Layer 3 — Earthquake Finder: irreducible inference

Covered in doc 04. Point processes, frozen-ETAS residual, the β(x,t) field, CSEP
evaluation. The only layer where ML does work physics genuinely cannot do — and the
only one where validation is against a noisy, contested signal.

**Do not start here.** Not because it is hard, but because a failure here is
unattributable until Layers 1 and 2 are known-good.

---

## The evolution ladder (per layer)

A capability level applies **within** each layer:

| Level | Capability |
|---|---|
| **L0** | Direct computation. No model. |
| **L1** | Fitted parametric model — physics-derived form, few parameters |
| **L2** | Classical statistical learning — GLM, GP, PCA |
| **L3** | Deep models — NN surrogates, neural operators, neural point processes |
| **L4** | Generative / foundation-scale |

**The discipline: a layer may not advance a level until the level below is
validated and *beaten on held-out data*.**

Targets:

| Layer | Current | Target | Ceiling |
|---|---|---|---|
| Ephemeris | L0 | **L0** | L0 — stays here |
| Harmonics | L0 | **L2 → L3** | L3 |
| Earthquake | — | **L3** | L3 |

L4 is not on the roadmap for any layer. Nothing in this problem justifies it, and
claiming otherwise would be the same overreach the project is built to avoid.

---

## "Models developing each other" — which patterns are real

The proposal's phrasing. Four candidate patterns, and they are not equally sound:

**✅ Surrogate enables search.** A fast Layer-2 surrogate makes millions of Layer-3
experiments feasible. Legitimate, and the main practical payoff of 2a.

**➖ Layer 2 outputs become Layer 3 features.** True, but that is just the pipeline.
Not models developing each other.

**⚠️ Layer 3 informs Layer 2 feature selection.** Tempting and **dangerous**. If
Layer 3 tells Layer 2 which features to compute precisely, and Layer 3 is then
evaluated on those features, **that is leakage** and it invalidates the result. Any
such loop must live strictly inside a cross-validation fold, or run on a discard set
never used for evaluation. Worth a written rule, because it is the kind of shortcut
that looks like efficiency.

**✅✅ Curriculum across the validation ladder — the strong one.**

Doc 05's difficulty ladder is **moonquakes → tremor → earthquakes**, ordered by
signal strength. That is exactly a **curriculum learning schedule**: pretrain the
transfer-function estimator where the effect is dominant and unambiguous, then
fine-tune where it is strong, then transfer to where it is marginal.

**The scientific de-risking sequence and the ML training sequence turn out to be the
same sequence.** The ladder we designed to make failures cheap also happens to be
the right order to train in. That is the version of "models developing each other"
worth building.

---

## Summary

1. **Layer 1 gets no ML.** It gets better algorithms. Say so plainly and move on.
2. **Layer 2 is the intermediate ML tier**, and its defining advantage is
   **computable ground truth** — build and validate it before touching Layer 3.
3. **Layer 2's flagship is the transfer function estimator** (2c) — an inverse
   problem whose output is a physical measurement, with symbolic regression (2d)
   attacking the region where the published theory admits it fails.
4. **Every layer climbs L0 → L3 on its own ladder**, and may not advance until the
   level below is beaten on held-out data.
5. **The validation ladder is the training curriculum.**
6. **Cross-layer feedback is leakage** unless fenced inside CV folds. Write it down.

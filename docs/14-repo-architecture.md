# 14 — Repository Architecture

Where earthquake forecasting lives, and why it is not here.

---

## Verdict

**Yes, a separate repository — at the start of Phase 3, not before.**

And name it **`EarthquakeForecastModule`**, not `EarthquakePredictionModule`. See
§4.

---

## 1. The chain forks at PlanetaryHarmonics

The original dependency chain was drawn as a straight line:

```text
RustSPICE -> PlanetaryHarmonicsModule -> AstrologyCore -> applications
```

**That is wrong for earthquake work.** Earthquake forecasting needs tidal tensors,
Coulomb stress, and harmonic decomposition. It needs **nothing** from the zodiac,
houses, decans, divisional systems, or interpretive mappings.

The real shape:

```text
RustSPICE
  └─> PlanetaryHarmonicsModule
        ├─> AstrologyCore ──> Cosmic Cypher, Star Seer, Resonant Finder
        └─> EarthquakeForecastModule          ← no AstrologyCore dependency
```

This is not a stylistic preference. It follows directly from the claim discipline
in [00-framing.md](00-framing.md): **if the forecasting module depended on the
interpretive layer, it would inherit that layer's baggage** — and every reviewer
would be right to discount it. The dependency graph is part of the argument.

---

## 2. Why separate rather than a subdirectory here

**Different toolchains.** PlanetaryHarmonics is Rust, targeting WASM. Forecasting
is Python — PyTorch, `pyCSEP`, ETAS tooling, RECAST. Different CI, different
packaging, different dependency resolution. Forcing them into one repo means one of
them is always the awkward guest.

**Different artefacts.** Forecasting carries earthquake catalogues, model
checkpoints, and training runs. Those do not belong in a library repo. We already
hit a rejected push from a 199 MB build artefact; model checkpoints are worse and
they are *supposed* to exist.

**Different audiences.** PlanetaryHarmonics is a reusable library for four
downstream projects. Forecasting is one application among them, sibling to Star
Seer and Resonant Finder — both of which are already separate.

**Reputational firewall, in both directions.** PlanetaryHarmonics should read as
sober, citable geophysics. Forecasting makes claims that are contested by
construction. If forecasting overclaims, it must not contaminate the library; if
forecasting produces a null result, the library still stands on its own merits
(and per [12-build-plan.md](12-build-plan.md), a rigorous upper bound is itself
publishable).

**The repo boundary enforces the leakage rule.** [13-ml-stack.md](13-ml-stack.md)
warns that Layer 3 informing Layer 2 feature selection is leakage that invalidates
results. A physical repository boundary makes that discipline structural rather
than a matter of remembering: forecasting cannot quietly reach into feature
generation and tune it.

---

## 3. When to split — at Phase 3, not now

The build plan's phases determine the split point:

| Phase | Content | Home |
|---|---|---|
| **1** Deep moonquakes | Validation of this library | **PlanetaryHarmonics** ✓ done |
| **2** Tectonic tremor | Validation of this library | **PlanetaryHarmonics** ✓ done |
| **3** Earthquakes | *see correction below* | **split** |

The line is **not** "is it Earth?" — it is **"are we forecasting?"**

### ⚠⚠ Second correction, 2026-08-22: the whole seismology programme moved

The correction below argued Phase 3 splits by "is it forecasting?", keeping
measurement upstream. **That test was too clever.** The sharper question is:

> **Who else needs this?**

PlanetaryHarmonics serves four downstream projects. Star Seer does not need a USGS
earthquake parser. Cosmic Cypher does not need Parkfield LFE families. Resonant
Finder does not need a moonquake nest catalogue. **A module only one consumer needs
is application code, not library code** — regardless of whether it computes a
forecast or a measurement.

By that test the entire seismology programme belongs downstream, including the
moonquake and tremor validation. Those validated the library, but
validating-a-library-for-a-purpose is the application's job.

**Moved to EarthquakeForecastModule:** `apollo`, `parkfield`, `cascadia`, `comcat`;
all 14 research examples; the catalogue fetch scripts; and docs 01, 04, 05, 07, 08,
09, 12, 16 — including the research log.

**Stayed:** `tidal`, `fault`, `love`, `doodson`, `stats`, `phase`, `demod`, `field`,
`ephemeris`, `catalog` (an interface type, not a parser); the CLI and Python
bindings; `fetch-kernels.sh`, which every consumer needs; and the library docs.

Two benefits beyond tidiness. PHM becomes **reviewable as a library** — evaluating
the tidal tensor no longer means wading through an earthquake research log. And EFM
becomes a **coherent scientific narrative** rather than a stub with a handoff
document attached.

The correction below is retained as a record of the reasoning it replaced.

### ⚠ Correction, 2026-08-22: Phase 3 is not homogeneous

The line above is right; Phase 3 was assigned to the wrong side of it. Phase 3
contains both measurement and forecasting, and only the second half is forecasting:

| Work | Forecasting? | Home |
|---|---|---|
| ComCat ingestion | no — a catalogue module, like `parkfield`/`cascadia` | **PlanetaryHarmonics** |
| Magnitude-of-completeness analysis | no — catalogue quality tooling | **PlanetaryHarmonics** |
| GCMT focal mechanisms | no — physics input to `fault` | **PlanetaryHarmonics** |
| **P3.4 — `R(ω)` for ordinary crust** | **no — it is a measurement** | **PlanetaryHarmonics** |
| ETAS baseline, residual model, β(x,t), CSEP | yes | **EarthquakeForecastModule** |

**P3.4 is the third rung of the validation ladder, not the start of forecasting.**
Doc 05's ladder runs moonquakes → tremor → earthquakes, and §3 of this document
already noted that the ladder spans the split. Measuring `R(ω)` for ordinary crust
uses the *identical* code path as C4 and C5 — `doodson`, `fault`, `love`, the
block-shift null — with no ETAS anywhere in it.

Three consequences:

1. The three-site comparison (moonquakes, two tremor sites, earthquakes) is **one
   coherent result** and should not be fragmented across repositories.
2. Forecasting consumes `ph-core` through the submodule exactly as `parkfield` and
   `cascadia` already do, so nothing is duplicated by keeping ingestion upstream.
3. **The handoff improves materially.** EarthquakeForecastModule then starts with
   *"here is the transfer function for ordinary crust"* rather than *"go find
   out"* — a foundation instead of a to-do list. Its `HANDOFF.md` §11 currently
   lists the band prediction as the central open question; answering it upstream
   first turns that section into a result.

Phases 1 and 2 are the library proving it measures what it claims to measure,
against known answers. That is library work and belongs here. Phase 3 is where
ETAS, the β(x,t) field, CSEP evaluation, and catalogue ingestion begin — none of
which any other downstream project needs.

**Splitting now would be premature.** Cross-repo changes are expensive during rapid
iteration, and there is nothing yet to put in the second repo.

---

## 4. On the name

Per [00-framing.md](00-framing.md)'s terminology table: **avoid "prediction."**

In seismology, *earthquake prediction* has a specific and largely discredited
meaning — deterministic statements of time, place and magnitude. *Forecasting*
means probabilistic rate estimation, which is what CSEP evaluates and what this
project actually does.

Naming the repository `EarthquakePredictionModule` would signal crankery to exactly
the audience whose methods we are borrowing. The cost is zero and the signal is
strong.

**Recommended: `EarthquakeForecastModule`.**

---

## 5. The interface between them

If the repos are separate, the contract matters. Two consequences:

**PlanetaryHarmonics emits, forecasting consumes.** Outputs are feature vectors and
the measured transfer function, each carrying the provenance metadata required by
[06-engine-architecture.md](06-engine-architecture.md) §6 — reference frame, epoch
system, kernel versions, feature tier (A/B/C), units. A feature whose frame is
ambiguous is worthless in a statistical test, and across a repo boundary that
ambiguity is much easier to introduce.

**PlanetaryHarmonics needs Python bindings.** This is the practical design
implication and it is worth planning for now rather than retrofitting. Rust core,
`pyo3`/`maturin` bindings, so the forecasting stack can consume features natively
without a serialisation layer or a subprocess.

That binding is also what Layer 2's own ML work (doc 13 §2c — transfer function
estimation) will need, since that will be fitted in Python too. So it is not
speculative work for a future repo; it is needed here regardless.

---

## 6. What this changes now

Nothing blocking, but two things to carry forward:

1. **Plan for `pyo3`/`maturin` bindings** on `ph-core`. Needed for Layer 2 ML
   before the split ever happens.
2. **Keep Phase 3 concerns out of `ph-core`.** No ETAS, no catalogue-specific
   forecasting logic, no CSEP. Tidal physics, harmonics, and statistics only. The
   `catalog` module stays deliberately minimal for this reason.

Doc 05's validation ladder is unaffected — it spans the split, which is the point:
the same instrument, validated on moonquakes and tremor here, is what forecasting
picks up.

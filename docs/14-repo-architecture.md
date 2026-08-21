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
| **2** Tectonic tremor | Validation of this library | **PlanetaryHarmonics** |
| **3** Earthquakes | Forecasting | **EarthquakeForecastModule** |

The line is **not** "is it Earth?" — it is **"are we forecasting?"**

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

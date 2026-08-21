# PlanetaryHarmonicsModule — Research Documentation

Working notes and literature base for the PlanetaryHarmonics computational framework.

**Scope:** deriving astronomical, tidal, gravimetric, and multi-basis harmonic features
for statistical analysis of Earth-system events. This module is the *scientific
computation layer*. It produces numbers, not interpretations. Symbolic and
interpretive mappings live upstream in `AstrologyCore`.

## Index

| Doc | Contents |
|---|---|
| [00-framing.md](00-framing.md) | Project scope, terminology, claim discipline |
| [01-literature.md](01-literature.md) | Tidal triggering literature, study by study |
| [02-angular-encoding.md](02-angular-encoding.md) | Fourier vs. base-harmonic encodings; Doodson; d'Alembert |
| [03-tidal-tensor.md](03-tidal-tensor.md) | Gravimetric stress tensor, concentration nodes, physics ceiling |
| [04-ml-architecture.md](04-ml-architecture.md) | Point-process modelling, ETAS residual design, validation |
| [05-research-frontier.md](05-research-frontier.md) | Under-explored options and validation datasets |
| [06-engine-architecture.md](06-engine-architecture.md) | Rust/WASM design, RustSPICE integration |
| [07-research-log.md](07-research-log.md) | Dated log of literature passes and decisions |
| [08-hypotheses.md](08-hypotheses.md) | Generated research hypotheses, prioritised |
| [09-deep-dive-agenda.md](09-deep-dive-agenda.md) | Equations and methods to derive next |
| [10-rustspice-requirements.md](10-rustspice-requirements.md) | Data and API spec for the RustSPICE layer |
| [11-osint-access.md](11-osint-access.md) | Free literature access routes; no institutional credentials |
| [12-build-plan.md](12-build-plan.md) | Decision to build; phased scope |
| [13-ml-stack.md](13-ml-stack.md) | Layered ML stack, capability ladders, what earns a model |
| [14-repo-architecture.md](14-repo-architecture.md) | Where forecasting lives; the chain forks here |
| [15-earthquake-forecast-handoff.md](15-earthquake-forecast-handoff.md) | **Copy to the new repo as `HANDOFF.md`** — self-contained |
| [bibliography.md](bibliography.md) | Full citations, tiered |

## Dependency chain

```
RustSPICE  (modules/RustSPICE)   ephemeris + time + coordinate frames
  └─> PlanetaryHarmonicsModule   harmonic / tidal / gravimetric features
        └─> AstrologyCore        symbolic and interpretive layer
              └─> Cosmic Cypher, Earthquake Finder, Resonant Finder, Star Seer
```

## Status

**Building.** `crates/ph-core` scaffolded; Phase 1 is deep moonquake validation
(see [12-build-plan.md](12-build-plan.md)). Literature base and
methodology design first, code second.

**Current strategic goal** (see [08-hypotheses.md](08-hypotheses.md)): convert the
multi-basis harmonic framework from a correlation search into a physically
grounded inverse problem with analytic significance testing. Critical path is
items 1, 6 and 2 of [09-deep-dive-agenda.md](09-deep-dive-agenda.md).

**Validation ladder** (see [05-research-frontier.md](05-research-frontier.md)):
deep moonquakes (tidal forcing dominant) → tectonic tremor (effect strong) →
earthquakes (effect small). Known-answer tests before the hard problem.

# PlanetaryHarmonicsModule — Research Documentation

Working notes and literature base for the PlanetaryHarmonics computational framework.

**Scope:** deriving astronomical, tidal, gravimetric, and multi-basis harmonic features
for statistical analysis of Earth-system events. This module is the *scientific
computation layer*. It produces numbers, not interpretations. Symbolic and
interpretive mappings live in `AstrologyCore`.

**It holds only what more than one downstream project needs.** Catalogue parsers,
research examples and domain analyses live with the application that needs them.
The seismology programme — moonquake validation, two-site tremor replication, and
the earthquake band-prediction test — moved to
[EarthquakeForecastModule](https://github.com/Jthora/EarthquakeForecastModule),
along with its research log. See [14-repo-architecture.md](14-repo-architecture.md).

## Index

| Doc | Contents |
|---|---|
| [00-framing.md](00-framing.md) | Project scope, terminology, claim discipline |
| [02-angular-encoding.md](02-angular-encoding.md) | Fourier vs. base-harmonic encodings; Doodson; d'Alembert |
| [03-tidal-tensor.md](03-tidal-tensor.md) | Gravimetric stress tensor, concentration nodes, physics ceiling |
| [06-engine-architecture.md](06-engine-architecture.md) | Rust/WASM design, RustSPICE integration |
| [10-rustspice-requirements.md](10-rustspice-requirements.md) | Data and API spec for the RustSPICE layer |
| [11-osint-access.md](11-osint-access.md) | Free literature access routes; no institutional credentials |
| [13-ml-stack.md](13-ml-stack.md) | Layered ML stack, capability ladders, what earns a model |
| [14-repo-architecture.md](14-repo-architecture.md) | Where forecasting lives; the chain forks here |
| [bibliography.md](bibliography.md) | Full citations, tiered |

## Dependency chain

```
RustSPICE  (modules/RustSPICE)   ephemeris + time + coordinate frames
  └─> PlanetaryHarmonicsModule   harmonic / tidal / gravimetric features
        └─> AstrologyCore        symbolic and interpretive layer
              └─> Cosmic Cypher, Earthquake Finder, Resonant Finder, Star Seer
```

## Status

Library. `ph-core` (79 tests), `ph-features` CLI, Python bindings via `ph-py`,
and a WebAssembly/TypeScript surface via `ph-wasm`. Validated by the seismology programme downstream: deep moonquake
periodicities recovered to better than 0.21%, and tidal response replicated at two
independent tremor sites. Literature base and
methodology design first, code second.

**Current strategic goal** (see [08-hypotheses.md](08-hypotheses.md)): convert the
multi-basis harmonic framework from a correlation search into a physically
grounded inverse problem with analytic significance testing. Critical path is
items 1, 6 and 2 of [09-deep-dive-agenda.md](09-deep-dive-agenda.md).

**Validation ladder** (see [05-research-frontier.md](05-research-frontier.md)):
deep moonquakes (tidal forcing dominant) → tectonic tremor (effect strong) →
earthquakes (effect small). Known-answer tests before the hard problem.

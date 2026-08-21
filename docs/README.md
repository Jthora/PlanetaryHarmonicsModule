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

## Dependency chain

```
RustSPICE  (modules/RustSPICE)   ephemeris + time + coordinate frames
  └─> PlanetaryHarmonicsModule   harmonic / tidal / gravimetric features
        └─> AstrologyCore        symbolic and interpretive layer
              └─> Cosmic Cypher, Earthquake Finder, Resonant Finder, Star Seer
```

## Status

Early research phase. No implementation committed yet. Literature base and
methodology design first, code second.

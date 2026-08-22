# 16 — Working Plan

Established 2026-08-21, after Phase 1's periodicity validation passed
([07-research-log.md](07-research-log.md)).

---

## Investigation findings that shaped this plan

**1. All Phase 1 kernels are available and small.** Confirmed at NAIF:

| Kernel | Purpose |
|---|---|
| `naif0012.tls` | Leap seconds |
| `de440s.bsp` | Ephemeris — **32.7 MB**, covers 1849–2150 (full `de440` is 114 MB and unnecessary) |
| `pck00011.tpc` | Body constants |
| `gm_de440.tpc` | GM values — the tensor is GM/d³ |
| `moon_pa_de440_200625.bpc` | Lunar orientation with libration |
| `moon_de440_250416.tf` | Lunar frame definitions |

Fetch with `./scripts/fetch-kernels.sh`. **Task 1 is unblocked.**

**2. Deep moonquake focal mechanisms are poorly constrained — and that is useful.**

Weber, Bills & Johnson (2009), *JGR* 114, E05001, and Weber & Knapmeyer
(NTRS 20120013649, **free**): Apollo's limited source/receiver geometry
*prohibits* determining fault parameters from first-motion polarities. Their
S/P-amplitude inversion reduces the fault-plane parameter space by only about half
on average; the best-constrained cluster eliminates 72% of the focal sphere.

Their method: **grid search over strike/dip/slip, with the failure criterion that a
linear combination of shear and normal stress best approximates a constant** at
event times. Data: 106 events across 25 clusters, 37 of them from A001.

**This reorders the plan.** Three consequences:

- **Coulomb projection must treat fault orientation as a *parameter*, not an
  input.** On the Moon we search for it; on Earth we fix it from GCMT. **One code
  path, two uses** — and building it against the harder case first means the Earth
  case is a simplification rather than a rewrite.
- **Weber 2009 is a published benchmark**, giving Phase 1 a *second* known-answer
  validation — this time of the stress projection, not just the timing analysis.
- They find **stress *rates*** matter for some clusters, independently supporting
  the `dΔCFS/dt` emphasis in [03-tidal-tensor.md](03-tidal-tensor.md).

So **task 4 (Coulomb projection) moves up**, and merges with task 3.

---

## Track A — Finish Phase 1  ◀ start here

| # | Task | Depends on |
|---|---|---|
| **A1** | Fetch kernels; wire `ph-core` to load them | — |
| **A2** | Compute lunar tidal tensor from real Earth + Sun geometry, in `MOON_PA` | A1 |
| **A3** | Tidal phase from the full quasi-periodic forcing (not trial-period folding) | A2 |
| ~~A4~~ | ~~Time-shift null on the pooled catalogue~~ — **does not work; see below** | — |
| **A4′** | **Per-nest** phase test: is the preferred phase consistent within a nest and structured across nests? | A3 |
| **A5** | Coulomb projection with fault orientation as a searchable parameter | A2 |
| **A6** | Validate A5 against Weber 2009 per-cluster fits; use `nakamura_2005_dm_locations.csv` (already downloaded, unused) | A5 |

**A3 is the linchpin.** The time-shift null is degenerate against a single exact
frequency, so A4 cannot be done correctly without it.

> **A4 failed as specified, 2026-08-21.** The time-shift null is degenerate
> whenever the **catalogue** shares the forcing's period — not merely when the
> forcing is a single sinusoid. Deep moonquakes are locked near the anomalistic
> month and so is the forcing, so a shift rotates the phase cluster without
> diluting it. Measured: analytic Schuster p = 1.07e-89, empirical p = 0.699.
> See [07-research-log.md](07-research-log.md).
>
> **Pooled-catalogue tests cannot produce a falsifiable Phase 1 claim.** A5/A6
> (per-nest Coulomb projection against Weber's constraints) are now the critical
> path, not a bonus validation.

**Exit criterion:** the Coulomb projection reproduces Weber's per-cluster fault
constraints, and preferred phase is consistent within nests and structured across
them.

## Track B — Prepare the split

| # | Task | Depends on |
|---|---|---|
| **B1** | `pyo3`/`maturin` bindings on `ph-core` | — |
| **B2** | CLI emitting feature vectors as CSV/Parquet with full provenance | — |
| **B3** | Create `EarthquakeForecastModule`; add `HANDOFF.md` ([15](15-earthquake-forecast-handoff.md)); add `ph-core` as a submodule | A5 |

B1 is needed by Layer 2's own ML regardless of any split
([13-ml-stack.md](13-ml-stack.md) §2c), so it is not speculative.

## Track C — Phase 2, tremor

| # | Task | Depends on |
|---|---|---|
| **C1** | Ingest a Cascadia or Parkfield LFE catalogue | A3 |
| **C2** | Full exponential response `R₀exp(S_T/Aσ₀)` with M > 1 — tremor is non-linear | C1 |
| **C3** | Measure the response spectrum where the effect is strong | C2 |

## Track D — Test the band prediction

| # | Task | Depends on |
|---|---|---|
| **D1** | M2-vs-S2 and alias-analysis validity gates | A3 |
| **D2** | Free HW95/KSM03 catalogue data (Kudryavtsev & Cionco 2025, arXiv:2508.18111, ships in HW95 format) | — |
| **D3** | Broadband response spectrum — tides, hydrological, pole tide, LOD in one fit | C3, D2 |

**D3 is the project's central question.** Everything before it is instrument-building.

## Track E — Layer 2 ML

| # | Task | Depends on |
|---|---|---|
| **E1** | Transfer function estimator — parametric fit first, GP second | C3, B1 |
| **E2** | Symbolic regression on the `T ≈ t_a` gap where the published theory fails | E1 |
| **E3** | Ocean loading surrogate | **only if** the band prediction fails |

E3 is explicitly conditional. If the responsive band is years-to-decades, ocean
loading barely matters and this is wasted effort.

## Track F — Research, parallel

| # | Task |
|---|---|
| **F1** | PTA red-noise machinery — port `enterprise`-style methods |
| **F2** | IERS Conventions Ch. 7 — correct Love numbers, frequency-dependent corrections |
| **F3** | Obtain Weber 2009 full text (ResearchGate / AGU PDF) for exact cluster parameters |

---

## Critical path

```
A1 → A2 → A3 → A4        Phase 1 complete
          A2 → A5 → A6   Coulomb validated against Weber
                    B1, B2 → B3        split
          A3 → C1 → C2 → C3 → D3       the central question
```

**A3 gates almost everything.** It unblocks the null (A4), tremor (C1), and the
validity gates (D1).

## Decision points

1. **After A4** — if the known periods do *not* survive the null, stop and fix the
   pipeline. Nothing downstream is meaningful.
2. **After A6** — if Coulomb projection cannot reproduce Weber's constraints, the
   stress projection is wrong and Earth work must wait.
3. **After D3** — if the band prediction fails and short periods dominate, ocean
   loading (E3) returns to the critical path and the schedule changes materially.
4. **Before any Phase 3 work** — pre-commit to what gets published if the answer is
   null: a rigorous upper bound plus the moonquake and tremor validations. Deciding
   now removes the incentive to keep scanning until something crosses p < 0.05.

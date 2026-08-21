# 02 — Angular Feature Encoding

How to turn relative geometry between bodies into ML features. This is the
technical core of the module.

---

## The question

Given a relative angle θ (e.g. difference in ecliptic longitude between two
bodies), what is the right input representation for a predictive model?

Three candidate encodings:

**E1 — Fourier / circular embedding**
```
φ_F(θ) = [cos θ, sin θ, cos 2θ, sin 2θ, …, cos Nθ, sin Nθ]   ∈ ℝ^(2N)
```

**E2 — Base-b sector one-hot**
```
φ_b(θ) = onehot( floor(b·θ / 2π) )   ∈ ℝ^b     stacked over b ∈ B
```

**E3 — Base-b harmonic scalar / "orb kernel"**
```
φ_b(θ) = cos(bθ)         or       s_b(θ) = exp( κ_b (cos bθ − 1) )
```
(the latter a von Mises bump peaking at exact b-fold divisions)

---

## Why E1 is the right substrate

### E3 is a strict subset of E1 — and loses the phase

`cos(bθ)` is exactly one component of E1 at order *n = b*. So "injecting a series
of base-number angular harmonics as separate inputs" is the **cosine-only,
phase-locked-to-zero subset** of the Fourier encoding.

That missing sine term is not cosmetic. Without it the model can only express a
response symmetric about θ = 0, and **cannot learn a phase offset** — it cannot
represent "rate peaks 40° after conjunction."

This matters enormously, because phase offset is *precisely where the tidal
triggering signal lives*:
- Métivier et al. find events cluster at the **uplift** phase, not at zero
- Scholz et al. find mid-ocean ridge events peak at **low tide** — the opposite
  of the naive expectation

Dropping the sine term assumes the answer. **Keep the (cos, sin) pair.**

### E2 costs more and buys less

A base-*b* sector indicator, Fourier-expanded, is a sum of harmonics at orders
b, 2b, 3b, … with *fixed* amplitude ratios (the Dirichlet kernel). So E2 also lies
inside the Fourier span, but with frozen weights, and it is discontinuous at bin
edges. It spends *b* dimensions to buy less expressiveness than 2 dimensions of
Fourier.

Worse for gradient-based learning: adjacent sectors share no parameters, so the
model must *learn* that sector 1 and sector 2 are neighbours. Circular structure
that E1 gets for free must be discovered from data.

### The decisive advantage: the spectrum is a readout

Fit a linear or GLM layer on φ_F. The coefficient pair (aₙ, bₙ) at order *n*
yields:

```
amplitude  Aₙ = √(aₙ² + bₙ²)
phase      ψₙ = atan2(bₙ, aₙ)
```

**Aₙ is a direct, testable estimate of how much base-n structure exists in the
data.** You get "which base matters" as an *output* instead of pre-committing to a
base list as an *input*.

For a project whose entire premise is *"not just base 12"*, this is the whole
ballgame. If order 12 lights up, that is empirical support for the zodiacal
division. If order 7 lights up and 12 does not, that is a discovery. One-hot
binning cannot give you this.

### It is the Schuster test, generalised

Aₙ and ψₙ at order *n* are exactly the **generalised Schuster test** statistic at
harmonic *n*. The classical Schuster test used throughout the tidal triggering
literature (Tanaka et al. 2002) is the **n = 1 case**.

So the ML feature basis reduces, at first order, to the field's own canonical
statistic. That is a direct bridge to reviewers who already know the method, and
it means every learned coefficient has a published significance test attached.

### Classical aspects are the low-order subspace

| Aspect | Angle | Harmonic order |
|---|---|---|
| Conjunction / opposition | 0° / 180° | n = 1, 2 |
| Trine | 120° | n = 3 |
| Square | 90° | n = 4 |
| Sextile | 60° | n = 6 |
| Zodiacal division | 30° | n = 12 |

A model carrying n = 1…64 **strictly contains the classical system as a subspace**,
plus everything it never considered. Nothing is discarded; the tradition becomes a
testable special case.

---

## When E2/E3 would actually win

Be fair to the alternative. E3 is preferable when:

1. **The true response is genuinely sharp.** If a real resonance has a very narrow
   orb, Fourier needs many orders and may ring (Gibbs). Roughly N ≈ 3√κ orders to
   represent a von Mises bump of concentration κ. Astrology's 6–8° orbs are broad,
   so moderate N suffices — but a sharp effect would favour E3.
2. **Strong prior on specific bases.** Fewer parameters, better sample efficiency.
3. **Interpretability to practitioners.** "Trine strength" is more legible than a
   coefficient at order 3.

### Recommended resolution

Use **E1 as the substrate**. Add E3 orb-kernels only as an **explicit,
pre-registered hypothesis test** against the E1 baseline. If E3 beats E1 on
held-out likelihood, that is positive evidence for genuine sharp resonance —
itself a real finding worth reporting.

### Regularisation

Apply **group lasso over harmonic orders**, grouping each (cos nθ, sin nθ) pair.
This gives sparse harmonic selection while preserving phase freedom within each
retained order. The set of surviving orders *is* the answer to "which bases
matter."

---

## Multi-body terms: commensurabilities

For more than two bodies the relevant object is the set of integer combinations

```
Φ_k(t) = Σᵢ kᵢ · θᵢ(t)        kᵢ ∈ ℤ,  |kᵢ| ≤ K
```

with features `[cos Φ_k, sin Φ_k]`.

### The d'Alembert constraint — a structural filter

Impose:

```
Σᵢ kᵢ = 0
```

This is the **d'Alembert rule** from the classical disturbing-function expansion in
celestial mechanics. It enforces **rotational invariance**: no physical quantity
can depend on where you arbitrarily placed the origin of longitude.

Two consequences, and they are the same consequence:

1. It is required for physical correctness.
2. It **automatically annihilates every feature that depends on absolute zodiacal
   position** while retaining every feature based on relative geometry.

The constraint that makes the physics right is the constraint that removes the
indefensible features. It also collapses the combinatorial feature space
substantially, which helps the multiple-testing budget.

---

## Doodson numbers: the rigorous ancestor

The Darwin–Doodson harmonic development of the tide-generating potential
decomposes Sun–Moon–Earth geometry into a frequency basis. Each constituent is
indexed by **six integers** — the Doodson number — giving multiples of six
fundamental astronomical arguments:

| Symbol | Argument | Period |
|---|---|---|
| τ | lunar time | 24.84 h |
| s | Moon's mean longitude | 27.32 d (tropical month) |
| h | Sun's mean longitude | 365.24 d |
| p | longitude of lunar perigee | 8.85 yr |
| N′ | negative longitude of lunar ascending node | 18.61 yr |
| pₛ | longitude of solar perigee | 20940 yr |

This is a base-independent, multi-argument angular encoding system that already
exists, is validated by a century of tide prediction, and has exactly the shape
of the thing we are building.

### The long-period constituents are the "planetary" cycles

| Constituent | Period | Physical meaning |
|---|---|---|
| Mf | 13.66 d | lunar fortnightly (declination cycle) |
| Mm | 27.55 d | anomalistic month (perigee cycle) |
| Ssa | 182.62 d | solar semiannual |
| Sa | 365.26 d | solar annual |
| Node | 18.61 yr | lunar nodal regression |

Per Beeler & Lockner's nucleation argument (doc 01), **these long-period terms are
where correlation is most likely to be detectable** — the forcing period exceeds
the nucleation time. This is the mechanism-grounded reason to care about slow
orbital cycles, and it is exactly the regime classical astrology emphasises.

### ⚠ The novelty claim here was wrong — corrected fourth pass

This section originally proposed "a generalised Doodson expansion including
planetary arguments" as our novel contribution. **It was published in 1995.**

The **HW95 catalogue** (Hartmann & Wenzel 1995, *GRL* 22, 3553) contains **12,935
tidal waves, of which 1,483 are due to direct planetary effects** — Venus, Jupiter,
Mars, Mercury, Saturn. Scope: Moon to degree 6, Sun to degree 3, planets to degree
2. Based on DE200, covering 1850–2150, accurate to better than 1 nGal.

Its modern successor is **KSM03** (Kudryavtsev 2004), updated by **Kudryavtsev &
Cionco (2025)** (arXiv:2508.18111, free) with DE441 and **38 new/updated terms of
period > ~18 years**, including several new waves near the 18.61-year nodal period.

**This is a gift.** We do not derive the expansion — we *consume* it. An
authoritative, geodesy-grade basis of 1,483 planetary tidal waves with correct
amplitudes is exactly the multi-body harmonic feature basis this project wants, and
using it is a far stronger position than rolling our own.

### The corrected novelty claim

Not the expansion. Rather:

1. Using HW95/KSM03 planetary waves as an **ML feature basis tested against
   seismicity** — not previously done
2. **Broadband response spectroscopy** across six decades of forcing frequency
   ([08-hypotheses.md](08-hypotheses.md) §12)
3. A **global β field** — published work is single-fault
4. The **moonquake → tremor → earthquake** validation ladder

The d'Alembert constraint and Fourier-encoding arguments above stand unchanged —
they remain the right way to *use* the basis, and to extend it to non-tidal angular
features the catalogues do not cover.

---

## Implementation notes

- **Never feed raw degrees to a network.** Discontinuity at 0°/360° destroys
  gradients.
- Compute cos/sin via **angle-sum recurrence** from n = 1, not N separate
  trig calls — O(N) multiply-adds instead of O(N) transcendentals.
- All angles in the **ecliptic of date** unless a doc says otherwise; record the
  frame in output metadata.
- Feature count for the multi-body basis grows fast. With the d'Alembert
  constraint and |kᵢ| ≤ 3 over 10 bodies it is manageable; budget for group-lasso
  selection rather than emitting everything.

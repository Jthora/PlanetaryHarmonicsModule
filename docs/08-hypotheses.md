# 08 — Generated Hypotheses

Ideas produced by the first literature pass. These are candidate research
contributions, not established results. Each states the idea, why the literature
supports it, and what would have to be true for it to fail.

---

## 1. The harmonic spectrum as an inverse method for fault friction ★

**The idea.**

Beeler & Lockner showed that a fault's correlation with periodic stressing depends
on whether the forcing period exceeds the **characteristic nucleation time**. That
is not just a caveat — it is a statement that a fault acts as a **low-pass filter**
on periodic stress.

So a rate-and-state fault has a **transfer function**: amplitude response R(ω) and
phase lag φ(ω) to sinusoidal forcing. Which means the harmonic spectrum we are
already computing (doc 02) should show:

- Strong response at low harmonic orders / long periods
- A **knee** at the nucleation frequency
- Suppressed response above it

And crucially — **the location of the knee estimates the nucleation time, hence
in-situ frictional properties (aσ, D_c).**

**Why this is the most valuable idea here.**

It converts the multi-basis harmonic framework from *exploratory feature
engineering* into a **physically-grounded inverse problem**. We would no longer be
asking "does base N correlate?" but "what does the shape of the response spectrum
tell us about the fault?"

That is a real geophysical measurement, publishable independently of any
forecasting result, and it makes every harmonic order physically meaningful.

**Failure condition.** If the observed spectrum is flat, or the knee is not
resolvable given catalogue size, the inversion is unconstrained. Power analysis
(doc 04 §7) determines this before we invest.

**Deep dive needed:** derive R(ω) and φ(ω) analytically. See doc 09 §1.

---

## 2. Phase-versus-frequency as an anti-artifact test ★

**The idea.**

Rate-and-state predicts not just *whether* seismicity responds, but a **specific
phase lag that varies with frequency**. That is a two-dimensional prediction:
a curve φ(ω), not a single number.

A spurious correlation — multiple testing, catalogue bias, hydrological aliasing —
has **no reason whatsoever to reproduce the predicted lag-versus-frequency
relation.** It might accidentally produce a significant amplitude at one order. It
will not accidentally trace a theoretical phase curve across many orders.

**Why this matters.**

It is a far stronger test than amplitude significance, and it is essentially free
once the Fourier encoding is in place — the phase is already in the coefficients.

This directly answers the Holt & Newman failure mode: they tested amplitude-like
significance across thousands of intervals. Testing *structure across orders*
rather than *significance at any order* is qualitatively harder to fake.

**Recommend making this the project's primary evidentiary standard.**

**Failure condition.** Requires enough events per harmonic order to estimate phase
with useful uncertainty. Phase estimates are noisier than amplitude estimates.

---

## 3. Use the tidal potential catalogue as the prior over harmonic orders ★

**The idea.**

We have been treating "which base numbers matter" as an open search. But the
**Cartwright–Tayler–Edden** tidal potential catalogue already tells us the true
amplitude of the astronomical forcing at each frequency — several hundred
constituents with measured coefficients.

So: **weight the prior over harmonic orders by the actual tidal potential
amplitude at that frequency.**

**Why this is powerful.**

1. It collapses the search space from "all bases" to "physically forced
   frequencies," which massively reduces the multiple-testing burden — the single
   biggest statistical threat to the project.
2. It gives every feature a physical amplitude, so effect sizes become
   *predictions* rather than free parameters. A constituent with 10× the forcing
   amplitude should show a larger response; that ratio is testable.
3. It preserves the exploratory goal: if an **unforced** base lights up beyond its
   prior, that is precisely the anomaly worth investigating — and now it stands out
   against a principled background instead of being one of 28,782 tests.

This reframes "not just base 12" from a scan into a **structured comparison
against known forcing.**

**Failure condition.** If real response is dominated by non-tidal periodicities
(hydrological, catalogue artifacts), the prior misleads. Mitigate by fitting the
unforced orders too and comparing.

---

## 4. The sensitivity field β(x,t) as a tomographic inverse problem

**The idea.**

Beaucé's result is regional and retrospective. Generalised, **β(x,t)** — the
coupling between tidal stress and seismicity rate — is a continuous field to be
estimated from sparse point observations (earthquakes). That is structurally a
**tomography / spatial inverse problem**: Gaussian process or spline field with a
smoothing prior, estimated inside the point-process likelihood.

**The valuable part: resolution.**

Inverse problems come with **resolution kernels**. We can compute *where we have
enough events to constrain β and where we do not*, and publish that map alongside
the estimate.

An honest map of "here is where we can say something, and here is where we cannot"
is more scientifically valuable — and far more credible — than a global field with
uniform false confidence. It also directly addresses the objection that celestial
correlation studies over-claim spatial coverage.

**Failure condition.** β may vary faster than the event rate can track, making the
field under-determined everywhere. Resolution analysis reveals this honestly
rather than hiding it.

---

## 5. Model the detection function instead of truncating at Mc

**The idea.**

Standard practice truncates the catalogue at the magnitude of completeness,
discarding most events. Instead, put a **detection probability function
P_det(m, x, t)** directly into the point-process likelihood.

This is routine in astronomy (selection functions) and under-used in seismology.

**Why it matters here.**

Beeler & Lockner say we need 10⁵–10⁶ events. Mc truncation throws away the
majority of the catalogue — exactly the small, shallow events where Métivier found
the anomaly is **largest**. Modelling detection instead of truncating could
recover an order of magnitude more data in the regime where the signal is
strongest.

It also converts the Mc-drift problem from a confounder into a modelled quantity.
Detection improved over the catalogue span; unmodelled, that trend aliases into
long-period features and manufactures signal (doc 04 §6b). Modelled, it is handled.

**Failure condition.** A misspecified detection function introduces bias worse than
truncation. Validate against Mc-truncated results on subsets where both are viable.

---

## 6. Deep moonquakes yield lunar interior structure as a by-product

Beyond validating the pipeline (doc 05 §1), deep moonquake tidal analysis is an
established route to constraining **lunar Love numbers and deep interior
structure**.

So Phase 1 of the validation ladder is not merely a test harness — it can produce
an independent scientific result. That matters strategically: it establishes
methodological credibility in a domain where tidal forcing is uncontroversial,
*before* any terrestrial claim is made.

---

## 7. Pipeline self-test with synthetic null features

**The idea.**

Before trusting any positive result, construct a **synthetic feature with identical
spectral and autocorrelation properties to a real celestial feature but no physical
link to seismicity**, and run the complete pipeline on it.

The pipeline must return null. If it does not, the pipeline is broken — and we
learn that cheaply, before publishing.

Extend to a small suite: shuffled-phase features, features from a fictitious body
on a plausible orbit, features from a real body evaluated at wrong epochs.

This is cheap, decisive, and almost never done in this literature. It should be a
CI test, not a one-off.

---

## 8. Ocean tidal loading is the compute bottleneck, and probably the signal

Cochran et al.'s factor-of-3 result was driven by **ocean tidal loading**, not
solid Earth tide. Loading is a **spatial convolution** of ocean tide height with
Farrell Green's functions — far more expensive than the five-number degree-2 solid
tide, and unavoidable if we want the settings where the effect is largest.

Design implication: the architecture in doc 06 optimises the *cheap* part. The
expensive part is loading. It needs its own precompute strategy — likely
constituent-wise loading coefficients per location, computed once and stored, so
runtime cost reduces to harmonic synthesis.

Flagged as the most likely source of unpleasant surprise in the implementation
schedule.

---

## Priority for next pass

| # | Idea | Value | Cost | Priority |
|---|---|---|---|---|
| 1 | Rate-and-state spectrum inversion | Very high | Medium | **1** |
| 2 | Phase-vs-frequency test | Very high | Low | **1** |
| 3 | CTE amplitude prior | High | Low | **2** |
| 7 | Synthetic null self-test | High | Low | **2** |
| 5 | Detection function | High | Medium | 3 |
| 4 | β tomography + resolution | High | High | 3 |
| 8 | Ocean loading precompute | Required | High | 3 |
| 6 | Moonquake interior inversion | Medium | Medium | 4 |

Ideas 1, 2 and 3 together would convert this project from a correlation search
into a physically-grounded inverse problem with a built-in artifact test. That is
the strategic goal of the next research pass.

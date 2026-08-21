# 03 — Gravimetric Stress Tensor and Concentration Nodes

Addresses the question: *can we track a 3D geocentric "node of gravimetric
concentration" over time, and does its motion or trajectory relate to seismicity?*

Short answer: **yes for the Sun–Moon system, where it is rigorous, computable and
cheap. No for inter-planetary "gravitational centrality nodes," which are not a
stress mechanism.** The correct home for planetary influence is orbital
perturbation of the lunisolar tide. Details below.

---

## 1. The tidal tensor

For a body of mass M at geocentric distance d, the tide-generating potential at
surface position **r** (|**r**| = R ≪ d) expands as

```
V(r) ≈ (GM/d) Σ_{n≥2} (R/d)^n Pₙ(cos ψ)
```

with ψ the angle between **r** and the body direction. The dominant n = 2 term is a
**degree-2 spherical harmonic**. (The Moon's n = 3 term is small but non-zero;
everything else is negligible.)

The tidal tensor at Earth's centre, in a frame with **ẑ** toward the body:

```
T = (GM / d³) · diag(−1, −1, +2)
```

Trace-free: **extension along the body axis, compression transverse.** Five
independent components for a symmetric trace-free 3×3.

### Two properties that make this cheap

**Linear superposition.** `T_total = Σ_bodies T_body`. Summing the tidal tensors of
all bodies costs a handful of flops per body.

**Five numbers describe the entire global field.** At degree 2, the complete
worldwide tidal field at an instant is five coefficients. Evaluating it at any
surface point is a synthesis, not a recomputation. This is the architectural basis
for global bulk processing (see [06-engine-architecture.md](06-engine-architecture.md)).

---

## 2. The concentration node, done rigorously

The rigorous version of "node of gravimetric concentration":

> **The principal axis of the summed tidal tensor** — the eigenvector of the
> largest eigenvalue of `T_total`. Where it pierces the surface gives two
> antipodal **tidal concentration points**.

These are real, well-defined, and trace genuine trajectories:

- The **sublunar point** sweeps ±28.5° latitude, its range modulated by the
  **18.61-year lunar nodal cycle**
- The **subsolar point** sweeps ±23.4° annually
- At **syzygy** (new/full Moon) the solar and lunar principal axes align →
  spring tides, maximum combined amplitude
- At **quadrature** they are 90° apart and partially cancel → neap tides

**Testable features that follow directly:**

1. Angular distance from a location to the instantaneous concentration axis
2. Rate of change of that distance (the "trajectory" question, made precise)
3. Latitude of the sublunar point relative to its nodal-cycle extremes
4. Eigenvalue ratio of `T_total` — a scalar "tidal focus" measure
5. Angular velocity of the principal axis across the surface

Feature 2 is the direct answer to *"do earthquakes happen more often during a
specific motion or transition of the nodal position?"* — and it connects to
rate-and-state friction, where **dCFS/dt matters as much as CFS amplitude**.

---

## 3. What actually loads a fault

The tidal tensor is not the endpoint. The physically meaningful quantity is
**Coulomb failure stress** resolved on the fault plane:

```
ΔCFS = Δτ + μ′ · Δσₙ
```

- Δτ — shear stress resolved on the fault in the slip direction
- Δσₙ — normal stress (positive = unclamping)
- μ′ — effective friction (Cochran et al. found best fit at **μ = 0.4**)

This requires **fault geometry** (strike, dip, rake), available from **GCMT** focal
mechanisms or regional stress-field models.

Most naive studies use scalar tide height and find nothing. Projecting the full
tensor onto fault geometry is where the signal lives — and it is precisely what
"gravimetric stress tensor" should mean in this project.

**Do not omit ocean tidal loading.** In coastal and subduction settings it
frequently **dominates** the solid Earth tide, and it drove the strongest positive
result in the literature (Cochran et al. 2004). Model solid tide and ocean loading
separately; `SPOTL` (Agnew) is the standard tool.

Required outputs per event/cell:
- ΔCFS and **dΔCFS/dt** at event time
- Tidal phase angle at event time (for Schuster-type tests)
- Constituent decomposition (M2, S2, N2, K1, O1, Mf, Mm, Ssa, Sa)
- Solid-tide and ocean-loading contributions, separated

---

## 4. Planetary contributions — the honest ceiling

Tidal effect scales as **M / d³**. Computed:

| Body | M/d³ (SI) | Relative to Moon |
|---|---|---|
| Moon | 1.30 × 10⁻³ | 1 |
| Sun | 5.9 × 10⁻⁴ | ~0.46 |
| Venus (closest) | 8.9 × 10⁻⁸ | ~7 × 10⁻⁵ |
| Jupiter (closest) | 9.4 × 10⁻⁹ | ~7 × 10⁻⁶ |

Consistent with the standard figure that Venus at maximum is ~10⁴ times weaker
than Sun and Moon combined, and Jupiter roughly a further order below that.

**Consequence:** given that the *lunar* tide produces only a few percent rate
modulation globally, **direct planetary tidal triggering is far below detection.**
This must be stated plainly in any output of this project. Claiming otherwise is
the fastest way to be dismissed.

### The critical error to avoid

> **A uniform gravitational field produces no stress.**

Earth is in free fall. Only the *gradient* — the tidal tensor — produces
deformation. Therefore:

- The **solar-system barycentre** position relative to Earth is **not** a stress
  mechanism
- **Lagrange points** are force nulls in a rotating frame, not stress features, and
  say nothing about the gradient at Earth
- A "gravitational centrality node between planets passing through Earth" does
  **not** transfer stress

This is the exact point where this line of reasoning usually goes wrong. Being
explicit about it saves months and protects the credible parts of the work.

---

## 5. Where planetary influence legitimately enters

Two real channels, both computable, both Tier B.

### 5a. Orbital perturbation of the lunisolar tide

Planets perturb the orbital elements of the Moon and of Earth's orbit. Those
perturbations **modulate the large lunisolar tide.** The effect enters through the
dominant tide rather than through a planet's own negligible tide.

This is already encoded in the higher-order Doodson and Brown lunar-theory
arguments. It is established celestial mechanics, requires no new physics, and is
free from any decent ephemeris — which RustSPICE provides.

**This is the correct home for the planetary intuition.** Small, but principled,
non-zero, and undismissable.

### 5b. Principal-axis deflection — a genuinely novel feature

Adding planetary tidal tensors to `T_total` perturbs the **direction** of the
principal axis by a small angle δ, not just its magnitude.

Direction is a *normalised* quantity. Near **quadrature**, where the solar and
lunar tensors partially cancel, the denominator shrinks and a tiny perturbation
produces a disproportionately larger angular deflection:

```
δ ~ (T_planet / T_lunisolar) / sin(separation)
```

So **angular deflection of the combined principal axis may be a more sensitive
planetary probe than amplitude perturbation.** Cheap to compute — it falls out of
the eigendecomposition already being performed.

**Honest caveat, and it is a real one:** the amplification is largest exactly where
total tidal amplitude is smallest, i.e. where triggering is *least* likely. The
sensitivity gain works against the signal. Worth computing and testing — it is
nearly free — but do not expect it to rescue the planetary hypothesis.

### 5c. Pole tide and length-of-day

Genuinely real, physically sound, and under-explored:

- **Polar motion / Chandler wobble** (~433 d) changes the centrifugal potential,
  producing the **pole tide** — real stress at the ~10⁻⁹ level, comparable to some
  minor tidal constituents
- **LOD variations** connect to Bendick & Bilham (2017)

Both are cheap to include and are Tier A physics (with Tier C interpretation).
IERS provides the series.

---

## 6. Summary of what to compute

| Quantity | Tier | Notes |
|---|---|---|
| `T_total` (5 components) | A | Sum over all bodies, degree 2 (+ lunar degree 3) |
| Principal axis direction & eigenvalues | A | Concentration node |
| Sub-body point tracks | A | Trajectory features |
| ΔCFS, dΔCFS/dt on fault geometry | A | Requires GCMT mechanisms |
| Ocean tidal loading | A | Often dominant; use SPOTL |
| Constituent decomposition | A | Doodson-indexed |
| Pole tide, LOD | A | From IERS |
| Planetary orbital perturbation of lunisolar tide | B | Via ephemeris |
| Principal-axis angular deflection | B | Novel; cheap; expect small |
| Direct planetary tidal stress | C | ~10⁻⁵ of lunar; document as below detection |

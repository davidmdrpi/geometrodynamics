# The congruence, its Jacobian, and the three cases that were being conflated

`geometrodynamics/viz/congruence.py` · probe `experiments/closure_ledger/congruence_probe.py` ·
renderer `scripts/geometrodynamics_v50_congruence.py`

## Why not a bright point

A singularity in general relativity is a failure of evolution and geodesic
completeness. Infinite curvature is one possible signature, not the definition.
Drawing "a singularity" as a glowing dot assumes the answer before the geometry
has said anything.

The object that does not assume it is a **congruence with a deforming
cross-section**. Integrate the geodesic-deviation equation in TT gauge,

```
F̈ = ½ ḧ F ,    F(0) = I ,    J = det F
```

and watch `J`, the cross-sectional area of the bundle in units of its initial
area. `J → 0` is a caustic **of the map** — neighbouring trajectories arrive at
the same place, several histories land on one rendered point. It says nothing
about the metric.

For the axisymmetric field `h = diag(h₊, −h₊)` the gradient stays diagonal in
the geodesic-polar frame, `F = diag(a, b)`, so `a` is the principal stretch
along `ê_d` and `b` the one along `ê_ψ`. Those two are the only lines the
renderer draws.

## The focusing is entirely shear

The transverse Raychaudhuri equation is *the same statement* as the deviation
equation here:

```
dθ/dt = −θ²/2 − 2σ² − R_ab u^a u^b        residual 6.7e-15
```

with `θ = d(ln J)/dt`. That residual is an **algebraic identity holding to
round-off, not an accuracy check** — substituting `θ = A + B` and
`σ = (A − B)/2` into the right-hand side gives exactly `−A² − B²`, which is what
the deviation equation produces, so it is zero symbolically. The `n_radial`
ladders elsewhere are what test the integration.

The content is elsewhere: the **Ricci term is identically zero**, because `h` is
trace-free — `ä/a` and `b̈/b` are exact opposites, so the matter term cancels
term by term rather than approximately. None of the focusing is matter; all of
it is shear-squared, which is *second order* in the amplitude. That is the shape
of every threshold below: a weak wave barely focuses and a strong one does so
abruptly.

## Two things that had to be got right first

### `ḧ` from the wave operator, not from a time difference

`h = sin²d · q` and `q` solves the spin-2 wave equation, so `ḧ = sin²d · ∇²q` is
available exactly at every step. The first draft instead formed `ḧ` from a
three-sample time difference seeded with `h(−dt) = h(0)`. That seeding injects a
spurious velocity impulse `½ḣ(0)` at the first step — and because the `1/dt²` of
the difference cancels against the `dt` of the update, **the error does not
shrink with the timestep**.

| `n_radial` | operator `min J` | seeded difference `min J` |
| --- | --- | --- |
| 600 | −13.055 | 0.6273 |
| 1200 | −12.906 | 0.6279 |
| 2400 | −12.836 | 0.6281 |

Both columns converge. They converge to different answers, and the wrong one
reports no caustic at all. A refinement ladder — the usual defence in this
series — cannot see this class of error; only the operator form can.

### A compactly supported launch, or the far side moves early

A Gaussian bump is nonzero everywhere, so the congruence at the far side starts
deforming immediately. Probing at `d = 2.9499`, where the causal bound is
`d − w = 2.7699`:

| `n_radial` | Gaussian arrival | compact arrival |
| --- | --- | --- |
| 800 | 1.9992 | 2.6476 |
| 1600 | 2.0008 | 2.6976 |
| 3200 | 2.0011 | 2.7278 |

The Gaussian's arrival is **grid-converged** three digits inside the bound: an
analytic tail, not a numerical precursor, and refining cannot help. The compact
bump `(1 − x²)⁴` lands at `2.7697` against `2.7699` at the 1e-6 level, with its
residual earliness shrinking under refinement. `wave_constraints` found this for
the scalar; the spin-2 launch needed it too, and it is the default.

The same care is needed at the poles. The solver's grid is cell-centred, so its
last sample sits at `π − dd/2`; interpolating `h` itself *clamps* there instead
of vanishing. `h` and `ḧ` are therefore built as `sin²d ·` (interpolated `q`),
with the `sin²d` evaluated on the slice's own distances — so the spin-weight
condition holds exactly at the sampled points, at the one place this round
claims nothing is driven at all.

## The three cases, separated

1. **Ordinary focus.** `J` dips and recovers; the map stays invertible.
2. **Caustic.** `J` reaches zero; trajectories cross; the background is regular.
3. **Curvature singularity.** The geometry itself would have to fail.

**Case 3 is not reachable in this program, and that is a statement about the
program rather than about the wave.** The background of the whole series is a
fixed round `S²` with Gaussian curvature `1` everywhere at every time. There is
no Einstein equation in the loop and no backreaction, so there is nothing that
could diverge or terminate. A caustic is the strongest thing available here, and
reading one as a singularity is exactly the conflation the congruence was built
to prevent.

## Two rings, and they are not alike

The source ring and the converging ring are different events with different
thresholds, measured in a `1.2π` window:

| ring | threshold peak strain |
| --- | --- |
| source | 0.026 |
| converging (`d > π/2`) | 0.247 |

A factor of ten. The source is shaken hardest and from `t = 0`; the converging
ring has to arrive before it can focus. Closing the antipodal neck costs a
strain nothing physical would reach — and even at that strain the crossing only
**grazes** zero: it crosses and returns within a few thousandths, and the depth
of that excursion does not converge. That is the honest answer to "does focusing
drive the neck radius to zero": barely, and only just.

### Every threshold carries its window

The wave on `S²` is exactly periodic, so a fixed material point is driven over
and over by the same returning ring. The deviation equation is then a Hill
equation and the accumulated focusing keeps growing:

| window | converging-ring threshold strain |
| --- | --- |
| 1.2π | 0.255 |
| 2.0π | 0.116 |
| 3.0π | 0.082 |
| 4.0π | 0.063 |

A threshold quoted without a window is not a number about the wave. This is a
genuinely different mechanism from the single-pass refocus and the two must not
be reported as one figure.

## The neck is a ring, and spin weight is why

`h = sin²d · q` vanishes at both poles no matter what `q` does, so the tidal
field is identically zero **at** the antipode and the congruence there is never
driven. The thinnest cross-section therefore sits on a *ring* around the focus,
at geodesic radius `≈ 0.44 w` — the same ratio across a 3.3× range of pulse
width, to within the one-cell resolution of the angular grid, and independent of
amplitude. A spin-2 focus has no centre.

Three things this measurement had to get right, each of which it got wrong
first:

* **fixed peak strain, not fixed gain.** A narrower pulse carries more amplitude
  for the same gain, so a fixed-gain comparison varies two things at once and
  produced ratios scattered from 0.79 to 1.21.
* **searched in the antipodal cap, not in `d > π/2`.** Before the refocus the
  far half's minimum sits *exactly* on the cut at `d = π/2` — an artefact of
  where the region was truncated, not a neck. It reported `r/w = 8.7` and up.
* **timed to the refocus window.** `J` decreases secularly under the returning
  ring, so a minimum taken over a whole run is just wherever the run stopped.
  Between `t = π` and `1.3π` the location is stable; after `1.5π` accumulation
  takes over and it drifts.

## What the equations decided

Of the three outcomes on offer — passage, singular termination, or
finite-radius reconnection into a detached resonator — this program gives
**passage**, measured where the crossing is unambiguous. At the source ring `J`
crosses zero with slope `−17.877`, converged to `−17.836` under a halved
timestep (0.2%), and plunges to `−471` and stays. A tangency would have driven that
slope to zero; this one converges to a definite nonzero value. The solver's own
invariant is unmoved at `2.5e-14` across the crossing.

The other two were never available, and for different reasons worth keeping
apart:

* **Termination** needs the geometry to fail, and the background here is fixed
  and round with curvature `1` at every time.
* **Reconnection** needs the congruence to act back on something, and each
  material point's `F` is driven only by the external `h`, never by its
  neighbours.

So this is not "we looked and did not find them". It is "this program could not
have produced them", which is a different and more useful statement — and it
says exactly what a program that *could* produce them would need: metric
evolution with backreaction, which is the next thing to build rather than the
next thing to claim.

## Scope

`gain` is a strength dial and is reported as a peak strain everywhere it
matters. The deviation equation is exact in `ξ` and linear in `h`; at the strain
where the converging ring closes, the field is no longer a weak perturbation,
and that is stated rather than hidden. What is derived is the separation of the
three cases, the two thresholds with their window, the ring law, and the
passage.

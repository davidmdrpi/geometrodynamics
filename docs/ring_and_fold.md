# The ring reaches across, but only a fold crosses

> A **representation** of the scalar field on a fixed round `S²`, drawn in a
> vacuole. The tangential mixing `λ` is a modelling choice and is reported as one.

The scalar wave is not only its focal pulse. A ring leaves the source, thins to
a minimum at the equator, and then **grows** as it converges on the antipode. So
it is fair to ask whether that ring — rather than the pulse — is what first
reaches across the vacuole, and whether shrinking the gap or raising the energy
can make it intersect.

Those are two questions with different answers, and the useful result is that
they are controlled by **different knobs**.

## The ring is real, and it grows

| quantity | value |
|---|---:|
| equator height (`d = 1.507`) | `0.156` |
| peak height at the focus | `0.933` |
| growth | `5.97×` |
| against `1/√(sin d)` | `1.034 ± 0.071` |

All of that happens *before* the focal pulse. The law is the geometric factor
for a ring whose circumference is closing, and it holds until `sin d < 0.25` —
inside a pulse width of the antipode the ring is no longer thin compared with
its own radius and the law stops applying, which is reported rather than fitted
away.

## Reaching across — yes, and the ring gets there first

The threshold is exactly `δ / max|u|`, so raising the energy and shrinking the
gap both buy it. Shrinking the gap buys something extra: **lead**, the head start
the converging ring gets on the focal pulse.

| δ | ε | dwell | spans from `d` | lead |
|---:|---:|---:|---:|---:|
| 0.26 | 0.40 | 0.062 | 3.142 | 0.209 |
| 0.16 | 0.80 | 0.371 | 2.437 | 0.811 |
| 0.09 | 0.80 | 1.000 | 1.832 | 1.413 |
| 0.03 | 0.20 | 0.771 | 1.832 | 1.413 |

At the bottom of that table the ring is spanning the gap from just past the
equator, for the whole converging leg — not an instant at the focus but a
sustained state. **The intuition is right: the ring reaches across, and it does
so long before the pulse.**

## ...and still never intersecting

Swept over gap and gain with a real segment-intersection test:

| δ | ε | seam crossings | self-intersections |
|---:|---:|---:|---:|
| 0.26 | 10.0 | 38 | **0** |
| 0.09 | 10.0 | 114 | **0** |
| 0.03 | 10.0 | 346 | **0** |

A curve `r = f(σ)` with `f` single-valued puts exactly one radius at each
direction, so it is **embedded** by construction: two of its points cannot occupy
the same place. This is the winding obstruction seen from the side — a graph can
no more cross itself than it can wind — and neither knob touches it.

The detector is validated first, against a circle (`0`), a limaçon with an inner
loop, a lemniscate, and a folded loop (all `> 0`), so the zero is a result rather
than a broken counter.

## What does cross: tangential freedom

Let each material point move **sideways** as well as outward:

```
σ(σ₀) = σ₀ + λ ε ∂_σ u          r(σ₀) = R_mid · exp(ε u)
```

The map `σ₀ ↦ σ` folds where `∂σ/∂σ₀ = 1 + λ ε ∂²_σ u < 0`.

| quantity | value |
|---|---:|
| predicted `λε = 1/max(−∂²_σu)` | `0.012692` |
| found by bisection | `0.012692` |
| relative error | `1.8e-12` |
| folds first, from the antipode | `0.0157` |
| convergence drift under refinement | `6.6e-03` |

`λ = 0` is exactly the height field of v46, embedded forever. Past the threshold
the curve is multivalued in `σ`, stops being a graph, and self-intersects.

**And it folds first on the converging ring** — `0.0157` from the antipode, at
the moment of tightest focus, because `∂²_σ u` peaks at the steep front. The
intuition was right about *where* to look.

### Folding is necessary, not sufficient

| | count |
|---|---:|
| crossing without a fold | **0** |
| fold without a crossing | 8 |

A crossing always comes with a fold — that direction is the embedded-graph
theorem again. A fold need not cross, because the two branches of a folded map
can stay radially apart and never meet in the plane.

## The two knobs are orthogonal

| δ | span threshold | fold threshold |
|---:|---:|---:|
| 0.26 | 0.260 | `0.012461` |
| 0.12 | 0.120 | `0.012461` |
| 0.05 | 0.050 | `0.012461` |

**The fold threshold does not know the gap exists** — spread `0.0` across a
fivefold range of `δ`, while the spanning threshold scales with `δ` directly.

What it does scale with is the pulse:

| pulse width | `λε` | `λε / w²` |
|---:|---:|---:|
| 0.36 | `0.048944` | 0.378 |
| 0.24 | `0.021995` | 0.382 |
| 0.18 | `0.012692` | 0.392 |
| 0.12 | `0.005592` | 0.388 |

`λε ≈ 0.385 w²`, spread `3.7%`. Narrow fronts fold sooner, because folding is
about how sharply the front *turns*, not how tall it is.

**So reducing the distance between the shells changes when the wave arrives at
them. It changes nothing about whether it can cross itself.** That is the answer
to the question: more energy and a tighter gap both buy reaching, and neither
buys an intersection — only tangential freedom does.

## Scope

`λ` is a modelling choice, not derived from the scalar equation: it says how much
of the displacement is along the surface rather than across it. What *is* derived
is the threshold given `λ`, its independence from the gap, and its `w²` scaling.

The fold threshold depends on `∂²_σ u` at a focusing front, exactly the quantity
a coarse grid gets wrong, so it is converged by refining the angular sampling and
the radial solve **together** — refining `σ` alone does not converge, it samples
interpolation noise in a solve that has not kept up.

## Reproduce

```bash
python -m experiments.closure_ledger.ring_and_fold_probe
# Verdict: THE_RING_REACHES_BUT_ONLY_A_FOLD_CROSSES  (8/8)

python scripts/geometrodynamics_v48_ring_and_fold.py                 # animate
python scripts/geometrodynamics_v48_ring_and_fold.py --still out.png

python -m pytest tests/test_viz_slice_folding.py -q                  # 17 passed
```

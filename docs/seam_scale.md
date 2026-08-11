# The scaling at the seam is a choice

> Both gluings are **representation** choices, not derived boundary conditions.
> Nothing here makes the wave dynamically aware of the seam.

`docs/circle_slice_bulk.md` folded the bulk by **translation**: a radius past
`R_outer` came back as `r − gap`. That carries a radial *offset* across
unchanged — and the two boundary circles do not have the same circumference. A
feature emerging at `R_inner` keeps its full radial height while sitting on an
arc shorter by `R_outer/R_inner = 1.7027`.

**The emerging wave was not the same wave**, and the scaling that made it so was
never chosen deliberately. This picks the choice apart.

## Two rules

| | map | translation in |
|---|---|---|
| `translate` | `r → r − gap` | `r` |
| `conformal` | `r → r · (R_inner/R_outer)` | `ln r` |

They agree to first order in the excursion — their difference grows in exact
proportion to it, slope `0.4127` — so the earlier pictures were not wrong where
the wave was small. They part company exactly where the wave reaches the seam.

## What the choice changes

### The shape of what emerges

| rule | height in | height out | arc out | aspect distortion |
|---|---:|---:|---:|---:|
| `translate` | 0.1000 | 0.1000 | 0.0370 | **1.7027** |
| `conformal` | 0.1000 | 0.0587 | 0.0370 | **1.0000** |

The conformal rule scales the offset by the same factor as the boundary, so
height and arc length shrink together and the feature comes back a **faithful
scaled copy**. The translate rule pays the circumference ratio in full and
returns a caricature.

### Whether the radius survives at all

| rule | sheet −1 | −2 | −3 | −4 |
|---|---:|---:|---:|---:|
| `translate` | 0.220 | **−0.300** | **−0.820** | **−1.340** |
| `conformal` | 0.435 | 0.255 | 0.150 | 0.088 |

Arithmetic sheets march straight through `r = 0` into negative radius —
subtracting a fixed `gap` has nothing to stop it. Geometric sheets shrink by a
constant factor and accumulate at the origin without ever arriving.

This forces a matching choice for how the field makes a radius:

| law | | |
|---|---|---|
| `additive` | `r = R_mid + εu` | can be driven negative |
| `multiplicative` | `r = R_mid·exp(εu)` | positive by construction |

The conformal seam pairs with the multiplicative law. Ask for a conformal seam
with an additive law and it will raise rather than quietly take the logarithm of
a negative number.

### What a winding number would look like

Take a curve that genuinely winds — a ramp climbing exactly one period as `σ`
goes around, which on the conformal seam is a **logarithmic spiral**:

| seam | magnifications from `r₀ = 0.80, 1.00, 1.20` | spread |
|---|---|---:|
| `translate` | 1.6500, 1.5200, 1.4333 | `2.2e-01` |
| `conformal` | 1.7027, 1.7027, 1.7027 | `2.2e-16` |

A scale is start-independent; a shift is not. On the conformal seam a curve of
winding `w` returns to the same point of the quotient **magnified by
`(R_outer/R_inner)^w`** — two turns give `2.899 = 1.7027²`. So a conformal
gluing turns topological charge into an *observable magnification*. The
translate gluing hides it.

## What the choice does not change

The winding number. Rebuilt on a conformal seam with a multiplicative radial law
— a different identification, a different sheet structure, a different notion of
size — it is still identically zero:

| seam | radial law | unsigned | signed | winding |
|---|---|---:|---:|---:|
| `translate` | `additive` | 274 | `+0` | `+0` |
| `conformal` | `multiplicative` | 268 | `+0` | `+0` |
| `conformal` | `additive` | 36 | `+0` | `+0` |

`ρ(σ)` comes from a single-valued function on the circle **whichever coordinate
the seam translates in**, so its degree is zero either way. That is worth
having: the earlier negative result was not an artefact of an arbitrary scaling,
and no choice of gluing rescues it.

What the conformal rule adds is that the winding you cannot have would have been
*visible* — as a factor of `1.703` per turn. The obstruction is unchanged and its
cost is now legible.

## Scope

The conformal rule is preferred here on grounds of **consistency** — it returns a
scaled copy and cannot produce a negative radius — not because any dynamics picks
it out. Neither rule makes the seam something the wave can feel, and neither
gives a height field anywhere to wind. `translate` remains the default so v46 is
reproducible.

## Reproduce

```bash
python -m experiments.closure_ledger.seam_scale_probe
# Verdict: THE_GLUING_SETS_THE_SCALE_BUT_NOT_THE_TOPOLOGY  (6/6)

python scripts/geometrodynamics_v47_seam_scale.py                 # animate
python scripts/geometrodynamics_v47_seam_scale.py --still out.png

python -m pytest tests/test_viz_seam_scale.py -q                  # 15 passed
```

The renderer shows one wave at one gain folded both ways side by side, the
emerging feature close up under each rule, a curve that does wind coming back
magnified, and the sheet ladder with `r = 0` marked.

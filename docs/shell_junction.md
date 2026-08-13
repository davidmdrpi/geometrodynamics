# Can a detached shell do the throat's exotic work for it?

`geometrodynamics/shells/junction.py` · probe `experiments/closure_ledger/shell_junction_probe.py`

## Scope of every claim below

**Einstein gravity, Darmois–Israel thin shells, spherical symmetry, vacuum
bulk.** `G = 1`. Nothing here bounds thick shells, non-spherical
configurations, modified gravity, or non-vacuum interiors. The dimension is a
parameter: `D = 4` is the regression case and `D = 5` (Tangherlini) is the one
this program actually cares about.

## The hope, and why it needs three answers

A wormhole throat needs negative surface energy. The hope is that a *detached*
closed shell in the bulk, glued with the opposite orientation, might supply the
exotic-looking restoring stress while itself being ordinary matter.

That is three questions wearing one coat, and they do not agree:

| | observable | answer here |
| --- | --- | --- |
| 1 | does the shell itself require exotic matter? `σ = −S^τ_τ` | **yes**, if it connects |
| 2 | does it support the throat? the gradient of `ΔV` | **yes**, and from ordinary matter |
| 3 | is it stable? the stiffness `V''(b₀)` | **yes**, with a negative `β²` at the throat |

## The orientation is derived from the gluing, not assumed

An earlier draft carried `ε = ±1` as a free flag, which made the central result
conditional on a hand-set sign. Here each side of the surface is a **branch** —
which radial half of its region is retained — and `ε` follows:

* the `−` side is approached from within along the `−→+` normal, so `ε₋ = +1`
  for INNER (`r ≤ R`) and `−1` for OUTER (`r ≥ R`);
* the `+` side is entered leaving the surface, so `ε₊ = +1` for OUTER and `−1`
  for INNER.

With `σ = −(D−2)(ε₊β₊ − ε₋β₋)/(8πG R)` there are therefore **four** gluings,
not two:

| `−` branch | `+` branch | `η = ε₊ε₋` | what it is | `σ` |
| --- | --- | ---: | --- | --- |
| INNER | OUTER | `+1` | ordinary bubble | either sign |
| OUTER | OUTER | `−1` | **minimal surface** | `< 0` always |
| INNER | INNER | `−1` | **maximal surface** | `> 0` always |
| OUTER | INNER | `+1` | anti-bubble | either sign |

**`η = −1` alone decides nothing.** It covers two gluings whose forced signs are
*opposite*. The earlier framing of the result in terms of `η` was too coarse;
the sign is a property of the branch pair.

The machinery is checked against published answers first: a flat-interior bubble
in Schwarzschild carries ordinary matter whose rest mass is the bulk mass to
`1e-3`, and a `Z2` throat reproduces Visser's `σ = −√f/(2πGR)` to `6.9e-18`.

## What is actually forced

For a **minimal surface** — `r` increasing away on both sides, which is what a
throat is — the two terms add:

```
σ = −(D−2)(β₊ + β₋)/(8πG R) < 0
```

with `β± = √(f± + Ṙ²) ≥ 0` for any timelike shell. For a **maximal surface**
they add with the other sign and `σ > 0` always. Both are identities in every
`D`, and neither is violated once in 40,000 random Tangherlini / de Sitter /
charged pairs across `D = 4, 5, 6`. The sweep checks the implementation; the
claims are identities.

### The dichotomy that follows

This is the real result, and it is sharper than "the oppositely-glued shell is
exotic":

> A detached surface that **connects** to the throat's asymptotic region does so
> through a minimal surface, and is necessarily exotic. A detached surface that
> is non-exotic by its gluing is a **maximal** surface, which caps off on both
> sides and therefore shares no bulk with the throat — it is non-exotic
> precisely because it is disconnected, and so cannot support anything.

Within Einstein–Israel spherical thin shells, exotic matter is relocated, never
removed.

## Observable 2, as a gradient rather than a force

An ordinary bubble screens mass, shifting the throat's potential by `2GΔμ/b`,
with `−∂ΔV/∂b = +0.024` at screened `μ = 0.6` and the acceleration contribution
half that, from `ḃ² + V = 0`.

**This is not an equilibrium-consistent force.** The continuation holds the
throat's rest mass fixed and omits its equation-of-state response, so it is the
gradient of the potential *shift* screening produces. It is reported and named
that way.

## Observable 3, as a stiffness

`ḃ² + V(b) = 0`; a static solution needs `V(b₀) = 0` and `V'(b₀) = 0` and is
stable when `V''(b₀) > 0`. A fixed global barotropic index admits no stable
static throat at all — `V''(b₀) = 2GM(n−1)/b₀³` with `n = 2 + 4w` in `D = 4`, so
static needs `w < −1/2` and stability `w > −1/4` — so `β² ≡ (dp/dσ)|₀` is left
free at the equilibrium, as usual.

`β²` is an **equation-of-state derivative parameter**. No sound-speed reading is
attached to its sign, which for exotic matter would not be meaningful.

Screening raises the critical `β²` monotonically:

| interior `μ` | `β²_crit` |
| --- | ---: |
| 2.0 | −1.0833 |
| 1.8 | −0.9463 |
| 1.6 | −0.8439 |
| 1.4 | −0.7648 |
| 1.0 | −0.6518 |

So the shell **does** enlarge the window, and never reaches `β² ≥ 0`.

Both stiffnesses are verified against direct RK4 integration of
`σ' = −(D−2)(σ + p)/R`, agreeing to 10 digits.

**They are stiffnesses, not normal-mode frequencies.** No kinetic metric has
been derived, so nothing here is a generalised eigenproblem, and the code names
them accordingly.

## Birkhoff, and what it is worth

The region between the two surfaces is Schwarzschild/Tangherlini with a constant
mass parameter, so the throat's data cannot depend on where the shell sits.

The measured version: at fixed screened mass, moving the shell from `a = 8` to
`a = 200` changes its surface density by a factor of **701** — genuinely
different shells — and the throat's `σ` does not change in its last bit. Spread
exactly `0.0`.

The vanishing `∂²V/∂a∂b` is **structural, not measured**: Birkhoff is imported
the moment the intervening region is written that way, and reporting the zero as
evidence would be circular. What it establishes is **no separation-dependent
coupling in this model** — not that every spherical trapped resonator is
impossible.

The consequence for what comes next: spherical symmetry has no radiative channel
here, the same `ℓ = 0` fact `wave_constraints` found for the scalar, so `ℓ ≥ 2`
internal modes are where a genuine throat–shell coupling would have to live.

## The scaling check

Rescaling every length and mass parameter together sends the stiffnesses as
`1/L²` (drift `0.0`). That is **dimensional bookkeeping only** — it does not
show that a fixed system has no dilation mode, which would need the kinetic term
this module does not derive.

## What is imported rather than derived

* Birkhoff's theorem.
* The Darmois–Israel formalism itself.
* `β²` as a free parameter at the equilibrium.

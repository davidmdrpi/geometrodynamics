# Can a detached shell do the throat's exotic work for it?

`geometrodynamics/shells/junction.py` · probe `experiments/closure_ledger/shell_junction_probe.py`

## The hope, and why it needs three answers

A wormhole throat needs negative surface energy. The hope is that a *detached*
closed shell in the bulk, glued with the opposite orientation, might supply the
exotic-looking restoring stress while itself being ordinary matter — negative
matter without negative matter.

That is three questions wearing one coat, and they do not agree:

| | observable | answer here |
| --- | --- | --- |
| 1 | does the shell itself require exotic matter? `σ = −S^τ_τ` | **yes**, if it is oppositely glued |
| 2 | does it support the throat? `F = −∂V_shell/∂b` | **yes**, and from ordinary matter |
| 3 | is the configuration stable? `V''(b₀)`, normal modes | **yes**, but only with a pathological equation of state at the throat |

A positive stiffness alone would mean "restoring" and would not establish that
the *shell* supplied the support, so all three are reported separately.

## The orientation is not a relabelling

For a shell at areal radius `R`, `K^θ_θ = (ε/R)√(f(R) + Ṙ²)` with `ε = ±1` the
sign of `dr` along the chosen normal. Taking both normals to point from the `−`
side to the `+` side,

```
σ = −(1/4πG R)·(ε₊β₊ − ε₋β₋) ,    β± = √(f±(R) + Ṙ²)
```

`η ≡ ε₊ε₋` is a statement about which manifold was built:

* **`η = +1`, aligned** — an ordinary bubble. `r` increases outward on both
  sides; the surface separates an inside from an outside.
* **`η = −1`, anti-aligned** — a *minimal surface*. `r` increases as you move
  away on **both** sides, so the surface is a throat in its own right and the
  region beyond is a second asymptotic region.

The machinery is checked against published answers before anything new is asked
of it: a flat-interior bubble in Schwarzschild carries ordinary matter whose
rest mass is the bulk mass to `1e-3`, and a `Z2` throat reproduces Visser's
`σ = −√f/(2πGR)` to `6.9e-18`.

## Observable 1 is a theorem

For `η = −1` the two terms **add**:

```
σ = −(β₊ + β₋)/(4πG R) ≤ 0
```

and `β± ≥ 0` for any timelike shell. So **every** anti-aligned shell carries
negative surface energy — whatever the bulk on either side, whatever its mass,
charge, cosmological constant or velocity.

Swept over 200,000 random Schwarzschild / de Sitter / Reissner–Nordström pairs:
**zero counterexamples**, worst `σ = −6.7e-04`. The aligned control is positive
`50.1%` of the time, so the sweep is capable of finding a positive `σ` when one
exists. But the sweep checks the *implementation* — the claim itself is an
identity, not a statistic.

**The same identity applies to the throat**, because a throat *is* a minimal
surface. No arrangement of bulk content can relieve it. So the answer to
"negative matter without negative matter" is **no**, and for a reason that is
topological rather than dynamical: an oppositely-glued detached shell does not
take the exotic matter away from the throat, it adds a second helping.

## What an ordinary shell can do

An **aligned** shell can be perfectly ordinary. Same regions, same radius, only
the gluing differs: `σ = +6.2e-05` aligned against `−9.5e-02` anti-aligned.

And it does act on the throat, by **screening mass**. With ADM mass `M₂` outside
and `M₁` between shell and throat, the shell's presence shifts the throat's
potential by `2G(M₂ − M₁)/b`, an outward force

```
F_shell = 2G ΔM / b²
```

matched to `1e-6`, growing with the screened mass, and exactly zero when there
is no shell. That is real support, from non-exotic matter.

## Observable 3, and what screening buys

Following the standard thin-shell analysis, `ḃ² + V(b) = 0`; a static solution
needs `V(b₀) = 0` and `V'(b₀) = 0`, and is stable when `V''(b₀) > 0`.

Fixing a global barotropic index `w` turns out to be too rigid to be
interesting. The algebra gives `V''(b₀) = 2GM(n−1)/b₀³` with `n = 2 + 4w`, so a
static solution exists only for `w < −1/2` and is then *always* unstable — the
existence and stability conditions are disjoint. So `β² ≡ (dp/dσ)|₀` is left
free at the equilibrium instead, which is the usual treatment, and `V''` is
linear in it.

Both stiffnesses — throat and shell — are verified against direct RK4
integration of the conservation law `σ' = −(2/b)(σ + p)`, to `~1e-6` relative,
on both signs of `β²`.

Screening raises the critical `β²` monotonically:

| interior mass | `β²_crit` |
| --- | ---: |
| 1.0 | −1.0833 |
| 0.9 | −0.9463 |
| 0.8 | −0.8439 |
| 0.7 | −0.7648 |
| 0.5 | −0.6518 |

So the shell **does** enlarge the throat's stability window. It never reaches
`β² ≥ 0`: the throat always needs `dp/dσ < 0`, an imaginary sound speed, on top
of its negative energy density.

Both normal modes can nevertheless be positive at once — `diag(0.151, 0.022)`
for a throat at `β² = −2` and an ordinary shell at `β² = +0.5`. Under a dilation
the spectrum scales as `1/L²` exactly and the smallest scaled eigenvalue stays
at `5.4e-04`, so that is not an artefact of whatever fixed the scale.

## The finding that shapes what comes next

**Birkhoff decouples the two surfaces exactly.** The vacuum between them is
Schwarzschild with a constant mass parameter, so the throat cannot tell *where*
the shell is — only that it is there.

The non-trivial version of that statement, which is what gets measured: at fixed
screened mass, moving the shell from `a = 8` to `a = 200` changes its surface
density by a factor of **701** — genuinely different shells, different rest
mass, different pressure, different stiffness — and the throat's `σ` does not
change in its last bit. Spread exactly `0.0`.

The off-diagonal Hessian entry `∂²V/∂a∂b` vanishes, but that is **structural,
not measured**: Birkhoff is imported the moment the region between is written as
Schwarzschild with a constant mass, and reporting the resulting zero as evidence
would be circular. It is labelled as such in the module and in the probe.

The consequence is what matters. Spherical symmetry has no radiative channel —
the same `ℓ = 0` fact `wave_constraints` found for the scalar — so two spherical
surfaces have nothing to talk through. **A genuine two-mode trapped resonator
cannot exist in spherical symmetry at all.** The `ℓ ≥ 2` internal modes are not
an optional later refinement; they are the only place such a coupling could
live.

## Scope

`G = 1`. Spherical symmetry, thin shells, static configuration plus linear
perturbation. No radiation, no backreaction on the shells' equations of state,
no `ℓ ≥ 2` structure. `β²` is a free parameter at the equilibrium in the
standard way rather than a global equation of state.

What is *derived*: the exoticity theorem for minimal surfaces, the screening
force, the stability window and how screening moves it, and the exact
decoupling. What is *imported*: Birkhoff's theorem, and the Israel junction
formalism itself.

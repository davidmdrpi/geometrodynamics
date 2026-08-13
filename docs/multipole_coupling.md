# ℓ ≥ 2: the coupling Birkhoff forbids at ℓ = 0, and what it costs to have it

`geometrodynamics/shells/multipole.py` · probe `experiments/closure_ledger/multipole_coupling_probe.py`

## 0. The answer, stated first

`shell_junction` found that two concentric surfaces in a vacuum spherical model
cannot talk, and **imported Birkhoff's theorem** to say why. That import was not
needed. For two concentric shells at `b < a`, each deformed by `δR = α R P_ℓ`,
the mutual stiffness of the gravitational interaction is

```
∂²U/∂α∂γ  =  G m_b m_a · ℓ(ℓ+1) · (b/a)^ℓ / (a (2ℓ+1)²)
```

The prefactor is `ℓ(ℓ+1)` — **the eigenvalue of the Laplacian on the sphere**.
So the decoupling of the previous round is not a special fact about spheres, and
not a separate theorem: it is the `ℓ = 0` case of a one-line multipole
statement, and it vanishes because the constant mode has zero eigenvalue.

The answer has a **second half that matters as much as the first**: the same
formula screens the coupling as `(b/a)^ℓ`. At `b/a = 0.4` the mutual stiffness
falls by a factor of **544** from `ℓ = 1` to `ℓ = 8`. The multipoles that carry
a spin-2 signal are precisely the ones separation suppresses hardest, so having
the channel and being able to use it are different statements.

## Scope

This is the **static Newtonian (Laplace) multipole problem** — the weak-field
limit of the junction problem, in which the interior- and exterior-regular
static solutions are `r^ℓ` and `r^{−(ℓ+D−3)}`. It is **not** a
Regge–Wheeler/Zerilli treatment on the Schwarzschild background, and nothing
here claims an `ℓ ≥ 2` quasinormal frequency, a radiative rate, or that any
particular shell would resonate. `G = 1`.

What is established is the multipole structure of the coupling and where
Birkhoff sits inside it.

## The verification that carries the result

The closed form is checked against **direct double integration over both
deformed surfaces**, which never expands in multipoles and so is free to
disagree:

| `ℓ` | closed form | brute force | rel. error |
| ---: | ---: | ---: | ---: |
| 0 | 0.000000000 | 0.000000000 | — |
| 1 | 0.017777778 | 0.017777863 | 4.8e-06 |
| 2 | 0.007680000 | 0.007680044 | 5.7e-06 |
| 3 | 0.003134694 | 0.003134716 | 7.0e-06 |
| 4 | 0.001264198 | 0.001264209 | 9.0e-06 |

at `b = 2`, `a = 5`. The two shells never touch, so there is no
coincident-point singularity to regulate — the control is clean.

## A trap that had to be caught first

**A pure `P₁` deformation is not a translation past linear order.** The second
variation of the area of `r = R(1 + αP_ℓ)` is

```
d²A/dα² / (4πR²)  =  [2 + ℓ(ℓ+1)] / (2ℓ+1)
```

exactly — `2, 4/3, 8/5, 2, 22/9, 32/11` for `ℓ = 0…5`, each matching its
rational to `1e-12` and the direct integration to `3.6e-08`. At `ℓ = 1` this is
`4/3`, **not zero**. A translation-invariance check built on a pure `P₁`
deformation would have reported a stiffness that is not there.

The resolution: a rigid displacement `d` is

```
r = R + d P₁ − d²/(3R) + (d²/3R) P₂ + O(d³)
```

so the true translation direction carries induced `ℓ = 0` and `ℓ = 2` pieces,
and those are what preserve the area. Three ways of asking separate cleanly:

| construction | `ΔA/A` at `d = 0.02` |
| --- | ---: |
| exact displaced sphere | `4e-16` |
| `O(d²)` family with induced `ℓ = 0, 2` | `2.1e-08` |
| pure `P₁` | `2.7e-04` |

The exact sphere is area-preserving to round-off at every displacement; the
truncated family is preserving to `O(d⁴)`, beating the pure dipole by `O(d²)` as
it should. **Translation invariance is exact — the naive test was not testing
it.**

## The shear response is an input, not a derivation

A perfect-fluid surface layer is `S_ij = diag(−σ, p, p)`: it resists **area**
change and nothing else. Its `ℓ ≥ 2` restoring force comes only from that area
cost and from gravity. Making a shell resist *shape* change at fixed area
requires an elastic modulus `μ_s` that no equation of state supplies, so it is
carried as an explicit parameter and never fitted.

This is stated because the previous round's conclusion — that `ℓ ≥ 2` is where
the coupling lives — is only half an answer. The coupling is there. What a shell
*does* with it depends on a constitutive choice that spherical symmetry never
had to make.

The elastic shear stiffness carries `(ℓ−1)(ℓ+2)`, which vanishes at `ℓ = 1`:
shape response cannot restore a translation, as it should not.

## What this closes, and what it does not

**Closes:** why the spherical model decoupled, in a form that needs no imported
theorem, and with the `ℓ`-dependence made explicit in both directions — the
coupling exists for every `ℓ ≥ 1` and is geometrically screened.

**Does not close:** the dynamics. No quasinormal frequency, no radiative rate,
no Regge–Wheeler treatment, and no claim that a trapped `ℓ ≥ 2` resonance
exists. The constitutive gap is named, not filled.

# Where the two-shell coupling starts, in a static Newtonian model

`geometrodynamics/shells/multipole.py` · probe `experiments/closure_ledger/multipole_coupling_probe.py`

## 0. Scope, because the headline depends on it

This is the **static Newtonian (Laplace) two-shell model** — the weak-field
analogue of the junction problem, with interior/exterior static solutions `r^ℓ`
and `r^{−(ℓ+D−3)}`. What it establishes is the **shell-theorem / multipole**
structure of the coupling.

**Birkhoff's theorem is a GR result and remains what `shell_junction` (PR #249)
relies on.** Nothing here replaces it; the `ℓ = 0` statement below is its
Newtonian analogue. Not Regge–Wheeler/Zerilli, no quasinormal frequencies, no
radiative dynamics. `G = 1`.

## 1. The answer, stated first

In this model the **monopole mutual stiffness vanishes** while higher angular
multipoles couple, with the coupling **suppressed geometrically by separation**.

For two concentric shells at `b < a`, each deformed by `δR = α R P_ℓ`:

```
∂²U/∂α∂γ  =  G m_b m_a · ℓ(ℓ+1) · (b/a)^ℓ / (a (2ℓ+1)²)
```

The prefactor is `ℓ(ℓ+1)`, the eigenvalue of the angular Laplacian — so **the
`ℓ = 0` decoupling is that zero eigenvalue**.

Verified against **direct double integration over both deformed surfaces**,
which never expands in multipoles:

| `ℓ` | closed form | brute force | rel. error |
| ---: | ---: | ---: | ---: |
| 0 | 0.000000000 | 0.000000000 | — |
| 1 | 0.017777778 | 0.017777863 | 4.8e-06 |
| 2 | 0.007680000 | 0.007680044 | 5.7e-06 |
| 3 | 0.003134694 | 0.003134716 | 7.0e-06 |
| 4 | 0.001264198 | 0.001264209 | 9.0e-06 |

## 2. Where the coupling starts — and it is not `ℓ = 1`

An earlier draft of this module concluded that *"everything `ℓ ≥ 1` couples"*.
That is **wrong as a statement about physical modes**, and the error was
checking translation invariance of the **area** rather than of the **mutual
energy**. Run on the energy, the two disagree:

| construction | mutual coupling |
| --- | ---: |
| **rigid translation** (exact displaced spheres) | `8.3e-13` — zero |
| pure `P₁` **shape** deformation | `1.78e-02` |

A rigidly displaced inner sphere leaves the mutual energy at exactly
`−G m_b m_a / a` — **Newton's shell theorem** — held to `1e-15` out to `d = 2.5`
while the inner stays entirely inside. The translation mode does **not** couple.
A pure `P₁` shape deformation is a different object, and it does.

So the honest ordering is:

* `ℓ = 0` decouples by the **vanishing eigenvalue**;
* the `ℓ = 1` **translation** mode decouples by the **shell theorem**;
* genuine coupling **starts at `ℓ = 2`** — which is what PR #249 guessed, and
  this establishes.

## 3. The same trap, twice

Both errors are the same mistake: **a pure `P₁` deformation is not a translation
past linear order.** The area second variation of `r = R(1 + αP_ℓ)` is

```
d²A/dα² / (4πR²)  =  [2 + ℓ(ℓ+1)] / (2ℓ+1)
```

exactly — `2, 4/3, 8/5, 2, 22/9, 32/11` for `ℓ = 0…5`, matching each rational to
`1e-12` and direct integration to `3.6e-08`. At `ℓ = 1` that is `4/3`, not zero.

A rigid displacement is instead
`r = R + dP₁ − d²/(3R) + (d²/3R)P₂ + O(d³)`, carrying induced `ℓ = 0` and
`ℓ = 2` pieces:

| construction | `ΔA/A` at `d = 0.02` |
| --- | ---: |
| exact displaced sphere | `4e-16` |
| `O(d²)` family with induced `ℓ = 0, 2` | `2.1e-08` |
| pure `P₁` | `2.7e-04` |

**The lesson: a zero-mode test has to be run on the quantity the claim is
about.** Translation invariance of the area does not decide whether `ℓ = 1`
couples; translation invariance of the energy does, and it says it does not.

## 4. The second half of the answer

The same formula **screens** the coupling as `(b/a)^ℓ`:

| `ℓ` | mutual stiffness | `(b/a)^ℓ` |
| ---: | ---: | ---: |
| 1 | 1.78e-02 | 0.400 |
| 2 | 7.68e-03 | 0.160 |
| 4 | 1.26e-03 | 0.0256 |
| 8 | 3.27e-05 | 0.000655 |

A factor of **544** from `ℓ = 1` to `ℓ = 8` at `b/a = 0.4`. The multipoles that
carry a spin-2 signal are precisely the ones separation suppresses hardest, so
having the channel and being able to use it are different statements.

## 5. The shear response is an input, not a derivation

A perfect-fluid surface layer is `S_ij = diag(−σ, p, p)`: it resists **area**
change and nothing else. Making a shell resist *shape* change at fixed area
requires an elastic modulus `μ_s` that no equation of state supplies, so it is
carried as an explicit parameter and never fitted.

PR #249's conclusion that `ℓ ≥ 2` is where the coupling lives is therefore only
half an answer. The coupling is there. What a shell *does* with it depends on a
constitutive choice that spherical symmetry never had to make.

The elastic shear stiffness carries `(ℓ−1)(ℓ+2)`, vanishing at `ℓ = 1`: shape
response cannot restore a translation, as it should not.

## 6. What this closes, and what it does not

**Closes:** where the coupling starts in this model, in both directions — the
monopole decouples by a vanishing eigenvalue, the translation mode by the shell
theorem, genuine coupling begins at `ℓ = 2`, and it is geometrically screened.

**Does not close:** the dynamics, and the GR statement. No quasinormal
frequency, no radiative rate, no Regge–Wheeler treatment, no claim that a
trapped `ℓ ≥ 2` resonance exists — and no replacement for Birkhoff, which
remains a GR theorem imported by PR #249.

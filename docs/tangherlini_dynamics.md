# The first evolved Einstein equations

**Module:** `geometrodynamics/tangherlini/dynamics.py`
**Probe:** `python -m experiments.closure_ledger.tangherlini_dynamics_probe` (7/7)
**Tests:** `tests/test_tangherlini_dynamics.py` (25)
**Renderer:** `scripts/geometrodynamics_v71_tangherlini_dynamics.py`

---

Every gravity result in this repository so far is stationary, weak-field or
linearized. `waves/backreaction.py` linearizes about the ESU; #176/#177 evolve a
field while the metric responds quasi-statically; `THESIS.md` has said, in the
same words across five rounds, that *"the strong-field endpoint (a horizon / a
resolved throat) is left for full numerical relativity."* Nothing in the tree
evolves the Einstein equations in time.

This round is the first that does — at the highest-symmetry 4+1 problem, in
horizon-penetrating coordinates. It gets two results and misses two targets, and
the miss is reported as a miss.

## 0 · the system, and the gauge

`D = 5`, spherical symmetry (`S³` angular sector), one minimally coupled massless
scalar. Vacuum is not an option: Birkhoff in `D` dimensions makes Tangherlini the
unique spherically symmetric vacuum, so there is nothing to evolve. The scalar is
the dynamical content.

Ingoing Eddington–Finkelstein with areal radius,

    ds² = −A(v,r) e^{2δ(v,r)} dv² + 2 e^{δ(v,r)} dv dr + r² dΩ_n² ,   n = D−2 = 3

Ingoing null cones `v = const` cross the future horizon, so the chart is
horizon-penetrating by construction — no excision and no singularity-avoiding
lapse is needed to reach it.

## 1 · derived, not recalled

The metric, connection, Ricci tensor and Einstein tensor are built in sympy for
**general `n`**, the massless-scalar stress tensor is added, and the independent
components come out. Three self-checks, all passing at **both** `n = 3` and
`n = 2` — the latter being the known `D = 4` system, which is what validates the
general-`n` derivation:

| check | `n = 3` (`D = 5`) | `n = 2` (the known `D = 4` system) |
|--|--|--|
| `rr` is the `δ` quadrature | yes | yes |
| `vr` is the `A` quadrature | yes | yes |
| Birkhoff: constant `m` solves `vr` | yes | yes |

What comes out:

    rr    →   ∂_r δ = (κ/n) · r · (∂_r φ)²
    vr    →   (r^{n−1} e^δ A)' = (n−1) r^{n−2} e^δ
    wave  →   2 r^n ∂_r(∂_v φ) + n r^{n−1} ∂_v φ + ∂_r(e^δ A r^n ∂_r φ) = 0
    vv    →   an independent equation containing ∂_v A — never used

**The `vr` result is the surprise.** It is not an ODE to be integrated alongside
`A`; it is an exact **quadrature**. Each slice costs three cumulative integrals
and no implicit solve, so the geometry is as accurate as the quadrature rather
than as accurate as an ODE stepper.

## 2 · the hierarchy, and why it takes no boundary condition

Given `φ(v, ·)` on one ingoing cone, with a regular centre: `δ` from the `rr`
quadrature (gauge-fixed by `δ(v,R) = 0`), `A` from the `vr` quadrature (fixed by
regularity at `r = 0`), then

    ψ ≡ ∂_v φ = −½ e^δ A ∂_r φ − (n/4) r^{−n/2} ∫₀^r s^{n/2−1} e^δ A ∂_s φ ds

fixed by regularity at `r = 0` again. Then `φ` is marched in `v` with RK4.

**No outer boundary condition exists and none may be imposed.** The three
constants are already spent — two on regularity, one on the gauge — so `ψ(v,R)`
is determined by the data on the cone. An early version froze `φ(v,R)`; the `vv`
residual at the outer edge then sat at `O(1)` and did not converge at all. That
was the diagnostic that found it.

## 3 · two exact solutions, before anything is claimed

**Tangherlini is an exact fixed point.** With `φ = 0` the `A` quadrature is
`∫2s ds = r²`, so `A = 1 − r_h²/r²` must come back at machine precision, not to
some tolerance that hides a wrong power of `r`:

| points | metric error | max abs `δ` | max abs `ψ` |
|--|--|--|--|
| 400 | `1.554e-15` | `0.0` | `0.0` |
| 1600 | `4.441e-15` | `0.0` | `0.0` |

**The exact flat-space mode is reproduced at second order.**
`φ = cos(ω(v−r)) J₁(ωr)/r` solves the flat `D = 5` wave equation exactly in these
coordinates, so `∂_vφ` is known in closed form and the `ψ` quadrature has
something real to be wrong about:

| points | flat-metric error | `ψ` relative error | rate |
|--|--|--|--|
| 400 | `8.9e-16` | `1.470e-03` | — |
| 800 | `4.4e-16` | `3.649e-04` | `2.010` |
| 1600 | `7.8e-16` | `9.078e-05` | `2.007` |
| 3200 | `1.3e-15` | `2.262e-05` | `2.005` |
| 6400 | `1.1e-15` | `5.644e-06` | `2.003` |

That is what pins the scheme's order; nothing else in the round does. It is
second and not fourth on purpose: the `r^{n/2}` weight from odd `n` leaves a
half-integer power in the origin quadrature, no integration by parts removes it,
and a fourth-order rule measures `2.5` there — so the whole scheme was taken to a
clean, uniform second order instead.

## 4 · the constraint the code never solves

The hierarchy **solves** `rr` and `vr` on every slice, so their residuals are
identically zero and testing them would be circular. `vv` is the one independent
component left over, it contains `∂_v A`, and the code never forms `∂_v A` for
any other purpose. Its residual therefore tests the evolution rather than
restating it:

| points | spacing | max abs `vv` residual | at radius | rate |
|--|--|--|--|--|
| 400 | `0.0501` | `1.5511e-04` | `4.01` | — |
| 800 | `0.0250` | `3.9070e-05` | `4.01` | **`1.989`** |
| 1600 | `0.0125` | `9.7862e-06` | `4.00` | **`1.997`** |
| 3200 | `0.0063` | `2.4478e-06` | `4.00` | **`1.999`** |

> **Constraint convergence: second order, `1.989 → 1.997 → 1.999`.**

This is the characteristic-scheme analogue of a Hamiltonian/momentum constraint
test, and the analogy is **stated rather than assumed** — it is not a
Hamiltonian/momentum pair, and calling it one would be wrong.

*One more thing this found.* A one-sided `∂_v A` capped the measured rate at
`1.05`; it is centred in `v` now.

## 5 · a regular centre forbids a trapped surface

The `vr` quadrature with a regular centre reads

    r^{n−1} e^{δ(r)} A(r) = (n−1) ∫₀^r s^{n−2} e^{δ(s)} ds

— a positive integrand over a positive interval. So

> **`A > 0` strictly, for `r > 0`, identically. No trapped surface can sit on a
> regular-centred ingoing null slice.**

The scan below is not the proof; it is the check that the code obeys the proof.
Four profile families are driven to amplitudes where `min A` falls to a few parts
in a thousand, and it never crosses:

| profile | amplitude | min `A` | at radius | trapped? |
|--|--|--|--|--|
| centred gaussian | 12 | `7.256e-02` | `1.041` | no |
| thin shell at `r = 2` | 5 | `1.513e-02` | `2.402` | no |
| **thin shell at `r = 2`** | **12** | **`5.627e-03`** | `2.469` | no |
| `r²e^{−r²/2}` | 12 | `2.703e-02` | `2.322` | no |
| oscillatory | 12 | `5.743e-03` | `2.575` | no |

**So horizon formation is not observable in this gauge with a regular centre**,
and the criterion has to be posed as the loss of central regularity rather than
as `A` changing sign.

*This is a statement about the chart, not about the physics.* Collapse still
happens; the trapped region is reached once the centre stops being regular, at
which point the slice carries a nonzero interior mass constant and this
quadrature no longer applies. It is the reason production characteristic-collapse
codes use **outgoing** null cones, or excise the centre.

*And the identity is exact while its representation is not.* `δ` is gauge-fixed
to zero at the outer edge and increases outward, so `e^δ ≤ 1` and underflows once
a slice spans more than about `700` e-folds — reached only at amplitudes far past
anything physical here. Numerator and denominator underflow together, so `A`
comes back `nan` at those points rather than as a wrong sign, which is the
failure mode that would actually matter. Tested for.

## 6 · a discrepancy found in passing, reported and not acted on

The Schrödinger-form potential for a minimally coupled massless scalar with
`ψ = r^{n/2}φ` is, at `n = 3`:

| | potential |
|--|--|
| derived here | `A[(ℓ(ℓ+2) + 3/4)/r² + (9/4) r_h²/r⁴]` |
| `tangherlini.radial.V_tangherlini` | `A[ℓ(ℓ+2)/r² + 3 r_h²/r⁴]` |
| difference | `3A²/(4r²)` — reproduced to `5.4e-16` |

**The flat limit settles which is which with no appeal to anything.** At
`r_h → 0` the regular solution is `φ = J_{ℓ+1}(ωr)/r`, so
`ψ = r^{1/2}J_{ℓ+1}(ωr)`, and Bessel's equation gives
`V → ((ℓ+1)² − ¼)/r² = (ℓ(ℓ+2) + 3/4)/r²` — the derived form, including the
`3/4`, matched to `4.3e-16`.

**Nothing was changed.** `V_tangherlini` is consumed by roughly fifty probes and
by several derived constants; replacing it is a decision about the repository's
published numbers and not a side effect of a dynamics round. *Caveat:* the
discrepancy is stated for a **minimally coupled massless scalar with
`ψ = r^{n/2}φ`**, which is what the existing docstring describes; a different
field or a different substitution has a different potential.

## 7 · what this round did not earn

| asked for | delivered? |
|--|--|
| constraint convergence | **yes** — second order, `1.999` |
| horizon formation | **resolved as a chart obstruction**, §5 |
| horizon persistence | **yes**, on a seeded background, exactly |
| perturbation spectrum | **no** |
| retarded outer→inner transfer function | **no** |

Two horizon-penetrating time-domain constructions were built for a test scalar on
a fixed Tangherlini background: a Kerr–Schild slicing `t̃ = v − r` of the same
chart, and a tortoise `(t, r*)` evolution using the derived potential. Both are
stable and **both converge** — the Kerr–Schild frequency is flat to four digits
under a four-fold refinement and across three extraction windows. They do not
agree with each other:

| | `ℓ = 1` |
|--|--|
| Kerr–Schild | `1.01622 − 0.36231i` |
| tortoise `(t, r*)` | `1.01876 − 0.26404i` |
| real parts | agree to `0.3%` |
| damping rates | **apart by `37%`** |

A frequency-domain shooting solve did not converge on the damping rates either —
the centrifugal tail makes `u → const` a poor outgoing condition at finite
radius.

**So no quasinormal frequency is reported**, and the transfer function is not
built, because it is a ratio of the same two signals and would inherit the same
unresolved error. Two converged numbers that disagree mean a systematic error in
at least one construction, and the arc has a standing entry for exactly this:
**a converged number is not a correct number.**

*First thing to chase:* the Kerr–Schild operator does not converge to the exact
flat-space mode at its inner cut — error flat at `1.07e-02` across a four-fold
refinement — which points at the excision boundary rather than at the operator.

## Scope

- **Classical, spherically symmetric, one massless scalar.** No matter model from
  the rest of the arc appears; no charge, no winding, no `S³` harmonics above the
  monopole in the nonlinear sector.
- **Second-order accurate**, stated as such and measured, not asserted.
- **Horizon persistence is shown only on a seeded background**, where it is
  exact. A dynamically formed horizon is not evolved, for the reason §5 gives.
- **The spectrum and the transfer function are not delivered.** They are named,
  the numbers that say so are recorded, and nothing is quoted from them.

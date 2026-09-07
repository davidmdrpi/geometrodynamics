# Round 10: what data does joint closure retain?

Pre-registered in [`joint_closure_composition_prereg.md`](joint_closure_composition_prereg.md)
at [`b78157a`](https://github.com/davidmdrpi/geometrodynamics/commit/b78157a)
(amendment A1), committed before any of the code or numbers below. Baseline
main is `22f77a3`, which includes #286 and #287.

**Result.** For two independently prepared triangles the joint reference
measure exists, is exactly a product, and is a function of `|D1 D2|` alone.
No module in the inspected machinery supplies any *other* joint weight rule
for two disconnected preparations. Independent composition therefore selects
no counting function, and the conditional obstruction of the freeze's P3 is
**not** invoked: its hypotheses are not established here.

| verdict field | value |
|---|---|
| reference composition | `INDEPENDENT_PHASE_PRODUCT_VERIFIED` |
| reduction of specified rules | `PRODUCT_STATISTIC_SUFFICIENT` (reference rule, through the **absolute** product) |
| additional physical rule | `JOINT_WEIGHT_RULE_UNSPECIFIED` |
| consequence for selection | `NO_SELECTION_FROM_INDEPENDENCE` |

## The inherited conditioning is chosen, not derived

Everything below rests on round 8's conditioning: product Haar with a separate
window in each triangle's **phase**. Round 8 established that this is
justified by the repository's phase axiom (`history/closure.py:11`) and is
**not** forced by the zero set — conditioning on the numerators `N_i` instead
has the same support and gives the uniform product measure. The probe
recomputes the `N` control (`|grad N| = |q|`, constant on the circle) so the
nonuniqueness stays visible. Nothing here promotes the phase choice to a
geometric necessity.

## The joint geometry

With `u_i = s_{Ai} a_i`, `w_i = -s_{Bi} b_i`, `q_i = u_i x w_i`,
`t_i = 1 + u_i.w_i` and `s_i = u_i + w_i`:

    N_i = x_i . q_i,   D_i = t_i + x_i . s_i,   theta_i = atan2(N_i, D_i).

`s_i` is orthogonal to `q_i` and `|s_i| = sqrt(2 t_i)`, so on the closure
circle `D_i = t_i + sqrt(2 t_i) cos psi_i` in arclength. The joint closure
locus is `Gamma_1 x Gamma_2`, the two differentials occupy orthogonal tangent
blocks, and for `F = (theta_1, theta_2)`

    sqrt(det(dF dF^T)) = |grad theta_1| |grad theta_2| = |q_1||q_2| / |D_1 D_2|,

so the coarea density is `|D_1 D_2| / (|q_1||q_2|)` per product arclength.
The measured Gram matrix is exactly block diagonal (off-diagonal `0.0`), its
finite-difference determinant matches the closed form to `2.5e-08` (worst case,
all eight configurations) at the finest frozen step with clean second-order convergence, and an independent
analytic route — `grad theta = (D grad N - N grad D)/(N^2 + D^2)`, which never
uses `N = 0` — agrees to `4.8e-16`.

Sector masses follow: with the prior `1/16` over all 16 joint sectors, a
sector mass is `W(t_1) W(t_2) / (16 |q_1||q_2|)` where
`W(t) = int |t + sqrt(2t) cos psi| dpsi`. That has a closed form,

    W(t) = t (4 psi_0 + 4 k sin psi_0 - 2 pi),  k = sqrt(2/t),
    psi_0 = arccos(-1/k)      for t < 2;        W(t) = 2 pi t   for t >= 2,

continuous at `t = 2` and equal to `pi + 4` at `t = 1`. It matches split
quadrature to machine precision and uniform `2048`-point grids to `1.5e-08`
in normalised sector probabilities.

### The punctures are exactly the two excluded geodesic legs

`D_i` vanishes on `Gamma_i` precisely at `x_i = -u_i` and `x_i = -w_i`, the
points the freeze excludes. There

    |dD/dpsi| = sqrt(2t) |sin psi_0| = sqrt(t(2-t)) = |q|

**identically**, so the arc with `|D|/|q| < eta` is exactly `|psi - psi_0| < eta`
and its coarea mass is exactly `eta^2` per puncture, independent of the
settings. Measured against that law the relative error is `3.3e-05` at
`eta = 0.02` and falls as `eta`; the excluded joint fraction, by
inclusion–exclusion over the two factors, stays below `2.5e-04`.

### Finite windows reach the same limit

Windows are evaluated directly, not by inserting the limiting density. In
`x = sqrt(1-z^2) r(psi) + z qhat` with area element `dpsi dz`, `N = z|q|` and
`dist(theta, pi Z) < epsilon` is exactly `|N| < |D| tan epsilon`; the `z`
boundary is root-solved per `psi`. Normalised sector probabilities approach
the coarea values quadratically in `epsilon` — `1.2e-04, 3.0e-05, 7.6e-06,
1.9e-06` at `epsilon = 0.04, 0.02, 0.01, 0.005` for the `(1.0, 1.3)` pair —
including the asymmetric width pairs `(e, 2e)` and `(2e, e)`. The worst
final-window discrepancy over all eight configurations is `3.3e-05`, against
the frozen gate of `2e-3`.

## What the factorization does and does not establish

The reference density is a function of `|D_1 D_2|`, hence also of the signed
product. That is a design consequence recorded in the freeze, not a
discovery. Three exact controls at `t_1 = t_2 = 1` show what it is worth:

| control | statistic match | reference density gap | cubic weight gap |
|---|---:|---:|---:|
| reflection `psi -> -psi` (pair statistic) | `2.2e-16` | `2.2e-16` | `6.7e-16` |
| `(1,2)` vs `(sqrt2, sqrt2)` (equal signed product) | `1.1e-15` | `1.1e-15` | `0.137258300203` |
| `(1/4, 1)` vs `(-1/4, 1)` (equal absolute product) | `5.3e-16` | `5.3e-16` | `0.005` |

All six configurations are regular and inside the attainable range
`[1 - sqrt2, 1 + sqrt2]`. The reference density is constant on every one of
these level sets. The factorwise cubic `Phi(D) = D^2 (1 - D/5)` is not:
`f(1)f(4) = 4 != 9 = f(2)f(2)`, and `Phi(d) - Phi(-d) = -2 d^3 / 5` makes it
sign-sensitive. **Sufficiency of the reference measure is not a sufficiency
theorem for any other weight.**

## Q2: no inherited rule composes two disconnected preparations

Two findings, both demonstrated rather than asserted.

**The generic closure checker is rank one on a union.**
`history/closure.py` sums every event and worldline phase into a single
`total_phase` and accepts the history when that one number is within tolerance
of `pi Z`. Handed the union of two disconnected preparations that is the
condition `theta_1 + theta_2 in pi Z`, which the freeze forbids as a
substitute for the two independent conditions. It is not hypothetical: two
sub-loops at `+pi/2` and `-pi/2` are each at the maximum possible distance
from closure, yet their union is accepted with mismatch `0.0`. Separately, at
the module's own default `sigma = 0.6` the worst attainable mismatch `pi/2`
still scores weight `0.0325 > 0.01`, so the phase gate cannot reject **any**
history on phase alone.

**The based-loop composition theorem needs a common base point.**
`history_action.py` proves `theta[g1 . g2] = theta[g1] + theta[g2]` for loops
based at the same `x`, where the holonomies lie in the `U(1)` generated by `x`
and commute (residual `3.3e-16`). Two independent triangles have distinct base
points on distinct spheres; their holonomies are generic `SU(2)` elements with
measured non-commutativity `1.94`. The theorem does not reach them.

| module | applies to disconnected pairs | supplies a joint weight rule |
|---|---|---|
| `history/closure.py` | yes | no — rank-one summed rule |
| `bulk/history_action.py` | no | no — common base point required |
| `bulk/closure_current.py` | no | no — single-pair measures only |
| `transaction/network.py` | no | no — presupposes a throat connection |
| `transaction/derived_network.py` | no | no — same connected loop |

This is a search over the five modules the freeze named, at the pinned
baseline. **Absence here is a repository gap, not a theorem that no BAM
completion can supply such a rule.**

## What is not established

- No `Phi` is selected. The freeze's P3 obstruction is **not** invoked: its
  hypotheses — physical completeness of the allowed rule family, the
  composite-scalar identification, and the domain hypothesis for the
  multiplicative classification — are none of them established here.
- The factorwise cubic construction of P1 remains an admissible extension. It
  is not a derived physical preparation.
- No Born rule, no operational source-local readout, no Hilbert tensor
  product. The existing [tensor-product construction](tensor_product_emergence.md)
  still imports quantization and cannot certify a classical composition rule.
- Five quantities are guaranteed by the product construction and are recorded
  as regression oracles, never as selection evidence: product marginals,
  three-factor associativity, joint window factorization, the vanishing Gram
  off-diagonal (`theta_i` depends only on `x_i` — that is what "independently
  prepared" means) and copy exchange (the density expression is symmetric).
  Each level-set control additionally reports its product-space chord
  separation, so a control cannot pass by comparing a configuration with
  itself; the three separations are `2.49`, `0.567` and `0.519`.

## Reproduce

```bash
python -m experiments.closure_ledger.joint_closure_probe   # 15 required checks
python -m pytest -q tests/test_joint_closure.py             # 36 tests
```

Archived report:
[`runs/20260907_joint_closure_probe`](../experiments/closure_ledger/runs/20260907_joint_closure_probe/probe.md).

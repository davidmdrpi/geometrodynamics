# Round 10: joint closure and sufficient data

Two independently prepared triangles under the inherited chosen phase conditioning. Not a Born rule, not an operational readout, not a Hilbert tensor product.

Public freeze: `b78157a`. Seed `2026090710`.

Conditioning inherited from round 8 is **chosen**, not derived.

| verdict field | value |
|---|---|
| reference composition | **INDEPENDENT_PHASE_PRODUCT_VERIFIED** |
| reduction of specified rules | **PRODUCT_STATISTIC_SUFFICIENT** |
| additional physical rule | **JOINT_WEIGHT_RULE_UNSPECIFIED** |
| consequence for selection | **NO_SELECTION_FROM_INDEPENDENCE** |

Reduction scope: the reference rule only, through the ABSOLUTE product; the separations in the level-set controls show the factorisation does not transfer to the cubic weight.

## Level-set controls at `t1 = t2 = 1`

| control | statistic match | reference density gap | cubic weight gap |
|---|---:|---:|---:|
| reflection `psi -> -psi` (pair) | 2.220e-16 | 2.220e-16 | 6.66133814775e-16 |
| `(1,2)` vs `(sqrt2,sqrt2)` (product) | 1.110e-15 | 1.110e-15 | 0.137258300203 |
| `(1/4,1)` vs `(-1/4,1)` (absolute product) | 5.274e-16 | 5.274e-16 | 0.005 |

## Repository audit for a joint weight rule

| module | applies to disconnected pairs | supplies a joint weight rule |
|---|---|---|
| `geometrodynamics/history/closure.py` | True | False |
| `geometrodynamics/bulk/history_action.py` | True | False |
| `geometrodynamics/bulk/closure_current.py` | False | False |
| `geometrodynamics/transaction/network.py` | False | False |
| `geometrodynamics/transaction/derived_network.py` | False | False |

Search scope: the five modules named in the freeze's Q2, at the pinned baseline; not an exhaustive search of BAM.

| criterion | pass |
|---|---|
| Q1 finite-difference Gram matches the closed Jacobian | True |
| Q1 independent analytic route agrees to 1e-10 | True |
| Q1 density is invariant under a common SO(3) frame change | True |
| Q3 level sets contain genuinely distinct histories | True |
| Q1 closed form matches split quadrature | True |
| Q1 uniform sector-probability grids agree at 2048 | True |
| Q1 finite windows converge to the coarea limit | True |
| Q1 punctures are exactly -u and -w with |dD/dpsi| = |q| | True |
| Q1 excised mass follows the two-term law with an eta^6 residual | True |
| Q1 the accepted window set is connected on the registered grid | True |
| Q1 joint excluded fraction is bounded and vanishing | True |
| Q3 reflection preserves the pair statistic and the density | True |
| Q3 same signed product separates the cubic weight | True |
| Q3 same absolute product, opposite sign separates the cubic | True |
| Q2 the generic closure rule is rank one on a union | True |
| Q2 on-closure holonomies are central, so composition is rank one | True |

Structural regressions (guaranteed by the product construction; not independent evidence):

- product marginals (structural): True
- three-factor associativity (structural): True
- joint window factorisation (structural): True
- joint Gram off-diagonal (structural): True
- copy exchange (structural): True

Passed 16/16 required checks.

No `Phi` is selected, no Born rule is derived and no operational source-local readout is constructed.

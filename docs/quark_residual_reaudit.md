# The quark residual sector, re-derived from the corrected operator

*The downstream re-derivation PR #271 deferred. Module
`geometrodynamics/qcd/residual_audit.py`, probe
`experiments/closure_ledger/quark_residual_reaudit_probe.py` (8/8), tests
`tests/test_quark_residual_reaudit.py`.*

---

## The result, first

PR #271 corrected the 5D radial scalar operator — `V_tangherlini` had been
short of the minimally coupled massless scalar by `3A²/(4r²)` — and reopened
the **lepton** sector's geometric closure. `γ = 22.5` survived as a requirement
of the locked surrogate but lost its derivation from the canonical
`R_OUTER = 1.26` barrier geometry.

The frozen v3 **quark** lock derives three residual knobs from the *same*
eigensolver. The obvious expectation was that they would break the same way.

**They do not. All three move toward their locked values.**

| residual | locked | legacy | corrected | legacy off | corrected off |
|--|--|--|--|--|--|
| `pinhole` = `Σ V_max[1..5]` | 22.25 | 22.008240 | 22.331195 | `−1.09%` | **`+0.36%`** |
| `transport` = `mean ⟨u_{ℓ₁}\|ΔV\|u_{ℓ₂}⟩` | 0.54 | 0.544746 | 0.543800 | `+0.88%` | **`+0.70%`** |
| `resistance` = `transport · ln(α₅/α₁)` | 0.14 | 0.140683 | 0.139965 | `+0.49%` | **`−0.02%`** |

The pinhole residual **changes sign**: the legacy sum undershot the fitted
22.25, the corrected sum overshoots it.

> Two numbers in the previous publication of this table need correcting. The
> transport residual was `+0.87%`, which reproduces as `+0.88%`. The resistance
> residual was `−0.43%`, which is what you get using the **locked** transport
> `0.54` in the formula instead of the derived one; with the derived transport,
> as the formula is written, the legacy value is `+0.49%`. Both conventions are
> carried in the module output.

---

## Why the two sectors part company: elasticity

One barrier feeds both sectors, and the correction moves it the same way for
both — `Σ V_max[1..5]` at `R = 1.26` goes `22.008 → 22.331`, a `+1.47%` shift,
and that single number is *simultaneously* the lepton `γ` and the quark
`pinhole`. What differs is how hard each Hamiltonian leans on it.

| sector | corrected residual | elasticity | linearised error | **measured error** |
|--|--|--|--|--|
| lepton (`m_μ`) | `−0.75%` vs 22.5 | `−17.471` | `+14.01%` | **`+15.16%`** |
| quark (`m_s`) | `+0.36%` vs 22.25 | `+4.798` | `+1.77%` | **`+1.79%`** |

A factor of `3.64` in stiffness. The lepton chain converts a sub-percent
geometric residual into a fifteen-percent mass error; the quark ladder converts
a comparable residual into under two percent.

> **A percentage agreement between a geometric quantity and a fitted knob means
> nothing until multiplied by the elasticity of what it feeds.**

That is the transferable lesson of this round, and it is why #271's headline
and this one point in opposite directions without contradicting each other.
Both sectors' geometric agreements were reported for years as sub-percent
near-equalities, which reads like a rounding detail. In one sector it was; in
the other it never was.

### A number this round corrects

The README row written in PR #271 read:

> *the residual improves to −0.75 %, but `d ln m_μ/d ln γ = −17.5` at the lock
> makes that a **9 %** muon error*

The elasticity is right and the residual is right; the `9 %` is not. It came
from a generic illustration in the #271 probe docstring — *"a half-percent
error in gamma is a nine-percent error in the muon"* — applied to a residual
that is not a half percent. Linearising `−0.75%` gives `+14.0%`, and the locked
block actually returns `+15.2%`.

#271's own A/B/C table carried the correct value the whole time, in the row
`A corrected R=1.26, gamma[1..5]`. The prose and the table disagreed and the
table was right. Fixed at both README sites.

---

## The composite is not the sum of its parts

Every substitution below goes into the **locked** quark Hamiltonian with
**nothing retuned**, `d`-anchored, max relative error over `{s, c, b, t}`.

| substitution | legacy-derived | corrected-derived |
|--|--|--|
| `pinhole` alone | `3.61%` | `3.40%` |
| `transport` alone | `2.10%` | `2.00%` |
| `resistance` alone | `1.70%` | `1.61%` |
| **all three** | `3.44%` | **`3.78%`** |
| *the fitted v3 lock* | `1.61%` | `1.61%` |

**The ordering reverses.** Each corrected residual is individually closer to
its lock, yet the corrected triple scores *worse* than the legacy triple. The
legacy set enjoys a partial cancellation: its pinhole error alone costs `3.61%`
and the other two walk it back to `3.44%`. The corrected set has no such luck —
`3.40%` alone becomes `3.78%` together.

Neither composite comes near the fitted lock's `1.61%`.

> **"Each knob is derived to within 1 %" and "the derived knobs reproduce the
> ladder" are different claims, and only the first was ever established** — under
> either operator. The correction did not create this gap. It exposed it.

This is the round's main caution and the sharpest successor target. A residual
set that is individually right and jointly wrong points at a *missing
correlation* between the three quantities, not at three independent near-misses.

---

## The elasticity, and the pinhole the ladder actually wants

| species | local `d ln m/d ln pinhole` at the lock | secant across the two derived pinholes |
|--|--|--|
| `s` | `+4.7976` | `+4.8188` |
| `c` | `−1.5064` | `−1.6322` |
| `b` | `−1.8610` | `−1.9865` |
| `t` | `−2.0130` | `−2.1380` |

*(Reported as two quantities, not one — the distinction PR #271's review
insisted on. The local derivative is the headline; the secant is the honest
description of the interval actually traversed.)*

Scanning the locked ladder over the pinhole **alone**, it is minimised at
`22.228022` (`1.348%`), not at the fitted `22.25` (`1.610%`) — the knob was
fitted jointly, so a single-axis scan lands elsewhere. Measured against the
ladder's *own* optimum rather than against a fitted knob:

| pinhole | distance from what the ladder wants |
|--|--|
| v3 lock `22.25` | `+0.099%` |
| legacy `22.008` | `−0.989%` |
| corrected `22.331` | **`+0.464%`** |

The correction roughly halves the miss.

---

## The one thing that genuinely improves: cross-sector `R_OUTER`

Each sector, read as a demand on `R_OUTER` through its own locked scalar:

| operator | lepton (`γ=22.5`) | quark (`pinhole=22.25`) | split | brackets `1.26`? |
|--|--|--|--|--|
| legacy | `1.28737` (`+2.17%`) | `1.27229` (`+0.98%`) | `1.179%` | no |
| corrected | `1.26788` (`+0.63%`) | `1.25645` (`−0.28%`) | `0.906%` | **yes** |

Under the legacy operator both sectors demand **more** than `1.26`, putting the
canonical value outside the bracket they define — 0.81 bracket-widths below it.
Under the corrected operator they **straddle** it: the quark sector wants less,
the lepton sector wants more, and `1.26` sits at 31% of the way across.

This is evidence *for* `R_OUTER = 1.26`, and it is **not** the evidence #271
reopened — that was a single-sector fixed point, this is a two-sector bracket.

It is also weak, and is recorded that way. Two numbers 0.91% apart bracket
anything between them; nothing here derives `1.26` rather than admitting it.
**Suggestive, not a restored derivation.** The single-sector fixed point stays
reopened.

---

## The last lepton claim that reads a radial eigenvalue

> README: *"Resistance = 7π/100 — selected over `4·(ω−1)` by `R_OUTER`
> bisection."*

| candidate | value | off the locked `0.217869` |
|--|--|--|
| `7π/100` (operator-independent) | `0.219911` | `+0.937%` |
| `4(ω−1)`, legacy | `0.218908` | `+0.477%` |
| `4(ω−1)`, corrected | `0.223306` | `+2.495%` |

Under the legacy operator the **rejected** candidate fitted *better* than the
selected one. The choice therefore rested entirely on the `R_OUTER` bisection —
which #271 reopened. Under the corrected operator the competitor degrades to
`+2.50%` and `7π/100` wins on proximity as well.

> **Conclusion survives; its stated reason does not.**

The downstream `ε = resistance/k₅⁴ = 7π/(100·5⁴)` bridge is unaffected: it
consumes the closed form, not `ω`.

---

## `N` still drifts

| residuals pinned to | `N` range | width |
|--|--|--|
| baseline (all free) | `[430, 494]` | `64` |
| legacy-derived | `[434, 480]` | `46` |
| corrected-derived | `[438, 496]` | `58` |

Pinning the residuals to *corrected* geometry does not stabilise `N` any more
than pinning them to legacy geometry did. The v3 lock's own conclusion is
unchanged: **`β` remains the model's last fit knob**, and `N = 466` is a
compensator, not a topological invariant. Only the digits move — including `N`
at PDG masses itself.

---

## The ledger

| claim | verdict | evidence |
|--|--|--|
| quark shell-index axioms (ε, η, χ, phase) in `k₅ = 5` | **EXACTLY INVARIANT** | expressible in `k₅` alone; reads no radial operator |
| transport operator `V_{ℓ₂} − V_{ℓ₁}` as an algebraic object | **EXACTLY INVARIANT** | the correction `3A²/4r²` carries no `ℓ`; verified below `1e-12` |
| `pinhole = Σ V_max[1..5]` vs the fitted 22.25 | **NUMERICALLY SHIFTED, AND IMPROVED** | `−1.09% → +0.36%` |
| `transport` = mean off-diagonal vs the fitted 0.54 | **NUMERICALLY SHIFTED, AND IMPROVED** | `+0.88% → +0.70%` |
| `resistance = transport·ln(α₅/α₁)` vs the fitted 0.14 | **NUMERICALLY SHIFTED, AND IMPROVED** | `+0.49% → −0.02%` |
| the derived residuals **reproduce** the quark ladder | **INTERPRETATION CHANGED** | never established under either operator: `3.44%` / `3.78%` against the lock's `1.61%` |
| `N = 466` as a stable integer of the residual-pinned fit | **NUMERICALLY SHIFTED** | still a compensator; drift width stays O(50) |
| lepton `resistance = 7π/100` selected over `4(ω−1)` | **INTERPRETATION CHANGED** | conclusion survives on proximity; the selector that chose it is reopened |
| `R_OUTER = 1.26` from a single-sector fixed point | **INTERPRETATION CHANGED** | reopened by PR #271 and **not restored here** |
| `R_OUTER = 1.26` bracketed by two independent sectors | **NEW, AND WEAK** | corrected operator straddles it at 0.31 across a 0.91% window |

---

## Not re-run, and why

The same exclusion PR #271 made, for the same reason. The four quark
shell-index axioms are expressible in `k₅ = 5` alone; the Hopf transport phase
`φ_h = π/k₅` is topological; the v4 flavor-CP layer inherits the v3 eigenvalues
by construction. None of them reads the radial scalar operator.

**Proximity is not dependence.**

## Still open

Nothing here *derives* `22.25`, `0.54` or `0.14`. Each remains a fitted knob
that a geometric quantity happens to land near, and the round has changed how
near, not what kind of claim it is.

The composite gap is the sharp successor: three residuals individually within
`0.7%` of their locks that jointly miss the ladder by `3.8%` is not three
independent near-misses. It is a missing correlation between them — most likely
in how the three enter the Hamiltonian, since `transport` appears both on its
own and inside `resistance`. Finding or refuting that correlation is the work
this round hands forward.

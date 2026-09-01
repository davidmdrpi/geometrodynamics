# Does Gauss–Bonnet reopen the throat?

**5/5 checks pass — the branch is narrowed, not closed.**

Frozen in `docs/gauss_bonnet_prereg.md` before this module existed.

> ## The verdict

> **Gauss-Bonnet is NARROWED, not closed: positive/string-sign alpha reinforces the violation, and a global matter-NEC solution needs a Gauss-Bonnet length of at least half the BULK radius.**

> Raychaudhuri does not care about the gravitational action, so flare-out still forces `R_kk < 0`; only *which tensor supplies it* changes. The branch's hope was that `α_GB H_kk` supplies the negative part geometrically. Instead `H_kk` has the **same sign** as `R_kk` at every neck, so the best-motivated coupling makes the violation worse.

> **An earlier draft of this probe said the branch was closed. It is not.** A sufficiently negative constant coupling *does* satisfy the matter NEC — see B3, which replaces a withdrawn claim that the derivative expansion had broken down. What is closed is the cheap version.

---

## B0 — validate the implementation where the answer is known

In `D = 4` the Gauss–Bonnet invariant is topological, so `H_ab` must vanish identically:

| metric | `R_kk` | `H_kk` |
|--|--|--|
| 4D Schwarzschild (vacuum) | `+5.81e-09` | `-9.07e-09` |
| 4D general A(r) = 1 - 2/r + 0.3/r^2 (non-vacuum) | `+3.15e-09` | `-4.66e-09` |
| 4D general A(r) = 1 + 0.2 r^2 (non-vacuum) | `-4.68e-10` | `-7.38e-10` |

> In D = 4 the Gauss-Bonnet invariant is a total derivative, so its variation H_ab vanishes identically. Getting zero for a general NON-VACUUM A(r) as well as for Schwarzschild shows the zero is a property of the implementation and the dimension, not of R_ab = 0.

## B1 — the decisive sign

```
R_kk = −3(N f″ − N′f′)/(N f)
H_kk = 12(f′² − 1)(N f″ − N′f′)/(N f³)
```

The lapse drops out at the neck for the same reason it did in P4: `N′` multiplies `f′`, and `f′(0) = 0` is what *makes* it a neck. So for any `N`:

```
R_kk = -393.343998 = −3f″/f₀
H_kk = -206292.667499 = −12f″/f₀³
H_kk/R_kk = 524.458664 = 4/f₀²  > 0
```

Against an independent Riemann tensor built by numerical differentiation, sharing no algebra with those closed forms:

| `s` | `R_kk` closed | `R_kk` numeric | `H_kk` closed | `H_kk` numeric |
|--|--|--|--|--|
| `0.000` | `-393.3440` | `-393.2924` | `-206292.67` | `-206170.92` |
| `0.085` | `-104.4537` | `-104.4674` | `-14547.43` | `-14555.99` |
| `0.198` | `-10.4993` | `-10.4995` | `-146.98` | `-146.99` |

Worst relative error `5.9e-04` (finite-difference limited).

> **The branch was invoked so that alpha_GB H_kk could supply the negative part geometrically and leave T_kk >= 0. Instead H_kk has the SAME sign as R_kk at every neck, so for the heterotic sign alpha_GB > 0 the Gauss-Bonnet term makes the violation WORSE by the factor (1 + 4 alpha_GB/f_0^2).**

> H_kk/R_kk = 4(1-f'^2)/f^2 = 4 mu/f^4 with mu = f^2(1-f'^2), the same quantity finite_mouth's P2 showed is continuous across the seam. Since mu = b^2 > 0 throughout, the reinforcement holds along the whole throat and not only at the neck.

## B2 — the coupling the matter NEC would demand

| case | `α_GB` | `8πG₅T_kk` at the neck | NEC? |
|--|--|--|--|
| heterotic sign, alpha = +b^2/4 | `+0.001906728` | `-786.6880` | **no** |
| alpha = 0 (pure Einstein) | `+0.000000000` | `-393.3440` | **no** |
| alpha = -b^2/8 (half-way) | `-0.000953364` | `-196.6720` | **no** |
| threshold alpha = -b^2/4 | `-0.001906728` | `+0.0000` | **yes** |
| beyond, alpha = -b^2/2 | `-0.003813456` | `+393.3440` | **yes** |

`alpha_GB <= -f_0^2/4 = -R^2 sin^4 a / 4`, i.e. `α_GB ≤ -0.001906728` here.

> **Heterotic string theory gives alpha_GB = alpha'/8 > 0. That is exactly the sign that deepens the violation here, so the best-motivated value of the coupling is the one that fails hardest.**

## B3 — a neck-only cancellation is not a wormhole

**Withdrawn first.** An earlier draft read `α_GB H_kk/R_kk = −1` as the derivative expansion breaking down. In `D = 5` the Lovelock tower already terminates at Gauss–Bonnet:

| Lovelock `k` | status in `D = 5` |
|--|--|
| 1 | dynamical |
| 2 | dynamical |
| 3 | identically zero |
| 4 | identically zero |

> **An earlier draft read alpha_GB H_kk/R_kk = -1 as 'the derivative expansion has broken down and the whole Lovelock tower contributes'. Both halves are wrong in D = 5: cubic Lovelock antisymmetrises six indices and is IDENTICALLY ZERO here, so there is no further tower; and exact Einstein-Gauss-Bonnet is a complete classical theory with second-order equations, not a truncation that can lose validity. Order-one is just order-one.**

**What replaces it.** The NEC must hold *along* the throat, and the neck is its easiest point:

| requirement | `α_GB` | `min T_kk` over throat | NEC everywhere? |
|--|--|--|--|
| neck only, `−b²/4` | `-0.001906728` | `-98.33` | **no** |
| whole throat, `−R²/4` | `-0.250000000` | `-0.000000` | **yes** |

The mouth sets `f_m⁴/(4b²) = R²/4` **exactly**, independent of `a`, so the global requirement is `131×` the neck-only one here — and the Gauss–Bonnet length must be `√|α| ≥ 0.50 R`, half the radius of the closed universe rather than a short-distance scale.

> The NEC has to hold along the throat, not only at the neck. Pointwise it demands alpha <= -f^4/(4 mu), which is WEAKEST at the neck; at the mouth f_m^4/(4b^2) = R^2/4 exactly, independent of the mouth angle. So a global solution needs alpha <= -R^2/4, i.e. a Gauss-Bonnet length of at least R/2 -- half the radius of the closed universe, not a short-distance scale.

> **NARROWED, NOT CLOSED. Positive/string-sign alpha reinforces the violation. A sufficiently negative constant coupling can satisfy the matter NEC, but only at a Gauss-Bonnet length comparable to the whole closed universe. Global existence and stability of such a solution are untouched here.**

## B4 — the ledger

| claim | verdict | evidence |
|--|--|--|
| the Lanczos implementation is trustworthy | **VALIDATED IN D = 4** | H_ab vanishes identically for Schwarzschild AND for two general non-vacuum A(r), so the zero is not an artefact of R_ab = 0 |
| Gauss-Bonnet can supply the negative null stress | **NO -- IT HAS THE SAME SIGN** | H_kk/R_kk = 4 mu/f^4 = 524.4587 > 0 at the neck, positive along the whole throat since mu = b^2 > 0 |
| the lapse can help in Gauss-Bonnet | **NO -- IT DROPS OUT AT THE NECK** | N' multiplies f', and f'(0) = 0 is what makes s = 0 a neck; the same structure as P4 |
| the heterotic coupling alpha' /8 > 0 helps | **NO -- IT MAKES IT WORSE** | the best-motivated sign deepens the violation by (1 + 4 alpha_GB/f_0^2); the matter NEC would need alpha_GB <= -0.001906728 |
| alpha H_kk/R_kk = -1 means the expansion has broken down | **WITHDRAWN -- WRONG IN D = 5** | cubic Lovelock antisymmetrises six indices and is IDENTICALLY ZERO in five dimensions, so there is no further tower; and exact EGB is a complete classical theory, not a truncation that can lose validity |
| a negative coupling that works at the neck suffices | **NO -- THE NECK IS THE EASIEST POINT** | the NEC needs alpha <= -f^4/(4 mu) pointwise, weakest at the neck; alpha = -b^2/4 leaves T_kk as low as -98.3 elsewhere on the throat |
| a global matter-NEC solution is available cheaply | **ONLY AT A COSMOLOGICAL GAUSS-BONNET LENGTH** | the mouth sets f_m^4/(4b^2) = R^2/4 exactly, independent of the mouth angle, so alpha <= -R^2/4 and sqrt|alpha| >= R/2 -- 131x the neck-only requirement here |

**What this narrows.** The classical-geometry escape from the source audit. Gauss-Bonnet is the only branch keeping both a classical geometry and a traversable throat, and it is not closed: a sufficiently negative constant coupling does satisfy the matter NEC. What is closed is the cheap version -- the string-motivated sign fails, and a global solution needs the Gauss-Bonnet length to be half the radius of the closed universe rather than a short-distance scale.

**What remains.**

- **0 negative-coupling EGB** — not refuted -- alpha <= -R^2/4 satisfies the matter NEC along the throat, at a Gauss-Bonnet length of half the bulk radius; global existence and stability are open
- **1 accept the horizon** — the Tangherlini branch N(0) = 0 as the particle, abandoning MTY traversability
- **2 ghost** — a wrong-sign field, with its stability problem
- **3 quantum stress** — Casimir-type support, so the geometry is no longer classical
- **4 reinterpret** — particle exchange needs no traversable throat

**Not refuted here.** Negative-coupling constant EGB itself, which this round NARROWS rather than closes. Also dilatonic alpha(phi) L_GB, where the scalar's own stress enters and known 5D solutions exist, and f(R). The Lovelock tower is not a separate premise in D = 5: it terminates at Gauss-Bonnet.

# Does Gauss–Bonnet reopen the throat?

**5/5 checks pass — and the branch closes.**

Frozen in `docs/gauss_bonnet_prereg.md` before this module existed.

> ## The verdict

> **Gauss-Bonnet does NOT reopen the throat: it reinforces the Einstein violation.**

> Raychaudhuri does not care about the gravitational action, so flare-out still forces `R_kk < 0`; only *which tensor supplies it* changes. The branch's hope was that `α_GB H_kk` supplies the negative part geometrically. Instead `H_kk` has the **same sign** as `R_kk` at every neck, so the best-motivated coupling makes the violation worse.

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

## B3 — and even then, outside the theory's own regime

| `α_GB` | `α_GB H_kk / R_kk` |
|--|--|
| `-0.000019067` | `-0.0100` |
| `-0.000190673` | `-0.1000` |
| `-0.000953364` | `-0.5000` |
| `-0.001906728` | `-1.0000` |

The required Gauss–Bonnet length is `√|α| = 0.043666`, exactly `0.5000 × b` — tied to the throat radius rather than to a separate short scale.

> alpha_GB H_kk/R_kk = 4 alpha_GB/f_0^2, which is exactly -1 at the threshold: the 'correction' equals the term it corrects. A derivative expansion truncated at Gauss-Bonnet has stopped meaning anything there, and the whole Lovelock tower would contribute. The required Gauss-Bonnet length is f_0/2, i.e. tied to the throat radius rather than to a separate short scale.

## B4 — the ledger

| claim | verdict | evidence |
|--|--|--|
| the Lanczos implementation is trustworthy | **VALIDATED IN D = 4** | H_ab vanishes identically for Schwarzschild AND for two general non-vacuum A(r), so the zero is not an artefact of R_ab = 0 |
| Gauss-Bonnet can supply the negative null stress | **NO -- IT HAS THE SAME SIGN** | H_kk/R_kk = 4 mu/f^4 = 524.4587 > 0 at the neck, positive along the whole throat since mu = b^2 > 0 |
| the lapse can help in Gauss-Bonnet | **NO -- IT DROPS OUT AT THE NECK** | N' multiplies f', and f'(0) = 0 is what makes s = 0 a neck; the same structure as P4 |
| the heterotic coupling alpha' /8 > 0 helps | **NO -- IT MAKES IT WORSE** | the best-motivated sign deepens the violation by (1 + 4 alpha_GB/f_0^2); the matter NEC would need alpha_GB <= -0.001906728 |
| a negative coupling of the right size would work | **ONLY OUTSIDE THE THEORY'S OWN REGIME** | at threshold alpha_GB H_kk/R_kk = -1.0, so the correction equals the leading term and truncating Lovelock at Gauss-Bonnet is unjustified |

**What this closes.** The classical-geometry escape from the source audit. Gauss-Bonnet was the only branch keeping both a classical geometry and a traversable throat, and for the natural D = 5 higher-curvature term with its best-motivated sign it fails.

**What remains.**

- **1 accept the horizon** — the Tangherlini branch N(0) = 0 as the particle, abandoning MTY traversability
- **2 ghost** — a wrong-sign field, with its stability problem
- **3 quantum stress** — Casimir-type support, so the geometry is no longer classical
- **4 reinterpret** — particle exchange needs no traversable throat

**Not refuted here.** Dilatonic alpha(phi) L_GB, where the scalar's own stress enters and known 5D solutions exist; the full Lovelock tower; and f(R). Those are separate premises, and this round tested constant-coupling EGB only.

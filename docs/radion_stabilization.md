# Radion stabilization from the primordial EM-capped throat: V_eff(φ), α, and the radion mass (PR #226)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. PR #225 left exactly one dynamical
> problem: *deriving ρ\* from the EM cap's own equations*. This PR
> executes it: the radion effective potential of the primordial
> EM-capped throat, the α it selects, and the radion mass. The
> companion probe machine-checks every claim (~2 s).

## 0. The answer, stated first

**The Hopf charge stabilizes the radion.** The committed throat is a
*dyon* — the #55 cap puts a fixed electric flux through it (charge
without charge), the #58 topology puts one Hopf charge on it (Σc₁ = 0
pair creation) — and under the #225 dilaton coupling e^{−√3φ}F² its
two flux energies carry **opposite radion charges**:

```
U_el  ∝ e^{+√3 φ}     (fixed charge: weaker coupling → cheaper field)
U_mag ∝ e^{−√3 φ}     (fixed flux:   weaker coupling → dearer flux)
```

With the #222 primordial frame factor e^{aφ} (a = 1/(2√3)), the
effective potential

```
V_eff(φ) = U_el,0·e^{pφ} + U_mag,0·e^{rφ},    p = 7/(2√3),  r = −5/(2√3)
```

has one exponent of each sign: **a minimum exists**, in closed form,
with a coefficient-free mass identity. Without the Hopf charge every
committed term has a positive exponent and the radion runs away to
decompactification (α → 0): the **no-go is proven**, and the #58
topological charge is exactly what saves the vacuum. The landing:
**α\* = √(5κ/28)** — order one (0.4226) at the symmetric Dirac point;
the observed α requires the electric–magnetic cap asymmetry
**κ = 2.98×10⁻⁴**, the new, sharply-quantified open problem.

## 1. The committed cap, imported

The #55 equilibrium (re-read from its committed ledger at machine
zero — continuing the #225 α-dependent holdout):

- the throat is the inner boundary that **caps the Coulomb field**:
  U_EM = αℏc/(2R), finite, with U/mc² = α/2 exactly;
- E(R) = A/R + B·R² with **A = αℏc/2** (a genuine function of α),
  equilibrium R\* = (A/2B)^{1/3}, stable (E″ = 6B > 0);
- the committed anchor R\* = λ_C (the Compton anchor of the matter
  sector).

## 2. The radion charges of the two flux energies

The #225 dilaton coupling is a dielectric ε(φ) = e^{−√3φ}. Radial
quadrature with the throat cap (machine-checked):

- **fixed electric charge**: the Gauss law pins D = q/4πr²
  independently of φ, so U_el = ∫D²/2ε ∝ **e^{+√3φ}**;
- **fixed magnetic/Hopf flux**: topology pins F, so
  U_mag = ∫εB²/2 ∝ **e^{−√3φ}**;
- the Dirac partner ratio U_mag/U_el = 1/(4α²) exactly (g = 2π/e) —
  the definition of κ = 1;
- the **primordial radius** (#222 forced the throat primordial —
  fixed in 5D-frame lengths; the #225 Einstein frame converts
  L_E = e^{−aφ}L₅) adds e^{+aφ} to both.

Exponents: p = √3 + a = 7/(2√3) ≈ 2.02 and r = −√3 + a = −5/(2√3)
≈ −1.44 — **one of each sign**.

## 3. V_eff(φ): the minimum, the no-go, the mass identity

- **The dyonic minimum exists** — numeric minimization = closed form
  to 10⁻⁸ — with V″ > 0.
- **The two-exponential theorem** (machine-checked symbolically and
  numerically): for V = Pe^{pφ} + Qe^{rφ},

```
m_φ² = V″(φ*) = −p·r·V_min = (3 − a²)·V_min = (35/12)·V_min
```

  — independent of every coefficient.
- **The no-go**: with κ = 0 (no Hopf charge) and the 5D-frame
  cohesion (exponent +a), every term has a positive exponent — V is
  monotonic and the radion runs away to φ → −∞ (α → 0,
  decompactification). **The topological charge of #58 is the
  stabilizer.**

## 4. α at the minimum

```
α*² = (κ/4)·(√3−a)/(√3+a) = 5κ/28
```

- At the **symmetric Dirac point** (κ = 1, U_mag = U_el/4α² from
  g = 2π/e): **α\* = √(5/28) = 0.4226** — order one, the self-dual
  landing, 58× the observed coupling. The mechanism is derived and
  modulus-free (pure exponent-and-topology arithmetic, anchor-free).
- The **observed α** requires κ = 28α²/5 = **2.98×10⁻⁴** (verified by
  direct minimization). The #165 guardrail scan finds no
  closure-constant match for κ (nearest: α/k₅² at 2.1% — rejected):
  the smallness of the electric–magnetic cap asymmetry is real and is
  **the named open problem** — strictly sharper than #225's "derive
  ρ\*".
- **Cohesion robustness**: both derived provenance cases (5D-fixed
  tension → e^{+aφ}; Einstein-frame-fixed → e^{−2aφ}) tilt α\* by
  < 0.3% at the #55-scale cohesion; the fiber-wrapped case carries
  exactly the magnetic exponent and is absorbed into κ.

## 5. The radion mass

Closed chain at the minimum (machine-checked):

```
U_mag/U_el = (√3+a)/(√3−a) = 7/5,   E_min = (12/5)·U_el*,   E″(φ*) = 7·U_el*  exactly
```

Anchored (B4 one-anchor discipline — both committed anchors reported):

| anchor | E″(φ\*) per throat |
|---|---|
| #225 fiber (R_E = 2l_P/√α): E″ = (7/4)α\*^{3/2}m_P | **1.33×10¹⁶ GeV** at observed α (GUT-scale — the radion is heavy, no fifth-force tension); 0.48 m_P at the Dirac point |
| #55 Compton (R = λ_C): E″ = (7α/2)m_ec² | **13.1 keV** |

The canonical mass is m_φ² = 16πn·E″(φ\*)/m_P² with n the throat
number density — the one astrophysical input, stated not derived.

## 6. Arc consistency — nothing else moves

- The minimum condition fixes α(φ), i.e. ln ρ: its gradient in the
  #225 five-knob log space is **exactly the (1,1,0,0,0) cap row the
  rank audit reserved** — rank 2, ∇α annihilated to machine zero,
  the three α-decoupled flats unchanged.
- The #225 ledger re-reads pass (ρ\* = 2/√α; the model-unit
  conversion 23.41 l_P), and the α law closes the loop: the fiber
  sits at ρ = 3.08 (near-Planckian) at the Dirac point, 23.41 at the
  observed-α minimum.
- #222's ×210 anchor exclusion — the verdict that forced the throat
  primordial, supplying the e^{aφ} frame factor — is re-read from
  its committed ledger.

## 7. Honest scope

- **κ is not derived.** κ = 1 is the symmetric Dirac point; O(1)
  Wu–Yang and cap-geometry factors can shift it by O(1), not by
  3×10³. The observed α needs κ = 2.98×10⁻⁴ — quantified, open.
- V_eff is the Born–Oppenheimer energy of a single throat against
  the asymptotic radion; the canonical mass carries the throat
  density n.
- Tree-level throughout: no fiber Casimir energy, no loop-induced
  potential, no backreaction. #225's V_tree = 0 in vacuo is
  unchanged — V_eff here is *sourced by the throat*, which is the
  point.
- The quadratures use the flat capped-field model of #55; a full 5D
  field profile shifts the O(1) coefficients inside κ, not the
  exponents (the exponents are #225's machine-checked Weyl algebra).
- The two anchors (#55 Compton, #225 Planck-fiber) are reported
  separately; α\* itself is anchor-free.

## 8. What would falsify this

- A committed term with an exponent below r = −5/(2√3) — the minimum
  would move or vanish. (Checked: the derived cohesion cases lie in
  (−2a, +a), strictly inside.)
- The Hopf charge removed from the throat — but #58's Σc₁ = 0 is
  committed topology, and the no-go shows the vacuum decompactifies
  without it.
- A closure-constant match for κ — the asymmetry would be
  numerological, not dynamical. (Checked: none within 1%.)
- The stabilizer row landing off the reserved rank-audit direction —
  something else in the arc would have to retune. (Checked: exactly
  (1,1,0,0,0), ∇α annihilated.)
- The mass identity failing — the potential would not be the
  two-exponential form derived. (Checked: (35/12)V_min numerically
  and symbolically; E″ = 7U_el\*.)

## 9. Companion probe

`experiments/closure_ledger/radion_stabilization_probe.py` (T1–T9,
~2 s): the #55 ledger re-read; the dielectric quadratures; the
minimum, the no-go, and the symbolic mass identity; α\* closed form,
the Dirac point, κ_required and its guardrail; the anchored mass
scales; the rank-audit slot check.

**Verdict:**
`THE_HOPF_CHARGE_STABILIZES_THE_RADION_V_EFF_HAS_A_MINIMUM_WITH_ALPHA_STAR_EQUALS_SQRT_5KAPPA_OVER_28_ORDER_ONE_AT_THE_DIRAC_POINT_AND_M_PHI2_EQUALS_35_OVER_12_V_MIN_THE_OBSERVED_ALPHA_NEEDS_KAPPA_3E_MINUS_4_THE_NAMED_OPEN_PROBLEM`

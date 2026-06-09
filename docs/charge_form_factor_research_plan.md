# Finite-momentum charge form factor on the antipodal cavity (PR #146)

> **Framing (to avoid a category error).** QFT on the *fixed classical* throat
> geometry — geometry → fields, **not** quantum gravity. The form factor is
> that of the matter QFT's dressed charge on the classical antipodal cavity.

PR #145 computed the dressed charge at zero momentum transfer — `F(0) = c₁`
exactly, the Ward identity `Z₁ = Z₂` — and flagged the q ≠ 0 structure as its
lead open item. This PR supplies it: the charge form factor `F(q)`, the
finite-q Ward identity that protects its normalisation, and the charge radius
it defines.

## The form factor

The dressed particle's charge density `ρ(x)` on the cavity integrates to the
exact charge `c₁` (#145); the form factor is its Fourier transform about the
charge centroid,

    F(q) = ∫ ρ(x) e^{iq(x − x̄)} dx,    F(0) = c₁,

with the small-q expansion `F(q) ≈ c₁(1 − q² r_c²/2)` defining the charge
radius `r_c² = Var_ρ(x)`.

## The finite-q Ward identity: the Bethe sum rule

Current conservation at **every** momentum transfer is the Bethe sum rule

    Σ_m (E_m − E_n) |⟨m| e^{iqx} |n⟩|² = q²,

the finite-q generalization of the TRK sum rule that gave `Π(0) = 0` in #144
(TRK is its `q² → 0` coefficient): the double commutator
`[e^{−iqx},[H,e^{iqx}]] = 2q²` is V-independent. **Verified numerically to
~1e-4 across q ∈ [0.5, 10]** — the finite-q face of the #142 Ward identity.

## The dressed density and the exact charge

The one-loop dressed state `|k̃⟩ = √Z₂(|k;0⟩ + Σ a_nm|n;m⟩)` with cloud
amplitudes `a_nm = g_knm/(s₀ − s_nm)`:

- reproduces the #145 Dyson Z₂ **exactly** (`1/(1 + Σa²) = 1/(1 − Σ′)`,
  machine precision — the two one-loop pictures agree);
- its real-space charge density integrates to `c₁` **exactly** at every
  coupling — the #145 anchor, now in real space.

## The charge radius is geometric

`r_c = 0.2649` (tortoise units): the radius from the density variance equals
the radius from the small-q fall-off of `F(q)`, and the one-loop cloud moves
it only at the **1e-4** level (× coupling²). The charge radius is the **bare
cavity mode profile** — a finite, O(cavity-scale) geometric length. The throat
charge is not pointlike: it is spread over the cavity with **no UV
divergence** (the form-factor face of the #55 finite self-energy
`U_EM/(mc²) = α/2`), its size set by the classical geometry with the QFT
dressing a small correction on top — geometry → fields, the program's arrow.

## The counterfactual

A cloud carrying the wrong charge (`c′ ≠ c₁`) shifts the total dressed charge
away from `c₁` (computed: −0.003 to −0.014): the `F(0)` anchor — and the whole
form-factor normalisation — rests on exact charge conservation at the unitary
throat (`Σc₁ = 0`, #58/#141/#142/#145); the absorbing throat leaks charge
(#142) and loses it. The same single postulate as stable matter, gauge
invariance, the massless photon, and charge universality.

## Scope and epistemic ledger

- **Derived:** the Bethe sum rule at every q; the exact real-space charge;
  dressed-state Z₂ = Dyson Z₂; the geometric finite charge radius with a
  parametrically small cloud correction.
- **Modelled:** the cubic coupling (the #136 posture); the one-loop dressing
  truncation; the radial reduction (no relativistic recoil).
- **Input:** the boundary value `α(μ₀) ≈ 1/137` (#143).
- **Open:** relativistic recoil and the F₁/F₂ (electric/magnetic)
  decomposition — the magnetic form factor / g−2 connection is #62's
  territory; higher loops; the absolute normalisation (#133); the flavor
  residuals (#134).

## Reproduce

```bash
python -m experiments.closure_ledger.charge_form_factor_probe
# Verdict: CHARGE_FORM_FACTOR_FINITE_Q_WARD_ANCHORED_GEOMETRIC_CHARGE_RADIUS
```

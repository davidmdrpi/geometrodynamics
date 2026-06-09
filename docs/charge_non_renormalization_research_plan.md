# Charge non-renormalization: Z₁ = Z₂ on the antipodal cavity (PR #145)

> **Framing (to avoid a category error).** QFT on the *fixed classical* throat
> geometry — geometry → fields, **not** quantum gravity. The renormalization
> constants are those of the matter QFT on the classical antipodal cavity.

PR #144 completed the one-loop two-point sector: the matter self-energy Σ
(#136) carries the wavefunction renormalization **Z₂**, the photon vacuum
polarisation Π (#144) carries **Z₃**, and the gauge vertex (#141/#142) carries
**Z₁**. The renormalized charge is `e = (Z₂/Z₁)·√Z₃·e₀`, so the Ward identity
**Z₁ = Z₂** is the statement that *charge is not renormalized by the matter
sector* — only the universal photon factor `√Z₃` dresses the charge. This PR
computes Z₁ = Z₂ on the antipodal cavity, closing the renormalization-constant
triangle.

## Why this probe, now

With Π in hand (#144), all three renormalization constants are on the table for
the first time. The #142 Ward–Takahashi identity `Γ(p,p) = ∂S⁻¹/∂p` ties Z₁ to
Z₂ structurally; this probe realizes it numerically at one loop, the same
upgrade #144 performed for transversality.

## The charged model

A charged matter field χ (the odd-l Dirichlet tower, Hopf charge `c₁ = 1`,
#129/#58) coupled to a neutral field φ (the even-l Neumann tower) through the
cubic triple-overlap vertex `g·χ†χφ` (#136/#137), charge-conserving at the
vertex (the #141 structure). The external pole sits below the lowest decay
threshold (the #136 stable-particle kinematics), so all constants are real.

## The computation

- **Z₂ ≠ 1, two ways:** `Σ′(s₀)` from the analytic spectral sum and from a
  central finite difference of the #136 self-energy agree to ~1e-12;
  `Z₂ = 1/(1 − Σ′) ≈ 0.986` at `g = 1` — the dressed particle is genuinely
  renormalized.
- **The Ward identity, term by term:** the q = 0 photon insertion on the
  internal charged line doubles the propagator,
  `Λ(0) = Σ c₁|g|²/(s₀ − s_nm)² = −c₁·Σ′(s₀)` exactly (machine precision; the
  neutral line carries zero charge and contributes nothing). Equivalently
  **Z₁ = Z₂**.
- **The dressed charge is exact and universal:**
  `F(0) = Z₂(c₁ + Λ(0)) = c₁(1 − Σ′)/(1 − Σ′) = c₁` to machine precision —
  and across species `(l_χ, l_φ, g) ∈ {(1,0,1), (3,0,1), (1,2,0.7), (1,0,0.5)}`
  the Z₂ varies (0.9855–0.9963) while `F(0) − c₁ = 0` identically. Each
  sector's self-interaction cancels out of its own charge — why different
  generations (`k ∈ {1,3,5}`, #71) carry exactly the same unit charge — and
  charge renormalization collapses to `e = √Z₃·e₀`: **the running of α is
  purely the #144 vacuum polarisation.**
- **The counterfactual:** a charge-violating vertex (internal charge
  `c′ ≠ c₁`) breaks the cancellation — `F(0) − c₁` ranges from `−0.003`
  (c′ = 0.8) to `−0.014` (c′ = 0) and becomes coupling-dependent. The
  protection is exact charge conservation at the throat — `Σc₁ = 0` from the
  unitary antipodal mirror (#58/#141/#142); the absorbing throat leaks charge
  (#142) and loses it. The same single postulate as stable matter, gauge
  invariance, and the massless photon (#129/#130/#142/#144).

## Scope and epistemic ledger

- **Derived:** `Λ(0) = −c₁Σ′` term by term (Z₁ = Z₂, machine precision); the
  exact, universal dressed charge `F(0) = c₁`; `e = √Z₃·e₀`; the dependence on
  exact charge conservation at the unitary throat.
- **Modelled:** the cubic coupling magnitude (the #136 posture); the q = 0
  charge-operator insertion (no q ≠ 0 form factors).
- **Input:** the boundary value `α(μ₀) ≈ 1/137` (#143, the 137 problem).
- **Open:** q ≠ 0 form factors; higher loops; the absolute normalisation
  (#133); the flavor residuals (#134).

## Reproduce

```bash
python -m experiments.closure_ledger.charge_non_renormalization_probe
# Verdict: CHARGE_NON_RENORMALIZATION_Z1_EQ_Z2_DRESSED_CHARGE_EXACT_UNIVERSAL
```

# The dynamical absorber: S = S_field + S_absorber + S_coupling (PR #214)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. #213 derived the Feynman propagator
> from the frozen bulk using complete absorption as a *boundary
> condition*. This PR takes the decisive step: the absorber becomes an
> actual degree of freedom with its own action, absorption becomes an
> **outcome of the dynamics**, and the iε becomes a **computed damping
> rate**. The companion probe machine-checks every claim (~1 min).

## 0. What #213 imposed, and what this PR derives

#213's chain was: static bulk ⇒ time-symmetric field; closed bulk ⇒
every wavefront returns; *complete histories leave no free remnants* ⇒
G_F unique. The last step was a boundary condition — physically
motivated by the refocusing geometry, but still a condition *imposed on
solutions*, not a consequence of an equation of motion. The decisive
successor introduces the absorber itself:

```
S = S_field + S_absorber + S_coupling

S_field    = Σₙ ∫ (u̇ₙ² − ωₙ²uₙ²)/2 dt          (the conformal tower ωₙ = n/R)
S_absorber = Σⱼ ∫ (q̇ⱼ² − νⱼ²qⱼ²)/2 dt          (a bank of oscillators —
                                                the throat's internal continuum)
S_coupling = −Σⱼ gⱼ ∫ qⱼ Φ(ψₐ, t) dt ,          Φ(ψₐ) = Σₙ Yₙ(ψₐ) uₙ ,
                                                Yₙ(ψ) = sin(nψ)/(sinψ·√(2π²))
```

Everything is classical and Gaussian — which is exactly what makes the
elimination *exact*, not perturbative. The engine's `MouthState.modes`
(the mode banks every throat mouth already carries in
`geometrodynamics/transaction/`) are precisely this degree of freedom;
this PR gives it an action.

## 1. Exact elimination

Integrating out the absorber is Gaussian-exact: the field's effective
resolvent is

```
G_eff(Ω)⁻¹ = ω₀² − Ω² − Σ(Ω) ,      Σ(Ω) = Σⱼ gⱼ²/(νⱼ² − Ω²) ,
```

verified against direct inversion of the full coupled (N+1)×(N+1)
system at complex frequencies to < 10⁻¹⁰, with the coupled stiffness
positive-definite (the absorber is a *stable* physical system — no
runaway, no ghost). The continuum limit of the flat bank has the closed
form Σ_c(z) = (κ/2z)·ln[((hi−z)(lo+z))/((hi+z)(lo−z))], verified
against a 2×10⁵-oscillator sum.

## 2. The damping is derived, and it is the ε

The resolvent pole moves **below the real axis by a computed amount**:
Ω* = ω̃ − iγ/2. Three independent measurements agree at the design
point γ = 0.05:

| measurement | γ |
|---|---:|
| complex pole of the finite-bank resolvent | 0.0499 |
| live time-domain energy decay of the coupled system | 0.0499 |
| golden-rule formula γ = πg²ρ/(2ω₀²) | 0.0500 |

and the scaling law is exact: γ ∝ g² with log–log slope **1.9996**
across g = 0.05–0.4, the golden-rule ratio staying within 1% over the
whole grid. The ordering integrals of #213 now read

```
I₊(Δ) = −(1/2ω̃) · 1/(Δ − ω̃ + iγ/2)
```

— **the iε is filled in as a physical rate** (numerical quadrature of
the damped offer segment against the closed form: 3×10⁻⁸). In the
covariant normalization ε_derived = ω̃γ, and the Δ²-pole form agrees
with the true pole to O(γ²) (gap 3×10⁻⁵ at γ² = 8×10⁻⁴). As g → 0 the
damped kernel (i/2ω̃)e^{−iω̃|t|−γ|t|/2} converges monotonically to
#213's (i/2ω)e^{−iω|t|}: **the boundary condition was the zero-coupling
shadow of the dynamics.**

## 3. The geometry adjudicates the absorber's address

Where must the absorber live — at the antipode, or distributed through
the bulk? The per-mode damping rate is γₙ = π(gYₙ(ψₐ))²ρ/(2ωₙ²), so the
placement question is a question about the mode functions — and the
geometry answers it three ways:

**The antipode is impedance-matched.** Yₙ(π) = n(−1)ⁿ⁺¹/√(2π²): the
coupling *grows* linearly with n and never vanishes, **exactly
cancelling** the 1/ωₙ² kinematic suppression:

```
γₙ = g²/(4π)   —  the same rate for every tower mode.
```

A single point absorber at the antipode damps the entire tower at one
flat rate. This is the #166 focusing caustic in its absorber role: the
1/sinψ amplification at the antipode is precisely the impedance
matching.

**A generic point leaks.** At ψₐ = 1 the couplings oscillate as sin(n)
with a near-zero at n = 22 (|sin 22| = 0.009), and the rates fall as
1/n² — high modes are long-lived remnants and complete absorption
fails.

**A uniform distribution is forbidden by a selection rule.** Exactly:
∫₀^π sin(nψ) sinψ dψ = (π/2)δₙ₁ — a uniformly distributed absorber
couples *only to the tower ground mode* and can absorb nothing above
it.

**The live adjudication** (24 tower modes kicked by a localized source
pulse + a 3000-oscillator bank, exact normal-mode evolution to t = 600):

| placement | residual field energy E_f(600)/E_f(0) |
|---|---:|
| antipodal point | **0.014** (= e^{−γt} predicted 0.0136; per-mode residuals uniform, as the flat-rate law demands) |
| generic point (ψₐ = 1) | 0.95 (near-decoupled modes survive) |
| uniform distribution | 1.00 (everything above n = 1 untouched) |

## 4. What ε *means* on a closed bulk

A finite bank revives: the field energy returns near the Poincaré time
T_rec = 2π/dν, measured at 383 vs 342 predicted (N = 300) and 726 vs
684 (N = 600) — revival time linear in N (ratio 1.90). So on the closed
bulk, **finite ε is Poincaré recurrence**: the absorber holds the
energy for one recurrence time, not forever. The ε → 0⁺ limit of #213
is the continuum limit N → ∞ of the throat's internal spectrum taken
*before* t → ∞ — an order of limits with physical content, not a
formality.

## 5. Honest scope

- The absorber is Gaussian/bilinear **by design** — that is what makes
  the elimination exact. Anharmonic absorber dynamics (real
  registration/irreversibility) is the standing #209 open, not claimed
  here.
- The flat bank density is a modeling choice; the physical absorber is
  the throat's internal continuum. Deriving its spectral density from
  the throat geometry (the quasinormal spectrum of the Tangherlini
  mouth) is the named successor.
- The source point ψ = 0 shares the antipode's matching property
  (Yₙ(0) = n/√(2π²)): source and antipode are the two distinguished
  points of the zonal geometry; the antipode is the unique *other* one.
- Classical throughout; per-mode field; frozen geometry, no
  backreaction.
- The adjudication compares one point against one distribution at a
  fixed budget; a dense random cloud of point absorbers also works
  generically. The sharp statements are the flat antipodal rate law and
  the exact uniform selection rule.

## 6. What would falsify this

- A placement-independent decay rate — the impedance-matching claim
  would be empty. (Checked: three placements differ by factors of
  ~70×.)
- Antipodal per-mode rates that vary with n — the flat-rate law would
  fail. (Checked: per-mode residuals within [0.7, 1.03]× of the
  uniform prediction.)
- A γ that does not scale as g², or disagrees among the pole, the live
  decay, and the golden rule — ε would not be the absorber response.
  (Checked: slope 1.9996; three-way agreement at 1%.)
- Revival times not scaling with N — the recurrence reading of finite ε
  would fail. (Checked: linear, ratio 1.90/2.0.)

## 7. Companion probe

`experiments/closure_ledger/dynamical_absorber_propagator_probe.py`
(T1–T8; ~1 min) machine-checks: the resolvent identity against direct
matrix inversion (10⁻¹⁰) and absorber stability; the three-way γ
agreement and exact g² scaling; the antipodal limit, the uniform
selection integrals (10⁻⁸), the generic near-zero, and the live
three-way adjudication; the revival times and their linear scaling; the
O(γ²) pole-form agreement, the derived-ε ordering quadrature, and the
monotone g → 0 recovery of #213.

**Verdict:**
`EPSILON_IS_THE_ABSORBER_RESPONSE_RATE_AND_THE_ANTIPODE_IS_THE_IMPEDANCE_MATCHED_ABSORBER_OF_THE_CLOSED_BULK`

# The coupled 5D EKG weld: the scales lock, and the lock refutes self-sourcing (PR #222)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. This PR is the one place the
> geometry is deliberately *unfrozen*: the requested modification
> couples the soliton core's localized energy density to the local 5D
> Tangherlini metric, exactly to derive the cross-sector scale weld the
> #221 identifiability audit demanded — and to test what the
> frozen-background reading posits. The companion probe machine-checks
> every claim (~6 min).

## 0. The request, executed

**Option A (the core).** The action is modified so the soliton's energy
density directly sources the local 5D metric: the static spherical
Einstein–Klein–Gordon system

```
ds² = −α²(r) dt² + a²(r) dr² + r² dΩ₃²,     Φ = φ(r) e^{−iωt}
```

solved by shooting (the #210 solver generalized to n = D − 2 sphere
dimensions; n = 2 is *exactly* the #210 system, n = 3 is 5D). The
exterior is automatically Tangherlini — N = 1 − μ/r², with
μ(r) = r²(1 − 1/a²) the 5D mass function and **r_s = √μ_∞** the throat
radius the exterior geometry encodes. The background cavity scale R\*
(the wave potential's structure on the *solved* metric), the soliton's
R_RMS (the ρ³ 5D measure), and r_s all emerge from **one coupled
solve**.

**Option B (the verification layer).** The Israel junction at the core
boundary is machine-checked: the interior numeric extrinsic curvature
matches the exterior Tangherlini's on the shared hypersurface, with
jumps tracking the scalar tail (no shell in the smooth-matching limit),
and the mass the exterior reads equals the integrated interior energy.

## 1. The solver, benchmarked

- The n = 2 reduction reproduces the **Kaup point: M_max = 0.6327**
  vs the classic 0.633 (the #210 benchmark, same machinery).
- The tt Einstein (mass-function) identity μ′ = (2κ/3) r³ ρ holds
  pointwise on the stored 5D profile (median residual 3×10⁻⁶) — an
  independent implementation check of the reduced equations.
- μ is dx-converged (halving shift 2.5×10⁻⁵) and xmax-converged
  (< 10⁻⁵).

## 2. The 5D family and the critical-mass marginality

Along the main branch (φ_c = 0.0025 → 0.7) μ decreases monotonically;
the frequency minimum is passed near φ_c ≈ 0.7 and the spiral turns
(ω_phys rises again at φ_c = 1.0, 1.4). The dilute endpoint is the
measured headline:

```
mu(phi_c -> 0)  ->  mu_crit = 7.695     (linear in phi_c, difference ratio 2.0)
omega -> 1,   X = sigma_RMS * omega  ~  phi_c^(-1/2)   (slope -0.503)
```

**The 5D dilute branch has a critical mass with a scale-free size zero
mode**: the mass does not vanish as the star dilutes — the size runs
free at fixed mass. The reason is structural: the 5D Newtonian
potential is 1/r², *marginal* — kinetic and binding energies scale
identically — the very same marginality as the Tangherlini exterior
tails that shaped #221's cavity measurements. A 5D self-gravitating
scalar is never light compared to its own gravitational radius.

## 3. The weld — the audit's T5, answered

In the coupled system every length is expressed in the scalar-mass unit
1/m, and the gravitational-coupling convention **drops out exactly**:
the rescale (KF, φ) → (KF/4, 2φ) leaves every geometric observable
(ω_phys, μ, R_RMS) unchanged to machine zero. The family relations
μ(ω), X(σ/r_s), r_s·ω are therefore **derived, convention-free
functions** — the audit's free radial-unit rescale is forbidden by the
field equations. The weld exists; the question becomes what it says.

The Israel junction confirms the shared boundary: at the core edge
(scalar tail φ/φ₀ ~ 10⁻⁴) the extrinsic-curvature jumps are
[K^θ_θ] ~ 10⁻¹⁶ and [K^t_t] ~ 5×10⁻⁷ (each tracking the tail), and
μ(R_b) = μ_∞ to 10⁻⁶: the interior soliton core and the exterior 5D
bulk read one and the same set of scales off one physical boundary.

## 4. The cavity from the same metric

The test-wave potential on the solved background (FD-constructed by
definition; validated on vacuum Tangherlini against the **#215 closed
form to 2×10⁻⁶**) has a well produced by the soliton's lapse:

| φ_c | R\* (well centroid) | width | R\*/r_s | R\*/R_RMS |
|---|---|---|---|---|
| 0.5 | 1.83 | 1.13 | 0.78 | 0.39 |
| 0.18 | 2.19 | 1.66 | 0.89 | 0.27 |
| 0.05 | 2.69 | 2.67 | 1.00 | 0.17 |

R\* and R_RMS emerge from one system — **and they decouple toward the
dilute end**: R\* tracks the *throat* scale r_s while the matter cloud
spreads. The cavity belongs to the gravitational core, not to the
cloud — exactly the structure the #221 eigenhistory reading needs (its
cavity is the throat's), now seen in the coupled solve.

## 5. The confrontation: the lock excludes self-sourcing

The convention-free lock, measured over the entire family (spiral
turnaround included; q-channel spot checks in the same band):

```
r_s · ω  ∈  [1.53, 2.774]          (sup = sqrt(mu_crit))
```

The EM-cap primordial anchor (#210: r_s = α·λ̄_C) requires
r_s·ω = α = 0.0073 — **excluded by ×210**. Jointly: σ/r_s = 206.8
(conv B) is reachable only on the dilute branch, where
X = σ·ω = 206.8 × 2.774 = **574** vs the required 1.51 — excluded
×380. Stated physically:

> **A 5D self-gravitating scalar cannot be much lighter than the
> throat its own field creates.** The reading in which the electron's
> own energy density generates its throat is refuted at the fully
> coupled level.

What this forces is exactly #210's relocation, now derived rather than
posited: **the throat is primordial** — its bulk mass μ = r_s² a
geometric datum, not the particle's self-field — and the #221 cavity
is the *primordial throat's* (consistent with §4: R\* tracks r_s, not
the cloud). The identifiability program continues on the fixed-μ
bridge: the soliton as a perturbative dressing of the primordial
throat, solved with this exact machinery — the named successor.

## 6. Honest scope

- The stationary complex-scalar (boson-star) ansatz stands in for the
  #180 real ψ–φ–q structure; the #210 q-channel potential is
  spot-checked in 5D and lands in the same r_s·ω band — the conclusion
  is a property of the potential class (the #210 compactness
  argument), not of the Kaup choice.
- One-mouth trivial topology: no cross-cap, no Pin-odd winding. The
  refutation is of *self-sourcing*; the bridge version (fixed-μ
  primordial throat + soliton dressing) is the successor.
- Ground states only; the spiral is entered but not traced — deeper
  tracing can only move the compact edge of the band *inside*
  [1.4, 2.85], never to α.
- The cavity potential is for test waves (no probe-wave back-reaction);
  R\* definitional band (centroid/minimum/width) reported.
- Classical throughout; α imported (#184) for the confrontation only.
- This PR intentionally unfreezes the geometry; the result *confirms*
  the frozen primordial background as the only consistent reading of
  the anchor — the program's framing survives its own stress test.

## 7. What would falsify this

- The n = 2 reduction missing the Kaup point, or the constraint
  identity failing — the equations would be wrong. (Checked: 0.6327;
  3×10⁻⁶.)
- Convention dependence surviving in the geometric observables — no
  weld. (Checked: exact invariance, machine zero.)
- A shell at the core boundary (junction jumps not tracking the tail)
  — the interior and exterior would not share scales. (Checked:
  jumps ~ tail, mass matched to 10⁻⁶.)
- A family member with r_s·ω ≪ 1 — self-sourcing viable after all.
  (Checked: bounded below by 1.53 across the family and the potential
  class; the dilute end *raises* it to 2.774.)
- μ vanishing on the dilute branch (4D-like behavior) — no critical
  mass, and the confrontation would soften. (Checked: μ → 7.695 with
  difference ratio 2.0.)

## 8. Companion probe

`experiments/closure_ledger/coupled_5d_ekg_weld_probe.py` (T1–T9,
~6 min): the benchmarked general-n solver; the 5D family with the
critical-mass extrapolation and the scale-free zero mode; the exact
convention invariance; the Israel junction; the vacuum-validated
cavity potential with the R\*/R_RMS decoupling; the confrontation.

**Verdict:**
`THE_WELD_IS_DERIVED_AND_IT_REFUTES_SELF_SOURCING_THE_COUPLED_5D_EKG_LOCKS_RS_OMEGA_INTO_1P5_TO_2P8_NEVER_ALPHA_THE_THROAT_IS_PRIMORDIAL_AND_THE_CRITICAL_MASS_MARGINALITY_IS_MEASURED`

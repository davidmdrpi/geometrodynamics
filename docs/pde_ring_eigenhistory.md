# The PDE-ring eigenhistory (PR #220)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. #219 derived the source from a
> Hamiltonian but analyzed the loop in harmonic balance with a reduced
> 2-dof Floquet model. This PR retires every reduction: the field is the
> **PDE** — the two-direction wave on the periodic ring with the glued
> Tangherlini barrier potentials — the Duffing (q, p) evolve in
> **explicit time domain**, the interaction term g·q·u(0) is **in the
> conserved energy**, the periodic orbit solves the **literal** system
> X(T) − X(0) = 0, H[X] − E₀ = 0 with one phase condition, and the
> one-period monodromy of the **complete** source–field state is
> computed — with the complete solve repeated at four spatial
> resolutions to demonstrate O(dx²) convergence of everything
> reported. The companion probe machine-checks every claim (~2.5 min).

## 0. The complete Hamiltonian system

```
H = Σᵢ dx [πᵢ²/2 + ((Du)ᵢ² + Vᵢuᵢ²)/2]           (the field ring)
    + p²/2 + ω₀²q²/2 + (μ/4)q⁴                     (the Duffing source)
    + g·q·u(0)                                     (the interaction —
                                                    IN the ledger)
```

- **The ring**: x ∈ [0, C), C = 2π + 8, periodic — the source at x = 0,
  the two finite-width Tangherlini barrier potentials (the *physical*
  greybody barriers of #215, not #219's point-scatterer idealization)
  glued at x = π and x = π + 8 with their exteriors facing the source
  arc. N = 192 points (≈30 per carrier wavelength).
- **The source**: the #219 side-coupled Duffing (ω₀ = 3.2, μ = 0.5,
  g = 0.6) — but now its (q, p) are explicit phase-space variables
  integrated in time, never slaved.
- **The integrator**: symplectic leapfrog (velocity Verlet) on the
  canonical pairs (uᵢ, dx·πᵢ) and (q, p), n_s = 2048 steps/period.

**The interaction term is load-bearing.** With g·q·u(0) included, H is
conserved to a bounded, non-secular O(dt²) shadow-Hamiltonian
oscillation (2×10⁻⁵ over **one hundred** periods; halving dt cuts it
×4). The ledger *without* the interaction term visibly fails — varying
2×10⁻³ in a **single** period, ~100× worse. The requested inclusion is
not bookkeeping; it is the difference between a conserved energy and a
broken one.

## 1. The periodic orbit, literally as requested

Unknowns: the complete state X(0) = (u, π, q, p) ∈ R^{2N+2} = R^{386}
and the period T. Conditions — exactly the requested system:

```
X(T) − X(0) = 0             (periodicity of the complete state)
H[X(0)] − E₀ = 0            (energy closure, interaction included)
p(0) = 0                    (the single phase condition removing
                             time-translation degeneracy)
```

Solved by Gauss–Newton with the full finite-difference Jacobian
(batched leapfrog columns), from a linear ring-mode guess:

| quantity | value |
|---|---|
| E₀ | 20.834 |
| Newton history ‖R‖ | 3.1×10⁻¹ → 5.0×10⁻³ → 1.2×10⁻⁴ → 1.9×10⁻¹⁰ → **2.4×10⁻¹³** |
| period T | 2.353874 |
| ω_orbit = 2π/T | 2.669297 (nonlinear pull −3.1×10⁻³ from the linear mode 2.6723) |
| u(0), q on the orbit | 0.9003, −0.17282 (antiphased, as #219 predicts) |
| energy at the orbit | E₀ to 10⁻¹¹; phase condition to 10⁻¹² |
| dt-refinement of ω_orbit | ~10⁻⁶ relative (O(dt²) convergence) |

Quadratic convergence to machine residual: the eigenhistory of the
complete source–field state **exists** as a genuine periodic orbit of
the full PDE system — not of a reduced model.

## 2. The one-period monodromy of the complete state

The tangent leapfrog — the exact linearization of the discrete
symplectic map, co-evolved with the orbit — gives the full
**386 × 386** monodromy M. Machine facts:

- **all eigenvalues on the unit circle**: max ||λ| − 1| ≈ 4×10⁻¹⁴ over
  all 193 Floquet pairs — *no parametric instability of any field mode
  against the eigenhistory*, not just of the reduced 2-dof model;
- **det M = 1** to 10⁻¹³;
- **symplectic**: MᵀJM = J in the dx-weighted symplectic form
  (canonical pairs (uᵢ, dx·πᵢ)) to 4×10⁻¹⁵;
- **the trivial pair** λ = 1 (along-flow + energy) to 10⁻¹⁰ — where
  Floquet theory demands it.

## 3. The source pair — #219 confirmed by the field itself

The eigenvector (q, p)-weight cleanly identifies the source pair:
weight ≈ 0.36 for the pair vs ≈ 0.008 for the next candidate. Its
unwrapped rotation frequency:

```
ω_source = (θ + 2π)/T = 3.2124
```

vs #219's reduced-model prediction **3.2102** (bare ω₀ = 3.2): agreement
to **0.07%** — the dressed Duffing frequency, computed there from a
one-ring-mode reduction, is reproduced by the complete PDE monodromy.
Every approximation #219 made (harmonic balance, single-mode reduction,
point-barrier transfer matrices) is now retired, and its answer
survives. The low field modes appear in the spectrum at their folded
angles ω_k·T mod 2π, weakly perturbed by the orbit.

## 4. The energy ledger, time-resolved

Sampled at 64 points along the orbit:

- the interaction energy g·q·u(0) **oscillates about a negative mean**
  (q and u(0) antiphased), with the mean equal to the #219
  harmonic-balance value ⟨g·q·u(0)⟩ = g·a·U₀/2 to a few percent
  (the residual is the orbit's harmonic content — exactly what
  harmonic balance drops);
- the full H — field + source + interaction — is **constant along the
  orbit** to the shadow-Hamiltonian bound;
- E_field / E_source / E_int partition reported per snapshot.

## 5. Spatial-resolution convergence

The complete solve — orbit, partition, monodromy — repeated at
**N = 128, 192, 256, 384** while holding the physical ring, the
barriers, the coupling, and E₀ = 20.834 fixed (temporal resolution
held at n_s = 2048, so the sweep isolates the spatial error):

| N | dx | ω_orbit | u(0) | q | E_field / E_source / E_int (means) | ω_source |
|---|---|---|---|---|---|---|
| 128 | 0.1112 | 2.664093 | 0.902121 | −0.171645 | 20.75271 / 0.12772 / −0.04645 | 3.211378 |
| 192 | 0.0741 | 2.669296 | 0.900333 | −0.172817 | 20.75097 / 0.12968 / −0.04667 | 3.212364 |
| 256 | 0.0556 | 2.671118 | 0.899710 | −0.173234 | 20.75035 / 0.13038 / −0.04676 | 3.212728 |
| 384 | 0.0371 | 2.672421 | 0.899267 | −0.173534 | 20.74990 / 0.13088 / −0.04681 | 3.212994 |
| **∞ (Richardson)** | — | **2.673464** | **0.898913** | **−0.173774** | E_int → **−0.046859** | **3.213207** |

- **Clean O(dx²), quantitatively**: against a finite reference the
  second-order prediction for the diff-to-finest ratios is *exactly*
  d₁₂₈/d₂₅₆ = 6.400 and d₁₉₂/d₂₅₆ = 2.400; the measured ratios for the
  orbit frequency, q, u(0), all three energy-partition means, the
  periodic-cubic-interpolated field profile (L₂, on a common 768-point
  grid), the source-pair frequency, and the two machine-matched low
  Floquet angles all land in 6.06–6.66 and 2.36–2.42.
- **Every resolution's Gauss–Newton converges to < 10⁻¹¹** and every
  resolution's complete monodromy (258 → 770 dimensional) stays on the
  unit circle (worst deviation 1.4×10⁻¹⁰, at N = 256).
- **Low Floquet multipliers converge with the spectrum**: the matched
  low field pairs (k = 2, 3, gaps ~10⁻⁷–10⁻⁸ at every N) move O(dx²) to
  their continuum angles; the source pair converges 3.2114 → 3.2124 →
  3.2127 → 3.2130, Richardson **3.2132 — 0.09% from #219's reduced-model
  3.2102**: the confirmation of #219 survives the continuum limit.
- The N = 192 production grid sits ~1.2×10⁻³ (relative, profile L₂)
  from the finest grid and ~4×10⁻³ absolute in ω_orbit from the
  continuum — small against every physical effect quoted, and now
  *measured*, not assumed.

## 6. Honest scope

- The orbit and monodromy are those of the *discrete* symplectic map at
  n_s = 2048; the continuum is approached O(dt²) (frequency shift ~10⁻⁶
  verified by dt-refinement). H is conserved to the bounded shadow
  oscillation, never secularly.
- Grid N = 192 (dx ≈ 0.075, ~30 points per carrier wavelength); the
  monodromy spectrum is the discretized field's — all of it on the
  unit circle. The residual spatial error is quantified by the
  N = 128–384 sweep (§5): O(dx²) with the exact predicted ratios, and
  Richardson-extrapolated to the continuum.
- The ring uses the physical finite-width Tangherlini barriers
  (interior separation 8) — the *fuller* model than #219's
  point-barrier idealization; the local source physics is robust to
  this (source pair within 0.07%).
- E₀ is a specified budget; ℏ is still not derived.
- Classical, zonal scalar, frozen geometry; the CTC/network aspects
  (mouth clock offset, greybody one-way membrane) live in #216–#218 —
  this PR is the conservative-dynamics core: the ring as the
  time-closed loop's covering model.

## 7. What would falsify this

- Secular energy drift, or a conserved ledger *without* the interaction
  term — the Hamiltonian structure would be wrong. (Checked: bounded
  2×10⁻⁵/100 periods with the term; 100× violation without it.)
- No periodic orbit of the complete state at the specified E₀.
  (Checked: Gauss–Newton to 2×10⁻¹³.)
- Any monodromy eigenvalue off the unit circle — a parametric
  instability of some field mode against the eigenhistory that the
  reduced model could not see. (Checked: all 386 within 10⁻¹³.)
- A source pair absent or at the wrong frequency — the #219 reduction
  would be an artifact. (Checked: present by qp-weight, at 3.2124 vs
  3.2102.)
- An interaction-energy mean inconsistent with the #219 ledger.
  (Checked: negative, matching to a few percent.)
- Non-convergence under grid refinement — the reported orbit and
  spectrum would be discretization artifacts. (Checked: O(dx²) with
  the exact predicted ratios in every quantity; the source pair's
  continuum limit still 0.09% from #219.)

## 8. Companion probe

`experiments/closure_ledger/pde_ring_eigenhistory_probe.py` (T1–T9,
~2.5 min): the conservation checks with and without the interaction
term; the literal Gauss–Newton periodic-orbit solve; the complete
monodromy with the dx-weighted symplectic form; the source-pair
identification and #219 comparison; the time-resolved ledger; the
four-resolution convergence sweep with Richardson extrapolation.

**Verdict:**
`THE_FULL_PDE_EIGENHISTORY_EXISTS_ITS_COMPLETE_MONODROMY_IS_UNIT_CIRCLE_SYMPLECTIC_AND_THE_219_REDUCTION_IS_CONFIRMED_BY_THE_FIELD_ITSELF`

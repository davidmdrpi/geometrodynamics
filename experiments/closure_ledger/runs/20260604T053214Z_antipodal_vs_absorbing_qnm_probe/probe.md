# Antipodal vs absorbing throat quasinormal spectrum (PR #130)

**Run:** 2026-06-04T05:32:14+00:00

Computes the full frequency spectrum of the BAM cavity under the two throat boundary conditions — the BAM-native antipodal (PR #129) and the ordinary absorbing horizon — and contrasts them. The antipodal throat gives real, undamped normal modes (stable matter); the absorbing horizon gives complex, damped quasinormal ringdown.

- **Antipodal**: real ω (Im ≈ 0), undamped normal modes, Q = ∞, l-parity graded
- **Absorbing**: complex ω = ω_R − i|ω_I| (Im < 0), damped ringdown, Q ~ O(1)
- **Consequence**: stable matter (sharp real masses) needs the unitary antipodal throat (#64/#129)
- **Open**: idealised r* → −∞ horizon QNMs; GW coupling; absolute normalisation

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | compute & contrast antipodal vs absorbing throat QNM spectrum | **PASS** |
| T2 | `T2_eigenproblem_setup` | setup: cavity operator, shell wall, antipodal vs ingoing throat BC | **PASS** |
| T3 | `T3_antipodal_real_undamped` | antipodal ⟹ real ω (Im≈0), undamped normal modes (Q=∞) | **PASS** |
| T4 | `T4_absorbing_complex_ringdown` | absorbing ⟹ complex ω (Im<0), damped ringdown | **PASS** |
| T5 | `T5_quality_factor_contrast` | Q = ω_R/(2|ω_I|): antipodal ∞ (sharp), absorbing O(1) (width) | **PASS** |
| T6 | `T6_stable_matter_needs_unitary_throat` | stable matter (sharp masses) needs the unitary antipodal throat | **PASS** |
| T7 | `T7_scope` | scope: finite-cavity spectrum; idealised horizon QNM / GW open | **PASS** |
| T8 | `T8_assessment` | ANTIPODAL_THROAT_REAL_UNDAMPED_VS_ABSORBING_COMPLEX_RINGDOWN | **PASS** |

## The spectra: antipodal (real) vs absorbing (complex)

| l | antipodal ω (BC) | absorbing ω (ingoing) |
|---:|---|---|
| 0 | 1.181-0.000i, 3.293-0.000i, 5.439-0.000i (N) | 1.893-1.241i, 3.917-1.166i, 6.009-1.044i |
| 1 | 2.295+0.000i, 4.363+0.000i, 6.478+0.000i (D) | 1.995-1.126i, 3.973-1.140i, 6.044-1.035i |
| 2 | 1.436+0.000i, 3.443+0.000i, 5.532-0.000i (N) | 2.168-1.000i, 4.063-1.100i, 6.102-1.021i |
| 3 | 2.594+0.000i, 4.536+0.000i, 6.595+0.000i (D) | 2.389-0.896i, 4.185-1.050i, 6.182-1.002i |

Antipodal spectrum is real to `max|Im ω| = 0.0` (undamped); absorbing modes all have `Im ω < 0` (damped ringdown).

## Quality factor: sharp normal modes vs leaky ringdown

| l | antipodal ω | Q (antipodal) | absorbing ω | Q (absorbing) | τ (absorbing) |
|---:|---|---|---|---|---|
| 0 | 1.181-0.000i | inf | 1.893-1.241i | 0.763 | 0.806 |
| 1 | 2.295+0.000i | inf | 1.995-1.126i | 0.886 | 0.888 |
| 2 | 1.436+0.000i | inf | 2.168-1.000i | 1.084 | 1.0 |

Antipodal: `Q = ∞` (sharp, infinitely-lived). Absorbing: `Q ~ O(1)` (Lorentzian width `Γ = 2|ω_I|`, lifetime `τ = 1/|ω_I|`) — the thin cavity leaks fast into the throat.

## Verdict

**ANTIPODAL_THROAT_REAL_UNDAMPED_VS_ABSORBING_COMPLEX_RINGDOWN.** THE ANTIPODAL THROAT GIVES A REAL, UNDAMPED SPECTRUM; THE ABSORBING HORIZON GIVES COMPLEX RINGDOWN. PR #129 fixed the BAM throat BC as the antipodal, l-parity-graded condition (a unitary mirror); this probe computes the full frequency spectrum for both that BC and the ordinary absorbing horizon, and contrasts them.

THE EIGENPROBLEM. The separated radial equation −d²ψ/dr*² + V_l ψ = ω²ψ runs on the BAM cavity [R_MID+ε, R_OUTER] with the shell wall (Dirichlet) at R_OUTER and, at the throat, either the antipodal real l-parity BC (Neumann for even l, Dirichlet for odd l; self-adjoint) or the absorbing ingoing BC ∂ψ(throat) = −iω ψ(throat) (non-self-adjoint). The absorbing case is a quadratic eigenvalue problem in ω, solved by companion linearisation.

ANTIPODAL ⟹ REAL, UNDAMPED. The antipodal BC makes the operator self-adjoint, so its spectrum is real: Im(ω) ≈ 0 (numerically zero). These are infinitely-lived normal modes — sharp spectral lines, quality factor Q = ∞, zero width — and they are l-parity graded (even-l Neumann, odd-l Dirichlet). The unitary mirror of PR #129 conserves the mode energy.

ABSORBING ⟹ COMPLEX, RINGDOWN. The ingoing BC makes the operator non-self-adjoint; the eigenfrequencies are ω = ω_R − i|ω_I| with Im(ω) < 0 — damped quasinormal modes. The lowest l=0 mode is ≈ 1.89 − 1.24i: a ringdown of lifetime τ = 1/|ω_I| ≈ 0.8 and quality factor Q = ω_R/(2|ω_I|) ≈ 0.8, the thin cavity leaking fast into the throat. Energy is lost into the horizon.

THE DECAY CONTRAST. The quality factor Q = ω_R/(2|ω_I|) is infinite for the antipodal throat (a sharp δ-line) but O(1) for the absorbing horizon (a Lorentzian of width Γ = 2|ω_I|, finite lifetime τ = 1/|ω_I|).

STABLE MATTER NEEDS THE UNITARY THROAT. A matter state is a sharp mass — a stable or long-lived particle — only if its cavity mode has a real frequency. The absorbing throat gives every mode a width / complex mass (a decaying resonance); only the antipodal, unitary throat (PR #129) yields the real, stable spectrum that the BAM matter sectors (the lepton/quark bound states) require. The undamped-vs-ringdown distinction is the spectral face of the program's global CPT / unitarity (PR #64): BAM matter is stable precisely because the throat reflects antipodally rather than absorbing.

SCOPE. This computes the spectrum of the FINITE BAM cavity (the physically appropriate region) under the two throat BCs. It does NOT compute the idealised r* → −∞ horizon QNMs, the coupling to gravitational radiation, or the absolute mode normalisation. The absorbing case is the counterfactual — BAM selects the antipodal BC (PR #129); the contrast shows what is at stake. The Im(ω) < 0 sign is the e^{−iωt} decay convention.

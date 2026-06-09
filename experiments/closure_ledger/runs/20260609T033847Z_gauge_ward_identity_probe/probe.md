# Gauge Ward identity and current conservation audit (PR #142)

**Run:** 2026-06-09T03:38:47+00:00

Audits the gauge consistency of the #141 coupling. The matter current is conserved at the antipodal throat (real modes carry no charge flux), the Ward–Takahashi identity ties the gauge vertex to the matter propagator, and the photon stays massless — all from the unitary antipodal throat (#129), the same postulate as stable matter. An absorbing throat would leak charge and break gauge invariance. *(QFT on the classical throat, not quantum gravity.)*

- **Noether current**: j^μ = i(ψ*∂^μψ − ...), ∂_μ j^μ = 0
- **Conservation**: real modes ⟹ j^r = Im(ψ*∂ψ) = 0 (no charge flux, #129/#135)
- **Absorbing (counterfactual)**: complex modes ⟹ j^r ≠ 0 (charge leaks) ⟹ gauge invariance broken
- **Ward–Takahashi**: q_μ Γ^μ = S⁻¹(p') − S⁻¹(p) (vertex #141 fixed by propagator #135)
- **Masslessness**: q_μ Π^μν = 0 ⟹ photon massless, 1/q² protected (#42–#44)
- **One postulate**: gauge invariance = the gauge face of the unitary mirror (#129); only α input
- **Open**: α (#105/#108); higher-order Ward identities; running of α; scale (#133); flavor (#134)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | gauge Ward identity and current conservation audit | **PASS** |
| T2 | `T2_conserved_noether_current` | Noether current j^μ conserved; stationary mode ⟹ ρ static | **PASS** |
| T3 | `T3_current_conserved_at_antipodal_throat` | real modes ⟹ j^r = Im(ψ*∂ψ) = 0 (no charge flux; #129) | **PASS** |
| T4 | `T4_absorbing_throat_breaks_conservation` | absorbing complex modes ⟹ j^r ≠ 0 (charge leaks) ⟹ breaks gauge inv. | **PASS** |
| T5 | `T5_ward_takahashi_identity` | Ward–Takahashi q_μ Γ^μ = S⁻¹(p') − S⁻¹(p) (#141 ↔ #135) | **PASS** |
| T6 | `T6_transversality_photon_massless` | transversality q_μ Π^μν = 0 ⟹ photon massless (1/q² protected) | **PASS** |
| T7 | `T7_gauge_invariance_from_unitary_mirror` | gauge invariance from the unitary mirror (#129) — same as stable matter | **PASS** |
| T8 | `T8_assessment` | GAUGE_WARD_IDENTITY_CURRENT_CONSERVED_PHOTON_MASSLESS_FROM_ANTIPODAL_MIRROR | **PASS** |

## Current conservation: antipodal (j^r = 0) vs absorbing (j^r ≠ 0)

| mode | ω | max\|j^r\| | charge flux |
|---|---|---:|---|
| antipodal n=0 (real) | 1.166 | 0.0 | none (conserved) |
| antipodal n=1 (real) | 3.25 | 0.0 | none (conserved) |
| antipodal n=2 (real) | 5.367 | 0.0 | none (conserved) |
| absorbing fundamental (complex) | 1.893-1.122i | 0.01445 | leaks into horizon |

The real antipodal modes carry **zero** radial charge current (current conserved, the unitary mirror #129); the absorbing complex modes carry a nonzero current into the horizon (conservation broken). Gauge invariance requires the antipodal throat.

## Verdict

**GAUGE_WARD_IDENTITY_CURRENT_CONSERVED_PHOTON_MASSLESS_FROM_ANTIPODAL_MIRROR.** THE MATTER CURRENT IS CONSERVED AT THE ANTIPODAL THROAT, THE WARD IDENTITY HOLDS, AND THE PHOTON STAYS MASSLESS — ALL FROM THE UNITARY MIRROR, THE SAME POSTULATE AS STABLE MATTER. PR #141 built the minimal gauge–matter coupling; this probe audits its consistency.

THE CONSERVED NOETHER CURRENT. The global U(1)_Hopf phase symmetry gives the current j^μ = i(ψ*∂^μψ − (∂^μψ*)ψ), conserved on shell ∂_μ j^μ = 0. For a stationary cavity mode the charge density ρ = 2ω_n|ψ_n|² is time-independent, so conservation reduces to the vanishing radial current.

CURRENT CONSERVATION AT THE ANTIPODAL THROAT. The antipodal cavity modes are REAL (#135), so the radial charge current j^r ∝ Im(ψ_n*∂_rψ_n) = 0 exactly (verified). No charge flows through the throat: the charge is static and conserved, a stable charged particle. Current conservation in the gauge sector IS the zero-flux unitary-mirror property (#129).

AN ABSORBING THROAT WOULD BREAK IT. An absorbing horizon gives complex quasinormal modes (#130) whose radial current j^r = Im(ψ*∂_rψ) ≠ 0 (verified ≈ −0.014 at the throat) carries charge into the horizon. Current conservation fails, and with it gauge invariance — a charged black-hole-style throat is not gauge-consistent. So gauge invariance REQUIRES the antipodal (unitary) throat, the same requirement as stable matter.

THE WARD–TAKAHASHI IDENTITY AND PHOTON MASSLESSNESS. Current conservation implies q_μ Γ^μ(p,p') = S⁻¹(p') − S⁻¹(p), tying the gauge vertex (#141) to the matter inverse propagator (#135): the gauge coupling is fixed by the matter dynamics (the differential form Γ^μ = ∂S⁻¹/∂p_μ normalises the vertex to the charge). It also makes the vacuum polarisation transverse, q_μ Π^μν = 0, so the gauge correction generates NO photon mass: the 1/q² photon (#42–#44) is protected.

ONE POSTULATE, BOTH. Current conservation, the Ward identity, and photon masslessness all follow from the unitary antipodal throat (#129) — the same real, self-adjoint, zero-flux structure that gave the stable spectrum (#130), the unitary propagator (#135), the stable self-energy (#136), and the bounded vacuum (#138). Gauge invariance is the gauge face of the unitary mirror, not an extra assumption. Only the coupling strength α (#105) stays input.

SCOPE. Audits gauge consistency on the antipodal throat and ties it to the unitary mirror. It does NOT fix α (#105) or compute higher-order Ward identities / the running of α. The α (#105/#108), bulk-scale (#133), and flavor (#134) residuals stand.

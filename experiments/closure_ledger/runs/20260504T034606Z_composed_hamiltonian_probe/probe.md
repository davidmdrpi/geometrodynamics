# Composed-Hamiltonian probe — closure quantum × radial matrix elements

**Run:** 2026-05-04T03:46:06+00:00
**Grid:** N = 80

Tests whether the locked lepton surrogate can be reconstructed from (closure-quantum uplift) + (Tangherlini radial matrix elements) alone, with **no other fit parameters**. The closure quantum is β = 50.0·π (integer 4β/2π = 100 antipodal cavity quanta); the radial matrix elements come from the canonical Tangherlini eigensolver on the tortoise grid. Mass extraction is linear in eigenvalue: m_i / m_e = λ_i / λ_0 (the same convention `compute_knotted_lepton_spectrum` uses).

## Available radial matrix elements

| l | ω(l, 0) | ω² | ⟨u_l|V_l|u_l⟩ |
|---|---:|---:|---:|
| 1 | 1.0547 | 1.1124 | 0.2047 |
| 3 | 1.2191 | 1.4862 | 0.4859 |
| 5 | 1.3960 | 1.9487 | 0.7523 |

Locked-surrogate diagonal target: [6.876, 34.91, 694.69]. The τ row is dominated by 4·β = 200π ≈ 628.3, which the closure quantum supplies on its own. The k=3 row needs ~33 in addition to ω² ≈ 1.5; no single radial matrix element listed above is anywhere near that scale.

## Variant results

| variant | eigenvalues | mass ratios (μ/e, τ/e) | log err (μ, τ) | match? | B1 spread | D1 spread |
|---|---|---|---|---|---:|---:|
| `HC_1_omega_sq_plus_closure` | [1.112, 1.486, 630.267] | (1.34, 566.56) | (-2.19, -0.79) | no | 0.381 | 0.000 ⚠ |
| `HC_2_omega_sq_plus_closure_with_d1` | [0.884, 1.713, 630.269] | (1.94, 713.08) | (-2.03, -0.69) | no | 0.286 | 0.660 |
| `HC_3_action_base_plus_omega_plus_closure_with_d1` | [7.167, 7.997, 636.552] | (1.12, 88.82) | (-2.27, -1.59) | no | 0.286 | 0.660 |
| `HC_4_action_base_times_omega_plus_closure_with_d1` | [6.932, 9.394, 640.564] | (1.36, 92.41) | (-2.18, -1.58) | no | 0.373 | 0.223 |
| `HC_5_v_diag_plus_closure_with_d1` | [-0.051, 0.741, 629.072] | (0.00, 0.00) | (+inf, +inf) | tachyonic | 0.269 | 0.691 |
| `HC_6_action_base_plus_v_diag_plus_closure_with_d1` | [6.232, 7.024, 635.355] | (1.13, 101.95) | (-2.26, -1.53) | no | 0.269 | 0.691 |
| `HC_7_v_cross_plus_closure_with_d1` | [1.481, 2.377, 629.543] | (1.61, 425.16) | (-2.11, -0.91) | no | 0.305 | 0.609 |
| `HC_8_action_squared_plus_closure_with_d1` | [43.907, 58.680, 705.253] | (1.34, 16.06) | (-2.19, -2.34) | no | 0.381 | 0.038 ⚠ |

Observed: m_μ/m_e = 206.77, m_τ/m_e = 3477.23.

### Per-variant detail

#### `HC_1_omega_sq_plus_closure`

Minimal composition: H_ii = ω_i² + β·max(0, k_i−3)². No off-diagonal. Tests whether the closure quantum alone, acting on the radial-eigenfrequency diagonal, recovers the ladder.

- Diagonal: `ω_i² + β·max(0,k_i−3)²`
- Off-diagonal: `0`
- Eigenvalues: ['1.1124', '1.4862', '630.2673']
- Mass ratios (predicted) m_e:m_μ:m_τ = 1.0000 : 1.3359 : 566.5584 (observed 1 : 206.77 : 3477.23)
- log₁₀ ratio error: μ = -2.1897, τ = -0.7880
- B1 closure: residues ['0.8819', '0.7707', '0.7605'] π, circular spread 0.381385 rad
- D1 closure: residues ['0.0000', '0.0000', '0.0000'] π, circular spread 0.000000 rad
- ⚠ **Near-identity regime**: max |v_ij| (i≠j) = 0.0000. The small D1 spread is a structural artifact (near-identity eigenvectors × D1's zero diagonal), not a real closure mechanism.

#### `HC_2_omega_sq_plus_closure_with_d1`

HC_1 with D1 transport off-diagonal: H_ij = (ω_j²−ω_i²)·⟨u_i|u_j⟩ for i ≠ j (mirrored to symmetric).

- Diagonal: `ω_i² + β·max(0,k_i−3)²`
- Off-diagonal: `D1 = (ω_j²−ω_i²)·⟨u_i|u_j⟩`
- Eigenvalues: ['0.8839', '1.7134', '630.2686']
- Mass ratios (predicted) m_e:m_μ:m_τ = 1.0000 : 1.9385 : 713.0825 (observed 1 : 206.77 : 3477.23)
- log₁₀ ratio error: μ = -2.0280, τ = -0.6881
- B1 closure: residues ['0.8514', '0.8012', '0.7605'] π, circular spread 0.285612 rad
- D1 closure: residues ['1.8945', '0.1046', '0.0009'] π, circular spread 0.660102 rad

#### `HC_3_action_base_plus_omega_plus_closure_with_d1`

HC_2 with the S³ great-circle action shifted into the diagonal: H_ii = 2π + ω_i² + β·max(0,k_i−3)².

- Diagonal: `2π + ω_i² + β·max(0,k_i−3)²`
- Off-diagonal: `D1`
- Eigenvalues: ['7.1671', '7.9966', '636.5518']
- Mass ratios (predicted) m_e:m_μ:m_τ = 1.0000 : 1.1157 : 88.8164 (observed 1 : 206.77 : 3477.23)
- log₁₀ ratio error: μ = -2.2679, τ = -1.5927
- B1 closure: residues ['0.8514', '0.8012', '0.7605'] π, circular spread 0.285612 rad
- D1 closure: residues ['1.8945', '0.1046', '0.0009'] π, circular spread 0.660102 rad

#### `HC_4_action_base_times_omega_plus_closure_with_d1`

Multiplicative action_base: H_ii = (2π)·ω_i² + β·max(0,k_i−3)². Mirrors the locked surrogate's per-row weighting of the radial eigenfrequency.

- Diagonal: `2π·ω_i² + β·max(0,k_i−3)²`
- Off-diagonal: `D1`
- Eigenvalues: ['6.9318', '9.3944', '640.5642']
- Mass ratios (predicted) m_e:m_μ:m_τ = 1.0000 : 1.3553 : 92.4101 (observed 1 : 206.77 : 3477.23)
- log₁₀ ratio error: μ = -2.1835, τ = -1.5755
- B1 closure: residues ['0.8793', '0.7733', '0.7605'] π, circular spread 0.373312 rad
- D1 closure: residues ['1.9640', '0.0351', '0.0009'] π, circular spread 0.223460 rad

#### `HC_5_v_diag_plus_closure_with_d1`

Potential expectation: H_ii = ⟨u_i|V_i|u_i⟩ + β·max(0,k_i−3)². Replaces ω² by the bound-state mean potential.

- Diagonal: `⟨u_i|V_i|u_i⟩ + β·max(0,k_i−3)²`
- Off-diagonal: `D1`
- Eigenvalues: ['-0.0514', '0.7407', '629.0722']
- Mass ratios (predicted) m_e:m_μ:m_τ = 0.0000 : 0.0000 : 0.0000 (observed 1 : 206.77 : 3477.23)
- log₁₀ ratio error: μ = +inf, τ = +inf
- B1 closure: residues ['0.8461', '0.8065', '0.7605'] π, circular spread 0.268958 rad
- D1 closure: residues ['1.8896', '0.1096', '0.0009'] π, circular spread 0.691126 rad

#### `HC_6_action_base_plus_v_diag_plus_closure_with_d1`

HC_5 with the S³ great-circle action: H_ii = 2π + ⟨u_i|V_i|u_i⟩ + β·max(0,k_i−3)².

- Diagonal: `2π + ⟨u_i|V_i|u_i⟩ + β·max(0,k_i−3)²`
- Off-diagonal: `D1`
- Eigenvalues: ['6.2318', '7.0239', '635.3554']
- Mass ratios (predicted) m_e:m_μ:m_τ = 1.0000 : 1.1271 : 101.9544 (observed 1 : 206.77 : 3477.23)
- log₁₀ ratio error: μ = -2.2635, τ = -1.5328
- B1 closure: residues ['0.8461', '0.8065', '0.7605'] π, circular spread 0.268958 rad
- D1 closure: residues ['1.8896', '0.1096', '0.0009'] π, circular spread 0.691126 rad

#### `HC_7_v_cross_plus_closure_with_d1`

Full cross-shell potential expectation: H_ii = Σ_l ⟨u_i|V_l|u_i⟩ + β·max(0,k_i−3)². Tests whether the muon lift comes from sampling all three l-potentials in the bound state.

- Diagonal: `Σ_l ⟨u_i|V_l|u_i⟩ + β·max(0,k_i−3)²`
- Off-diagonal: `D1`
- Eigenvalues: ['1.4807', '2.3775', '629.5432']
- Mass ratios (predicted) m_e:m_μ:m_τ = 1.0000 : 1.6056 : 425.1608 (observed 1 : 206.77 : 3477.23)
- log₁₀ ratio error: μ = -2.1098, τ = -0.9127
- B1 closure: residues ['0.7949', '0.8577', '0.7605'] π, circular spread 0.305366 rad
- D1 closure: residues ['1.9026', '0.0965', '0.0009'] π, circular spread 0.609214 rad

#### `HC_8_action_squared_plus_closure_with_d1`

Quadratic action: H_ii = (2π)²·ω_i² + β·max(0,k_i−3)². Tests whether the larger prefactor on ω_i² lifts the muon eigenvalue without introducing new parameters.

- Diagonal: `(2π)²·ω_i² + β·max(0,k_i−3)²`
- Off-diagonal: `D1`
- Eigenvalues: ['43.9075', '58.6803', '705.2531']
- Mass ratios (predicted) m_e:m_μ:m_τ = 1.0000 : 1.3365 : 16.0623 (observed 1 : 206.77 : 3477.23)
- log₁₀ ratio error: μ = -2.1895, τ = -2.3354
- B1 closure: residues ['0.8818', '0.7708', '0.7605'] π, circular spread 0.381165 rad
- D1 closure: residues ['1.9935', '0.0057', '0.0008'] π, circular spread 0.038351 rad
- ⚠ **Near-identity regime**: max |v_ij| (i≠j) = 0.0251. The small D1 spread is a structural artifact (near-identity eigenvectors × D1's zero diagonal), not a real closure mechanism.

## Verdict

**No variant reaches factor-of-5 agreement with the observed lepton mass ratios.** The closure quantum lifts the τ row correctly (since 4β = 100·(2π) is geometric and dominates the τ diagonal), but the μ row cannot be explained from radial matrix elements alone — its lift in the locked surrogate (γ_pinhole ≈ 22.5, not a clean geometric integer) does not have a clean origin in this catalog.

### Where the gap sits

All composed variants share the same τ-row prediction once the closure quantum is included: λ_τ ≈ 4β + (radial diagonal for k=5), which evaluates to roughly 630–640 across the catalog. The locked surrogate's λ_τ ≈ 695 is reached when additional non-radial pieces (action_base, pinhole γ at k=5) are summed in. The dominant τ-mass contribution is geometric.

The μ-row prediction varies more across variants because the k=3 entry has no closure-quantum contribution (since max(0, 3−3)² = 0). The radial diagonal for k=3 is at most O(2π) ≈ 6.3 in the most generous variant (HC_4 with action_base × ω_i²: 2π · 1.49 ≈ 9.4). The locked surrogate's λ_μ ≈ 41 is reached only because the pinhole γ ≈ 22.5 lifts the k=3 row substantially. None of HC_1–HC_8 supplies an O(20) diagonal contribution at k=3 from radial matrix elements alone.

The μ-row gap is therefore **the structural finding** of this probe: closure-quantum uplift + radial matrix elements is **necessary but not sufficient** for the lepton ladder. The missing piece — a pinhole-like γ ≈ 22.5 active at k = 3, 5 — does not have an obvious natural origin in the canonical Tangherlini machinery on the (l, 0) ground modes.

### Closure-ledger side effect

Even though the mass ladder isn't reproduced, the probe's eigenvectors give well-defined closure spreads. The spreads are roughly the same across variants because all share the τ-dominated-by-closure-quantum block structure: τ is near-pure depth-5, so its eigenvector projection of B1's radial phase is dominated by Φ(l=5, n=0) ≈ 0.76π. The e/μ rows pick up small mixing-induced shifts. This confirms that the closure spread is set by the per-mode B1/D1 phases on the depth basis, not by the diagonal Hamiltonian structure.

**2 variant(s) sit in the near-identity D1 regime** (max |v_ij| < 0.10 for i ≠ j and D1 spread < 0.10 rad), flagged with ⚠ in the status table. Their tight D1 numbers are a structural artifact: when the diagonal of H dominates the off-diagonal, eigenvectors are approximately the standard basis; D1 has Φ_ii = 0; v^T D1 v is therefore tiny by construction, not because the operator structure closes the ledger. These variants do NOT count as evidence for closure.

### Conclusion

The locked lepton surrogate is **not** fully reconstructible from closure quantum + Tangherlini radial matrix elements. The closure quantum (β = 50π, integer 4β/2π = 100) cleanly supplies the τ row; the radial matrix elements give the e row approximately right (eigenfrequency O(1) across all variants); but the μ row is a factor of 5–30 too low because no radial matrix element on the (l, 0) ground basis matches the locked surrogate's pinhole γ ≈ 22.5 at k=3. That value is a fitted phenomenological lift in the existing surrogate, and this probe does not find a natural geometric replacement for it. Closing the lepton ladder structurally requires either a second closure-quantum channel active at k=3 OR a non-ground radial mode whose matrix elements can supply an O(20) diagonal contribution at the muon row.
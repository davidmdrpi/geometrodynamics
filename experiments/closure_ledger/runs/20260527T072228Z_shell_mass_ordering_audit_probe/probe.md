# Shell Hamiltonian mass-ordering and `n_part` audit

**Run:** 2026-05-27T07:22:28+00:00

Tests whether the PR #77 shell waveguide basis reproduces the qualitative quark ordering structure better than the v3 lepton-shaped basis, and whether the `n_part = 233` compensator (PR #76) shrinks or disappears.

- **Identification**: shell basis structurally better than v3 (cavity wavefronts vs throat traversals); within-generation inversion requires sign-flipping χ_n (PR #79 slot); n_part NOT yet resolved at PR #78
- **Next PR**: PR #79 — boundary stress tensor to derive χ_n; PR #80 — color algebra to identify H_couple structure
- **B4 caveat**: shell ω dimensionful; mass-ratio audit scale-free; absolute scale via single B4 anchor (PR #53)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_shell_kinetic_spectrum_n_varied` | shell kinetic ω²(l=1, n=3,4,5) — range factor 2.2 | **PASS** |
| T2 | `T2_coverage_gap_shell_vs_observed` | coverage gap: shell ×2.2 vs observed ×6.4e9 (~9 OOM) | **PASS** |
| T3 | `T3_uniform_chi_cannot_invert` | uniform χ cannot invert (best 2/3 blocks) | **PASS** |
| T4 | `T4_signflipping_chi_n_existence` | sign-flipping χ_n CAN invert (existence proof) | **PASS** |
| T5 | `T5_kinetic_only_compensator_residual` | kinetic-only compensator log-residual (large) | **PASS** |
| T6 | `T6_compensator_with_signflipping_chi_n` | χ_n improves ordering, modestly improves spread | **PASS** |
| T7 | `T7_comparison_with_v3_lepton_shaped_baseline` | shell structurally better; n_part not resolved at #78 | **PASS** |
| T8 | `T8_honest_scope_assessment` | PR #78 SHARPENS PR #79–#80 scope; does not close | **PASS** |

## T2: Coverage gap — shell kinetic vs observed mass²

| range | min mass² | max mass² | factor |
|---|---:|---:|---:|
| shell kinetic ω² (n=3,4,5) | 14.63 | 32.49 | 2.22 |
| observed (u → t) | 4.67 | 2.98e+10 | 6.39e+09 |

Coverage deficit: **9.5 orders of magnitude in mass²**. Shell kinetic alone is insufficient; PR #79–#80 must contribute.

## T3: Uniform χ cannot reproduce within-generation inversion

| χ | l1n3 (u<d?) | l1n4 (c>s?) | l1n5 (t>b?) | correct/3 |
|---:|:---:|:---:|:---:|---:|
| -2.0 | ✓ | ✗ | ✗ | 1 |
| -1.0 | ✓ | ✗ | ✗ | 1 |
| -0.5 | ✓ | ✗ | ✗ | 1 |
| +0.0 | ✓ | ✗ | ✗ | 1 |
| +0.5 | ✗ | ✓ | ✓ | 2 |
| +1.0 | ✗ | ✓ | ✓ | 2 |
| +2.0 | ✗ | ✓ | ✓ | 2 |
| +5.0 | ✗ | ✓ | ✓ | 2 |

Best uniform-χ correctness: **2 of 3 blocks**. Uniform splitter is structurally insufficient — PR #79's boundary stress tensor must produce a generation-dependent χ_n.

## T4: Sign-flipping χ_n — existence proof

Illustrative `χ_n` (signs matter, magnitudes are not derived): `{'n=3': -1.0, 'n=4': 5.0, 'n=5': 20.0}`.

| species | shell mass² (illustrative) | within-block correct? |
|---|---:|:---:|
| u | 13.6278 | |
| d | 15.6278 | |
| s | 17.6659 | |
| c | 27.6659 | |
| b | 12.4902 | |
| t | 52.4902 | |

| within-block ordering | correct? |
|---|:---:|
| l1n3 (u<d) | ✓ |
| l1n4 (c>s) | ✓ |
| l1n5 (t>b) | ✓ |

All 3 within-block orderings correct: **True**. Global mass ordering: `['b', 'u', 'd', 's', 'c', 't']` (correct = `['u','d','s','c','b','t']`).

## T5–T6: Compensator measure (log₁₀ residual after best single-scale fit)

| configuration | log₁₀ residual |
|---|---:|
| kinetic only (no χ, no coupling) | 3.276 |
| kinetic + sign-flipping χ_n (illustrative) | 3.257 |
| observed range (for reference) | 9.806 |

Residual reduction from χ_n: 0.019 (modest). The remaining residual must be supplied by PR #79's full χ_n derivation and PR #80's H_couple inter-mode mixing.

## T7: Comparison with v3 baseline

**Shell basis improvements over v3:**

  - distinct from throat sector (cavity wavefronts)
  - kinetic = ω²(l, n), not phenomenological β·k²·(2π)
  - Z₂ partition slot for within-generation inversion (T4)
  - 3 × 2 = 6 flavors structural (PR #69 match)

**Not yet resolved at PR #78:**

  - coverage gap: shell-kinetic range factor 2.2 vs observed 6.4e9
  - within-generation inversion needs χ_n, not uniform χ (T3)
  - n_part compensator NOT reduced by PR #78 alone
  - PR #79 (boundary stress tensor χ_n) and PR #80 (color algebra coupling H_couple) must populate the operator

## Verdict

**SHELL_BASIS_STRUCTURALLY_BETTER_N_PART_NOT_YET_RESOLVED.** SHELL BASIS STRUCTURALLY BETTER, N_PART NOT YET RESOLVED. The PR #77 shell waveguide basis is structurally cleaner than the v3 lepton-shaped basis for the quark sector. The kinetic operator is ω²(l, n) — the cavity-wavefront eigenfrequency squared — not the phenomenological winding cost β·k²·(2π) that v3 inherited from the lepton ladder. The Z₂ partition slot is the natural structural home for the within-generation mass-ordering inversion (u < d but c > s, t > b).

STRUCTURAL FINDINGS. (1) A UNIFORM χ·σ_z partition splitter cannot reproduce the inversion: with any single sign of χ, at most 2 of 3 (l, n) blocks have the correct +/− ordering. (T3 confirms this with a parameter scan.) (2) A GENERATION-DEPENDENT, sign-flipping χ_n CAN reproduce the inversion — concretely, χ_3 < 0 (so + is lighter at n=3, u < d) with χ_4, χ_5 > 0 (so + is heavier at n=4, 5, c > s and t > b). T4 establishes the existence of this structural slot; the magnitudes are illustrative, NOT derived. PR #79 (boundary stress tensor on the cavity wall) is the natural channel to derive χ_n from physics.

COVERAGE GAP. The shell kinetic spectrum spans only a factor of ~2.2 in mass² (ω²(n=5)/ω²(n=3) ≈ 32.5/14.6) while the observed quark mass² spans a factor of ~6.4·10⁹ (u to t). Shell kinetic alone is therefore insufficient by ~9 orders of magnitude. Sign-flipping χ_n improves the ORDERING but not the SPREAD. Closing the spread requires either a large n-dependent χ_n (PR #79) or strong inter-mode coupling H_couple (PR #79 + PR #80).

N_PART STATUS. PR #78 does NOT reduce or remove the v3 phenomenological n_part = 233 compensator on its own. The shell basis identifies the structural SLOTS that must populate the operator, but the values that span the observed range are not derivable at PR #78 — they require PR #79 (boundary stress tensor → χ_n) and PR #80 (color algebra → H_couple). The n_part question is therefore SHARPENED, not closed: PR #78 identifies the right machinery and the missing structural inputs.

COMPARISON WITH V3. The shell basis (this PR) is structurally distinct from v3 in four ways: (i) basis = cavity wavefronts not throat traversals, (ii) kinetic = ω²(l, n) cavity eigenfrequency squared not phenomenological β·k²·(2π), (iii) Z₂ partition slot is the right structural home for the inversion, (iv) the 3 × 2 = 6 flavor count matches PR #69. These are qualitative improvements at the STRUCTURAL level. Quantitatively, PR #78 alone does not win on mass accuracy or n_part reduction — both bases require substantial compensation to span the observed range. The shell basis advantage is that its compensation is structurally located (PR #79 + PR #80 slots) rather than absorbed into a single phenomenological integer.

HONEST SCOPE. PR #78 SHARPENS what PR #79 must produce: sign-flipping χ_n from the boundary stress tensor, with magnitudes large enough to push the n=5 block up by a factor of ~10⁹ in mass². It also IDENTIFIES PR #80's role: H_couple inter-mode mixing for the inter-generation hierarchy (and the color algebra it transforms under). PR #78 does not close n_part; it relocates the phenomenological content to structurally-named slots that PRs #79–#80 must populate from physics.

## What this leaves open

- **PR #79** — derive the generation-dependent `χ_n` from a boundary stress tensor on the cavity wall. Signs and magnitudes are the structural input.
- **PR #80** — identify the color algebra acting on `(l, n, p)` and populate `H_couple` with inter-mode mixing terms transforming under that algebra. This is the channel for the inter-generation hierarchy.
- **`n_part` audit** — re-run after PR #79 and PR #80 populate the operator. If the resulting spread is structurally accounted for, `n_part = 233` (v3) is replaced by a derived value; otherwise it persists as a residual phenomenological compensator with reduced scope.

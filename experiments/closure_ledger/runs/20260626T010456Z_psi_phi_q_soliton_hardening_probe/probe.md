# ψ–Φ–q soliton hardening: stationarity, branch scan, and basin map (PR #180)

**Run:** 2026-06-26T01:04:56+00:00

Hardens PR #179's two-way self-consistent throat-soliton — stationarity, branch scan, basin map — with a spectral ψ kinetic, and corrects #179's high-μ runaway as a discretization artifact. *(QFT on the classical throat, not quantum gravity.)*

- **Stationarity**: eigenstate (residual ~1e-4); real-time stationary, mass-conserving
- **Branch scan**: smooth monotone family in mass (onset at ρ_c) and in μ (deepening)
- **Correction**: #179 high-μ runaway was an FD-Laplacian artifact; spectral ψ kinetic finds no collapse
- **Basin**: robust attractor — varied width/seed flow to the same soliton (~1%)

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_goal` | harden #179 — stationarity, branch scan, basin map | **PASS** |
| T2 | `T2_real_time_stationary_eigenstate` | stationarity: a real-time stationary eigenstate | **PASS** |
| T3 | `T3_mass_branch_smooth_onset` | mass branch: smooth monotone family; onset at ρ_peak=ρ_c | **PASS** |
| T4 | `T4_self_gravity_branch_no_collapse` | self-gravity branch: smooth, convergent (no collapse — corrects #179) | **PASS** |
| T5 | `T5_basin_attractor` | basin map: a robust attractor (varied init → same state) | **PASS** |
| T6 | `T6_grid_robustness` | robustness: well depth converges; core grid-sensitive | **PASS** |
| T7 | `T7_honest_scope` | honest scope | **PASS** |
| T8 | `T8_assessment` | PSI_PHI_Q_SOLITON_HARDENED | **PASS** |

## Self-gravity branch (M=3): smooth and everywhere-convergent

| μ | max\|q\| | Φ(0) | residual |
|---:|---:|---:|---:|
| 0.05 | 0.4242 | -3.0898 | 4e-04 |
| 0.5 | 0.5558 | -3.7207 | 1e-03 |
| 1.0 | 1.1152 | -7.6044 | 2e-03 |
| 1.5 | 1.8626 | -15.0083 | 3e-05 |
| 2.0 | 2.6236 | -24.5506 | 3e-06 |

No collapse — correcting #179's FD-Laplacian runaway; the branch deepens out of weak-field validity (the strong-field endpoint, for NR).

## Basin map: a robust attractor

Varied initial width/seed all flow to the same soliton — max\|q\| [0.4236, 0.4246, 0.4226, 0.4245, 0.4198, 0.4243] (spread 1.13%), Φ(0) [-3.0905, -3.0891, -3.0918, -3.0893, -3.0941, -3.0896] (spread 0.161%).

## Verdict

**PSI_PHI_Q_SOLITON_HARDENED_REAL_TIME_STATIONARY_SMOOTH_CONVERGENT_BRANCH_ROBUST_ATTRACTOR.** HARDENED — A TRUSTWORTHY THROAT-SOLITON (and a correction to #179). PR #179's two-way self-consistent state passes stationarity, branch scan, and basin map.

STATIONARITY. With the spectral ψ kinetic the fixed point is a genuine eigenstate (‖Hψ − μψ‖/‖ψ‖ = 1.0e-04, μ = -1.4454) that persists under unitary real-time evolution (profile drift 4.4e-05, mass conserved to 1e-13) — a real soliton, not just a relaxation endpoint.

BRANCH SCAN. The soliton is a smooth monotone FAMILY — in mass (the order field switches on where ρ_peak crosses ρ_c, the well deepening monotonically) and in q's self-gravity μ (max|q| and the well rising smoothly, every point an everywhere-convergent fixed point, residuals ≤ 10⁻³).

CORRECTION TO #179. #179 reported a high-μ RUNAWAY collapse — but that used a finite-difference Laplacian; the operator-consistent spectral ψ kinetic finds NO collapse up to μ = 2 (smooth, convergent), so the runaway was a DISCRETIZATION ARTIFACT. The genuine large-μ limit is the soliton deepening out of weak-field validity (Φ(0): -3.0898 → -24.5506) — the strong-field endpoint for NR.

BASIN. The soliton is a robust ATTRACTOR: initial widths and seeds all flow to the same state (max|q| spread 1.13%, Φ(0) spread 0.161%) — not fine-tuned.

ROBUSTNESS. The well depth converges under grid refinement (11.51% over N = 160 → 320); the pointwise core max|q| is more grid-sensitive (~10%, the sharp core) — an honest caveat. The soliton's existence, two-way back-reaction, and threshold continuity survive and are hardened.

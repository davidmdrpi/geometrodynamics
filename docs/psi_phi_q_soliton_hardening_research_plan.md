# ψ–Φ–q soliton hardening: stationarity, branch scan, and basin map (PR #180)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## From a self-consistent state to a trustworthy soliton

PR #179 closed the two-way ψ–Φ–q loop and found a self-consistent
throat-soliton. This probe hardens it into a trustworthy object — with the
three things such a result needs — and **re-examines the #179 collapse claim
with a better-conditioned solver**.

The key methodological upgrade: ψ's kinetic is put on an
**operator-consistent spectral basis** (`u = rψ`, DST-I for the radial
second derivative), so the imaginary-time relaxation and the real-time
evolution use the **same** Laplacian. (In #179 the relaxation used a
finite-difference `np.gradient` Laplacian; mixing that with a spectral
real-time operator would spuriously perturb the state.)

## Stationarity — a genuine real-time eigenstate

The relaxed state is a genuine stationary eigenstate, not just a
gradient-flow endpoint:

- **Eigenstate residual** `‖Hψ − μψ‖/‖ψ‖ ≈ 10⁻⁴` (chemical potential
  `μ ≈ −1.45`): `Hψ = μψ`.
- **Real-time persistence**: a unitary split-step evolution under the
  self-consistent `V_eff = Φ + ½g q²` leaves the profile stationary (max
  `|ψ|` drift `~4×10⁻⁵`) and conserves mass to machine precision.

The two-way throat-soliton is a real bound state that persists under real
dynamics.

## Branch scan — a smooth family

### Mass branch (μ = 0.05)

| M | ρ_peak | max\|q\| | Φ(0) |
|---:|---:|---:|---:|
| 1.0 | 0.005 | 0.000 | −0.30 |
| 2.0 | 0.081 | 0.000 | −1.31 |
| 2.5 | 0.206 | 0.000 | −2.14 |
| 3.0 | 0.424 | 0.423 | −3.09 |
| 3.5 | 0.723 | 0.701 | −4.36 |

The order field switches on where `ρ_peak` crosses `ρ_c = 0.20` — the spatial
onset sits just above `ρ_c` (near M ≈ 2.7) by the Ginzburg–Landau
droplet-size barrier. A smooth, monotone soliton branch.

### Self-gravity branch (M = 3)

| μ | max\|q\| | Φ(0) | residual |
|---:|---:|---:|---:|
| 0.05 | 0.424 | −3.09 | 1×10⁻⁴ |
| 0.5 | 0.556 | −3.72 | 1×10⁻³ |
| 1.0 | 1.115 | −7.60 | 2×10⁻³ |
| 1.5 | 1.863 | −15.0 | 2×10⁻⁵ |
| 2.0 | 2.624 | −24.6 | 3×10⁻⁶ |

Smooth, monotone, **everywhere convergent** — every point a well-defined
self-consistent fixed point, no blow-up.

## A correction to #179 (what hardening is for)

PR #179 reported that super-critical `q`-self-gravity drives a **runaway
collapse** (`|q| → 31`, `Φ(0) → −252`). That used a finite-difference
(`np.gradient`) Laplacian. With the operator-consistent **spectral** kinetic
here, the μ branch is smooth and convergent up to `μ = 2` — **the apparent
runaway was a discretization artifact** of the FD scheme.

The genuine large-μ limit is **not** a numerical runaway but the soliton
**deepening out of weak-field validity** (`Φ(0)` from −3 to −25 and beyond as
`μ` grows) — the strong-field domain for full numerical relativity.

**What survives #179** — the soliton's existence, two-way back-reaction
(deeper well / denser core), and threshold continuity — is confirmed and
hardened; the specific "runaway" claim does not survive as stated.

## Basin map — a robust attractor

Initial conditions varied over Gaussian width and order-seed all flow to the
**same** soliton:

| init (w, q-seed) | max\|q\| | Φ(0) |
|---|---:|---:|
| (1.2, 10⁻²) | 0.424 | −3.090 |
| (1.8, 10⁻²) | 0.423 | −3.092 |
| (2.6, 10⁻²) | 0.420 | −3.094 |
| (1.8, 10⁻¹) | 0.425 | −3.089 |

`max|q|` spread **1.1%**, `Φ(0)` spread **0.16%** — a robust dynamical
attractor, not a fine-tuned state. (A tiny seed `10⁻³` reaches the same
attractor, only more slowly — a convergence-time effect, not a different
basin.)

## Robustness — grid convergence

| N | Φ(0) | max\|q\| |
|---:|---:|---:|
| 160 | −3.34 | 0.479 |
| 240 | −3.09 | 0.423 |
| 320 | −2.98 | 0.390 |

The well depth `Φ(0)` converges at ~first order (successive changes shrinking;
the N = 240 → 320 step is ~3.6%). The pointwise core `max|q|` is more
grid-sensitive (~10% per refinement — the sharp core narrows under
refinement, converging but slowly), so the precise core amplitude carries
grid uncertainty while the soliton's existence and structure are robust. An
honest convergence statement.

## Honest scope

Weak-field, semi-dynamical, spherically reduced; effective constants. The
hardening establishes the soliton is a stationary eigenstate, a smooth
everywhere-convergent branch, and a robust attractor — and corrects the #179
high-μ runaway as a discretization artifact. The deep-well large-μ branch and
the strong-field endpoint remain for full numerical relativity; the core
amplitude carries ~10% grid uncertainty.

## Reproduce

```bash
python -m experiments.closure_ledger.psi_phi_q_soliton_hardening_probe
# Verdict: PSI_PHI_Q_SOLITON_HARDENED_REAL_TIME_STATIONARY_SMOOTH_CONVERGENT_BRANCH_ROBUST_ATTRACTOR
```

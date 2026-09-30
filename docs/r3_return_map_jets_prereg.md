# Specification: constrained R3 return map and the leading frequency-shift coefficient

Date: 2026-09-29. Parent: `claude/geometrodynamics-qft-audit-vpktax` at `40d03da`
(the #315 review is open against it).
Publish this specification before any cubic-order computation for the ESU.
The measurement code is committed together with this document. Its hashes
are bound into the archive. Corrections are dated notes only.

## 1. Question

The review in #315 showed that the kicked-trajectory run (freeze `2e984ac`)
cannot decide R3. That run had three problems: restart kicks with no
error budget, finite-window estimator bias, and signed displacement used
as if it were action. This test replaces it with a direct local computation.

On the constrained clock section, the n=2 LRS tensor rotates at
frac(rho(I)) = theta(I)/2π, where theta(I) is in (0, π) and I is the action
of the reduced symplectic structure. The quantity measured is

    theta(I)/2π = theta0/2π + nu I + O(I^2),

with theta0/2π = .484666. The #315 linear control fixed the branch
rho0 = 1.484666. The refocusing target frac(rho) = 1/2 corresponds to
theta = π, so **nu > 0 is a shift toward the target and nu < 0 a shift away
from it**. Both labels are local statements, as #315 required. A sign alone
says nothing about a crossing, and nothing about closure.

## 2. Reduced system, section and symplectic structure

LRS restriction of `nonlinear_supported_tt` (kappa = a = 1, conformal time):
beta = x b0, with b0 = diag(1,1,-2)/sqrt6, and q = (q,0,0,0). The Lagrangian

    L = -3A'^2 + q'^2/2 + H x'^2/2 - V,   H = A^2 - q^2/6,
    V = -H r/2 + q^2 tr M^-1 /2 + 3A^4/2,  r = 2(2 tr M^-1 - tr M^2),

gives Euler–Lagrange equations identical to `conformal_rhs`. Its energy
T + V equals the Hamiltonian constraint. A pre-freeze test confirms both to
1e-13 at random states (`tests/test_r3_return_map.py`). The momenta are
p_A = -6A', p_q = q' and p_x = H x'.

- **Section:** q = 0 with q' < 0, with the constraint solved for q'.
- **Reduced coordinates:** z = (A, p_A, x, p_x), canonical for
  omega = dp_A ∧ dA + dp_x ∧ dx. On the section dq = 0, so the p_q term drops out.
- **Return map:** P(z) is the next section crossing, at conformal time
  roughly π. The return time depends on z, and both methods include that
  dependence exactly.
- **Fixed point:** z* = (1, 0, 0, 0), the breathing ESU.
- **Action:** I = (1/2π) ∮ (p_A dA + p_x dx) over an invariant circle. It is
  not the signed turning displacement.

## 3. Method 1: jets and Birkhoff normal form

`geometrodynamics/waves/jets.py` and `r3_return_map.jet_return_map`:

1. Expand the section coordinates as order-3 jets about z*, and solve the
   constraint for q' as a jet.
2. Integrate the reduced scalar equations with fixed-step RK4 at h = π/N.
   The last 16 steps use a jet-valued step s, found by Newton so that q = 0.
   This is the return-time correction, exact through third order.
3. Run N = 1024, 2048 and 4096. Primary: Richardson extrapolation of
   (2048, 4096). Convergence check: Richardson of (1024, 2048).
4. Normal form (`normal_form`):
   - diagonalise the linear part into u (multiplier 85.02), s (.01176),
     w (mu = e^{i theta0}) and wbar;
   - remove every quadratic term. All are nonresonant, and the removal
     includes the centre-manifold corrections in u and s;
   - read the resonant coefficient g of w^2 wbar in the w component.
   - Then theta(I) = theta0 + Im(g/mu)|w|^2, and I = |omega(q, qbar)| |w|^2 + O(|w|^4), so

         nu_M1 = Im(g/mu) / (2π |omega(q, qbar)|).

   - Re(g/mu) must vanish for a symplectic map. That is a check, not an
     assumption.

On a toy map, a hyperbolic block times an exact twist conjugated by nonlinear
symplectic shears, this procedure returns nu = beta/π to 1e-10. It does so at
theta0 = 2.5 and at 3.0, which is nearer the 1:2 resonance.

## 4. Method 2: invariant circles on the full system

`r3_return_map.esu_map` and `invariant_circle`:

- **Map evaluation.** P is evaluated with the full matrix system
  (`conformal_rhs`), DOP853 at rtol 1e-12, atol 1e-14, with event location.
  None of Method 1's equations, integrators or jets are used.
- **Circle solve.** Solve P(K(θ)) = K(θ + omega) for a Fourier circle K on
  31 grid points. The x-harmonic is fixed at c1 = a/2 (real). Use
  Gauss–Newton with finite-difference Jacobians (h = 1e-6) and at most 15
  iterations, stopping at a residual of 1e-12.
- **Amplitude ladder.** a ∈ {.004, .008, .016, .032}. The first circle starts
  from the linear circle given by Method 2's own finite-difference Jacobian.
  Each later circle starts from the previous one, rescaled.
- **Measured quantities.** omega and I, the latter spectrally.
- **Estimate.** With theta0 = arccos(tr M_T2(π)/2) from #310, form
  y_k = (omega_k - theta0)/(2π I_k). nu_M2 is the I = 0 intercept of a
  quadratic fit over all four points. Its uncertainty is the difference from
  a linear fit over the three smallest.
- **Integrator check.** At a = .008, re-evaluate the converged circle with
  Radau.

On the toy map the circles reproduce the exact rotation to 1e-12.

## 5. Checks

The archive records every check. Failing any one of them gives **UNRESOLVED**.

| id | requirement |
|---|---|
| C1 | fixed-point error <= 1e-10 (both methods) |
| C2 | \|theta_M1 - theta0\| <= 1e-9; \|theta_M2(FD) - theta0\| <= 1e-6 |
| C3 | Method 1 jet symplecticity (DP^T Ω DP - Ω through order 2) <= 1e-10 × max coefficient; FD Jacobian defect <= 1e-4 |
| C4 | \|Re(g/mu)\| <= 1e-5 \|Im(g/mu)\| |
| C5 | Method 1 image-constraint jet and section-q jet <= 1e-8 × max coefficient; Method 2 constraint residual <= 1e-10 at the fixed point and at every circle image |
| C6 | \|mu^k - 1\| >= .1 for k = 1..4; every quadratic divisor >= .1 |
| C7 | \|nu_M1(2048,4096) - nu_M1(1024,2048)\| <= 1e-6 \|nu_M1\| + 1e-9 |
| D1 | circle invariance residual <= 1e-10 at every amplitude |
| D2 | Fourier tail (\|k\| >= 13) <= 1e-10 |
| D3 | Radau invariance residual of the a = .008 circle <= 1e-9 |
| A1 | uncertainty_M2 <= 1e-3 \|nu_M1\|, and \|nu_M1 - nu_M2\| <= max(1e-4 \|nu_M1\|, 3 uncertainty_M2) |
| A2 | the remainder omega_k - theta0 - 2π nu_M1 I_k has log-slope >= 1.8 in I between consecutive amplitudes with \|remainder\| >= 1e-8 (predicted O(I^2)); at least one such pair |

A pre-freeze run of the finite-difference Jacobian at z* gave a symplectic
defect of 3.6e-6. The C3 bound of 1e-4 was chosen for that reason, before any
nu-bearing number existed.

## 6. Decision and continuation

- **Resolved** means every check passes and
  |nu_M1| >= 10 max(C7 difference, |nu_M1 - nu_M2|, uncertainty_M2).
- **SHIFT_AWAY_FROM_TARGET** (resolved, nu < 0): stop the small-amplitude R3
  search. Keep the narrow conclusion, namely that the LRS n=2 tensor's
  rotation moves away from the refocusing resonance near the breathing
  ESU. This does not exclude a large-amplitude turn, other sectors, or
  other closure conditions.
- **SHIFT_TOWARD_TARGET** (resolved, nu > 0): no crossing is claimed. A
  separate freeze would search for a crossing with multiple shooting or
  collocation, then test the full-state condition P^2(z) = z.
- **UNRESOLVED:** improve the error budget. No physical verdict is assigned.

A verified periodic orbit would still establish only a classical closed
history. It would not show that action is selected or quantized. That needs
a physical reason for the closure condition, and evidence about the family
of actions that results.

## 7. Assumptions, and how each is checked

| Assumption | How it is checked |
|---|---|
| P is smooth near z* | H, M > 0; q' = -sqrt3 ≠ 0 at the section, so the crossing is transversal (implicit function theorem for the return time) |
| z* is a saddle-centre fixed point | C1, C2, and the multipliers 85.02, .01176, e^{±i theta0}; no multiplier is 1 on the section |
| Nonresonance through order 4 | C6. The smallest divisor, \|mu^2 - 1\| ≈ .19, comes from the nearby 1:2 resonance |
| Symplecticity | C3, C4 |
| Nonzero twist | the "resolved" condition |
| Invariant circles exist at the measured amplitudes (KAM, Diophantine rotation) | assumed, not proved; supported numerically by D1, D2 and Newton convergence |

A reference suggested for saddle-centre normal forms (arXiv:1501.05935)
could not be retrieved from this environment and was not consulted. The
method here is the textbook order-3 Birkhoff reduction. It is validated on a
toy map with an exact answer, and on this system through C3–C7.

## 8. Pre-freeze disclosure

Before this freeze, only the following were run:

- toy-map tests of both methods;
- the reduced-versus-full equation check;
- ESU jets of order 1 and 2 (multipliers, the #310 trace to 5e-12, and the
  quadratic symplectic defect after Richardson, about 4e-10 on coefficients
  up to 1.6e4);
- the Method 2 fixed point and finite-difference Jacobian at z*.

No order-3 ESU jet, normal form or ESU invariant circle has been computed.
The kicked run's descriptive coefficient (#315) is known. It is not used by
any rule here.

Crashes in the committed code may be fixed after the freeze, under a dated
note. Any such fix must not change parameters, tolerances or rules.

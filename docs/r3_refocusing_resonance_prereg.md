# Preregistration: does nonlinearity carry the n=2 tensor mode onto the refocusing resonance?

Date: 2026-09-29. Parent: `main` at `3c8f0d7`.
Publish this freeze before implementing or running any amplitude-dependent
measurement below. Retain failed outcomes. Corrections are dated notes;
nothing here is edited after publication.

## 1. Why this test, and why it comes first

Step 4 (#313, #314) asks whether the dynamics can *select* discrete action
values. Well-posed evolution makes every finite-time readout a continuous
function of the preparation. #314 closed the plateau route (R1) in the
homogeneous sector: the response is approximately quadratic, with no
plateau. Two routes remained:

- **R2, discrete events.** The correct integer is the degree N of phi/|phi|
  on the closed slice S3 (pi_3(S3)=Z, no window boundary). With the
  background's oddness phi(-x)=-phi(x), Borsuk–Ulam makes N odd, and zeros
  occur in antipodal pairs with Delta N = ±2. A scratch check before this
  freeze found:
  - For a random odd cubic perturbation of the breathing background, N ran
    +1, -1, +1, -1, +1 through the collapse at eta=pi/4, and was never even.
  - The event times, divided by epsilon, were identical at epsilon =
    .001, .01 and .05.
  - A complex-structure (Hopf) perturbation produced no events.

  Absorbed action is continuous in the preparation while N jumps, so the
  relation "action = q·N" forces q=0. R2 therefore yields a conserved Z2
  (N mod 2) but no action quantum, and it is not pursued further. The
  scratch computation is motivation only. No verdict here depends on it.
- **R3, closed histories.** In a Hamiltonian system, histories that close
  on themselves after a fixed winding are isolated in amplitude. That is the
  one remaining mechanism for exact discreteness. This preregistration
  tests its cheapest necessary condition.

## 2. The registered target, and a correction to the proposal

The proposal to the user said to impose phi(eta+pi,x) = -phi(eta,-x) on
the full nonlinear system. For the scalar quartet, that sign flip is a
symmetry: O(4) together with the antipodal isometry. **For the metric it is
not.** Free linear gravitational waves refocus as h -> -h after conformal
time pi. On the n=2 tensor block of #310, Rf = P M(pi) = -I. But g-h and
g+h are different metrics, and the Einstein equations have no symmetry
h -> -h; the Bianchi IX potential contains an odd tr(beta^3) term. A twisted
period-pi identification is therefore not a consistent nonlinear boundary
condition for the tensor sector.

The consistent formulation replaces it with a closure condition, which is
reparametrisation-invariant:

> **Refocusing resonance.** Let rho be the tensor phase advance per
> oscillation period of the scalar clock, in units of 2 pi. A free conformal
> field has rho = omega_T/omega_q = 3/2 for n=2. At rho = 3/2 the tensor is
> inverted after each clock period, and the joint history closes after two
> clock periods. R3's candidate discrete histories are the nonlinear orbits
> with rho = 3/2 exactly.

Linearly, rho sits off this value. #310 gives tr M_T2(pi) = -1.99073, so
frac(rho0) = .48467 or .51533, depending on branch. WKB,
omega_T^2 = 8 + <2R^2/f> = 8.8285, suggests rho0 ≈ 1.4857 < 3/2. The
frozen code fixes the branch by continuous phase unwrapping and does not
rely on WKB. Closed rho = 3/2 histories exist at small amplitude only if
the nonlinear shift of rho has the right sign:

    rho(eps) = rho0 + c eps^2 + O(eps^4),   require sign(c) = sign(3/2 - rho0).

**That sign is the gate.** If it is wrong, the resonance is not reached
perturbatively, and R3's refocusing branch is closed in the only tensor
sector that can be computed exactly here.

## 3. System: the exact homogeneous n=2 tensor sector

Use `geometrodynamics/waves/nonlinear_supported_tt.py` (`conformal_rhs`,
`constraints`) unchanged, with kappa = a = 1. This is the exact homogeneous
Einstein–quartet system on S3: scale A, quartet q in R^4 acting as
phi = B(q)x, and unimodular left-invariant shape M = exp(2 beta).

- **Linear check.** Its shape equation linearises to
  b'' + (f'/f) b' + [8 + 2R^2/f] b = 0, with f = A^2 - |q|^2/6. That is the
  #310 tensor equation at n=2. The left-invariant anisotropies are n=2 TT
  harmonics.
- **Invariant subsystem.** Restrict to the locally rotationally symmetric
  (LRS, Taub) subsystem, which the equations preserve exactly:
  - beta = x·b0, with b0 = diag(1,1,-2)/sqrt6 (tr b0^2 = 1), so L = x'·b0;
  - q = (q0,0,0,0), so the quaternionic current vanishes and the momentum
    constraint holds identically.
- **Background.** A = 1, q0 = (sqrt3/2) cos 2 eta, x = 0. This is the
  breathing ESU, with period pi.
- **Linear facts computed before the freeze** (script: `experiments/closure_ledger/r3_prefreeze/linear_facts.py`):
  - The background returns after pi to 3e-13.
  - The homogeneous (A, A', q0, q0') block has multipliers 85.020, .011762,
    1 and 1 per pi. This is the Einstein-static instability plus the time
    translation and constraint directions.

## 4. Measurement

1. **Clock section.** Take q0 = 0 with q0' < 0, once per clock period. At the
   section, Q = 0 and H = A^2.
2. **Initial data at the section.** Set x = s_pol·eps with s_pol = ±1, and
   x' = 0, M, L from x, A' = 0, and A = 1 + s. Solve q0' from the Hamiltonian
   constraint on the negative root:
       q0'^2 = 2[3A'^2 - H ell/2 + H r/2 - 3A^4/2].
3. **Centre-manifold tracking** (straddle method). The Einstein-static
   direction grows 85-fold per period, so orbits are kept on the
   centre-stable manifold by repeated bisection.
   - Starting from the current section state, bisect a kick s to A, then
     re-solve q0' from the constraint.
   - Classify each trial by its runaway over at most 9 sections:
     - |A - 1| > 0.5 counts as expanding (+1) or collapsing (-1).
     - A chart failure counts as -1.
     - If there is no exit, classify by sign(A - 1) at the ninth section.
   - Bisect until the bracket is narrower than 4e-16 or 70 iterations have
     run. Accept the first 4 sections of the midpoint trajectory, then
     restart from the last accepted state.
   - Initial brackets are ±1e-2 for the first window and ±1e-6 after that.
   - Collect K = 48 accepted clock periods, discarding none.
4. **Tensor phase.** Use Theta(eta) = unwrap(atan2(-x'/3, x)), sampled at
   256 dense-output points per period. Let dTheta_k be its increment over
   clock period k. The rotation number is the weighted Birkhoff average
       rho_K = sum_k w_k dTheta_k / (2 pi sum_k w_k),
   with w_k = exp(-1/(t_k(1-t_k))) and t_k = (k+1/2)/K. It does not depend
   on how the angle is defined, and it converges faster than any power of K
   on smooth tori.
5. **Integrators.**
   - Primary: DOP853 with rtol 1e-12, atol 1e-14.
   - Secondary: DOP853 with rtol 1e-10, atol 1e-12, re-run independently,
     bisection included.
6. **Amplitude ladder.** eps in {.01, .02, .04, .08, .16}, times s_pol = ±1:
   10 runs per integrator. There is no further search over eps, K,
   tolerances or window lengths.
7. **Linear reference rho0.** From #310, 2cos(2 pi a) = tr M_T2(pi) gives
   `fl.monodromy('T', 2)`. rho0 is whichever of 1+a and 2-a lies nearer the
   primary rho at eps = .01, s_pol = +1.

## 5. Gates

Numerical gates. Failing any of them makes the verdict **UNRESOLVED**, with
the failure reported.

- **N1 (constraints).** Every accepted section satisfies
  |Hamiltonian residual| <= 1e-8.
- **N2 (orbit shadowing).** Every accepted state stays within |A-1| <= .1.
  The phase-plane radius sqrt(x^2 + (x'/3)^2) stays between .2 and 5 times
  eps.
- **N3 (convergence).** Define
  err(eps) = max(|rho_K - rho_{K/2}|, |rho_primary - rho_secondary|),
  where rho_{K/2} uses the first 24 periods. Then for every run:
  - Delta(eps) = rho(eps) - rho0 is resolved when |Delta| >= 100·err;
  - err <= 1e-5.
- **N4 (linear control).** At eps = .01, |rho - rho0| <= 1e-3.
- **N5 (scaling).** For each s_pol, the log-slope of |Delta| between
  consecutive resolved eps in {.01, .02, .04} lies in [1.7, 2.3]. At least
  one such pair must exist for each s_pol.

Decision gate, applied only if N1–N5 pass:

- Let c(eps) = Delta(eps)/eps^2 over resolved runs with eps <= .04 and both
  s_pol values. All must share one sign, sign(c). Otherwise the verdict is
  **UNRESOLVED**.
- **PASS** if sign(c) = sign(3/2 - rho0). The resonance is reached at finite
  amplitude. Report eps* ≈ sqrt((3/2 - rho0)/c) from the eps = .01
  coefficient.
- **FAIL** if sign(c) = -sign(3/2 - rho0). The refocusing resonance is not
  reached perturbatively in the n=2 LRS sector, and R3's refocusing branch
  is recorded as closed there.

The runs at eps = .08 and .16 are reported descriptively, including whether
rho crosses 3/2 anywhere in the ladder. They do not change the verdict.

## 6. What a verdict does not establish

- **PASS** does not establish selection. Three further steps would be
  needed, none implied by this test:
  1. Locate the closed rho = 3/2 orbit and its action.
  2. Show that the result survives other polarisations and the
     inhomogeneous n >= 3 sectors.
  3. Supply a physical reason why histories must close. The antipodal-return
     premise of the program is an assumption, not a derived result.
- **FAIL** is confined to the LRS n=2 tensor at small amplitude. It says
  nothing about large amplitude or other sectors, except that this is the
  only exactly computable one.
- The Einstein-static instability means the tracked orbits exist only on a
  centre manifold. Generic data run away. Any selection built on them would
  also require explaining how physical data reach that manifold.

## 7. Outputs

- `geometrodynamics/waves/r3_resonance.py`
- `experiments/closure_ledger/r3_resonance_probe.py`
- `experiments/closure_ledger/runs/20260929_r3_resonance/r3_resonance.json`,
  which binds the SHA-256 of both source files
- `tests/test_r3_resonance.py`
- `docs/r3_refocusing_resonance.md`, which reports whichever verdict results

## 8. Pre-freeze disclosure

Before this freeze, only the following were computed:

1. The background return error and homogeneous-block multipliers quoted in
   section 3, at eps = 0.
2. tr M_T2(pi) and the WKB mean, both already archived in #310.
3. The R2 scratch computation summarised in section 1 (`r3_prefreeze/r2_degree_*.py`).

No amplitude-dependent frequency, rotation number or nonlinear orbit of the
tensor sector has been computed.

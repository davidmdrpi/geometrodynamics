# Antipodal shift with normal reversal: pre-registration

Baseline: `d77d915bcab363dc4e75c31987dcdc4d93910302` (merged PR #352). This file is published before
the implementation and measurements. Analytic predictions below are derived
in advance; exploratory outcomes are not promoted to predictions afterwards.
Seed: 2026111201.

## Question and scope

Can a spatially non-orientable throat preserve angular circulation while
its two mouths respond? Compare the antipodal normal reversal
A(phi,theta,u)=(phi+pi,theta+pi,-u) with the reflective control
R(phi,theta,u)=(phi+pi,-theta,u). Time is preserved in both.

The A surface quotient is a torus, not a Klein bottle. Its normal line
reverses around a quotient loop. R reverses surface orientation. Gluing
monodromy is not connection holonomy without a transport prescription.
Periods below are coordinate periods, not metric lengths.

Use a real field
F(phi,theta,u,t)=a(phi,u,t) cos(n theta)+b(phi,u,t) sin(n theta).
The two quadratures are from one real angular harmonic, not an added
charged matter field. This is a finite-dimensional Galerkin/normal-lattice
model, not a solved Einstein constraint or Maxwell boundary problem.
Neither time reversal nor an antilinear internal lift is inserted.

## Model fixed before running

Periodic phi grid N=16 (primary), 24 (replicate); radius 1. Two normal
samples u=+a,-a, joined with stiffness kappa=0.75. Angular harmonics n=1,2.
Mass squared 0.25. Symmetric mouths have waist coordinates q_+,q_-,
mass M=1 and support stiffness K_s=4. All field coordinates use the
orthonormal discrete phi basis and normalized angular basis.

H = 1/2 sum(p_a^2+p_b^2)
  + 1/2 sum_c c^T [L_phi+n^2+mu^2+kappa L_normal+V(q)] c
  + sum_s [p_qs^2/(2M)+K_s q_s^2/2],

V_s(phi)=2 exp[-2 h(phi) q_s], h(phi)=0.6+0.4 cos(2 phi).
The common pi-periodic profile makes H compatible with both controls.
q_s is an outward-opening scalar, so A swaps q_+ and q_- without a sign.
The signed outward normal displacements are q_+ and -q_-; distinguish
these from independent translational recoil degrees of freedom.

The response force is derived from this H:
M qddot_s=-K_s q_s+sum_phi h V_s (a_s^2+b_s^2).
No prescribed torque, gain, damping or advanced signal is used.

## Predictions

P1. Jacobians: det(A_surface)=+1, det(A_throat)=-1;
det(R_surface)=det(R_throat)=-1. A^2=R^2=identity on the cover torus.
A_* d_theta=+d_theta and R_* d_theta=-d_theta.
Uniform nonzero theta flow descends only for A. R can still carry a
nowhere-zero full vector field (e.g. d_phi); zeros of its angular component
are not compulsory stagnation points.

P2. For F(Ax)=eta F(x), eta=+1 or -1, a Fourier mode with phi index m and
normal parity epsilon_u is allowed iff
(-1)^(m+n) epsilon_u=eta.
Verify both eta sectors, harmonics n=1,2, all resolved real Fourier modes,
and compare a separately assembled grid deck with this rule. Frozen q=0
frequencies squared:
4 sin^2(pi m/N)/(2pi/N)^2+n^2+mu^2+2+kappa(1-epsilon_u).
Frequency convergence is tested against the continuum m^2 formula at
N=16,24,32,48 for fixed m=1,2,3; it is an approximation, not an exact continuum.

P3. For the canonical angular momentum
L_theta=n sum(a p_b-b p_a),
A preserves L_theta, R reverses it, and both are canonical (time kept).
Thus A allows nonzero L_theta on its fixed space, whereas R-fixed data
have zero total L_theta. Antiperiodic eta=-1 does not change that statement.
This is angular momentum of a real scalar mode, not electric charge or
spin-1/2. Normal orientation alone does not make conjugate electric charges.

P4. On A-fixed data the dynamically responding mouths stay equal in their
outward-opening coordinates and the per-face angular momenta are equal.
Signed normal loads are opposite. Test the seam defect and global
L_theta conservation (<1e-9 relative, using nonzero-signal runs);
the equality of face quantities is a separate check, not inferred from a
cancellation in the total. Normal loads are generalized waist loads, not
an ADM or bulk momentum-flux balance.

P5. The frozen two-face stiffness for a mode is
Z_m(Omega)=[[lambda_m+kappa-Omega^2,-kappa],
            [-kappa,lambda_m+kappa-Omega^2]],
lambda_m=4 sin^2(pi m/N)/(2pi/N)^2+n^2+mu^2+2.
The deck selects its symmetric or antisymmetric eigenchannel according
to P2. Compare the projected inverse to an independent full-grid solve at
Omega=0.3,1.4 for m=0,1, n=1,2 and both eta.
Changing kappa shifts only the odd-normal channel. These are susceptibilities
of this normal lattice, not a continuum throat Dirichlet-to-Neumann map.

P6. Mechanical response is quadratic in small field amplitude; the field
correction to a calculation with q frozen at zero is cubic. Measure slopes
over amplitudes 0.05,0.1,0.2,0.4 on [0,4]; accept slopes 2 +/-0.15 and
3 +/-0.2. Report failure rather than retuning the range.

## Runs and numerical checks

For each N and n use eta=+1 and m=0,1, normal parity dictated by P2.
Start a=A g, b=p_a=0, p_b=omega A g, ||g||=1, with q=p_q=0.
These are exact travelling modes only for the frozen Hamiltonian.
Amplitudes 0.05,0.1,0.2,0.4; integrate to t=40 using velocity Verlet
with h=0.01, recording every 10 steps. Report extrema over every draw/run,
not just medians. Required seam and Noether checks: <1e-9 relative;
maximum relative energy excursion <2e-3.

Numerical validation: highest-amplitude primary n=1,m=0 trajectory;
compare h=0.01,0.005 and independent DOP853 (rtol=1e-10, atol=1e-12)
on [0,40]. Require relative final-state difference <1e-2 against DOP853,
and reduction under halving the step. Compare the field gradient, mouth
force and canonical vector field by finite differences in tests.

## Measurements, not predictions

Mouth excursions, deformation of the circulating pattern, whether each
face current retains its sign, and response amplitudes in the two normal
channels. A conserved total angular momentum does not establish uniform
pointwise angular speed, one-way transport at every site, stability for
arbitrary time, or spontaneous generation of circulation.

## Interpretation and verdicts

Report separate geometry, spectral, angular momentum, responding-mouth,
susceptibility and numerical-validity checks. Missing checks are INCOMPLETE;
numerically failed checks are UNRESOLVED_NUMERICALLY. The orienting wave,
support, normal discretization and metric are assumed. Electric charge,
spin lift, time extension, Born law and Einstein dynamics remain NOT_DERIVED.
A successful result would establish compatibility of normal reversal,
circulation and mechanical response in this model, not a derivation of QFT.

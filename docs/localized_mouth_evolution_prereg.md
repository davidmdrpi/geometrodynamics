# Prospective short-time handle evolution and crossing experiment

Date: 2026-09-23. Parent: #307, `19a674d8625469c869c4f373efe062c7ab55c1e1`.
Publish before implementing or measuring this evolution. Retain failed outcomes.

## Questions and scope

Evolve the certified four-scalar initial data with gravitational backreaction.
Does a timelike test worldline cross the evolved neck? Does the covariant
matter momentum balance hold? Does a finite momentum impulse persist as the
window around a crossing shrinks? These are separate questions. No detector,
reset, projection, kick, packet-counting rule or quantization is introduced.

Test geodesics have negligible stress and cannot establish reciprocal exchange
with a self-gravitating absorber. A pass for test-worldline crossing or matter
balance is not a pass for a discrete reciprocal measurement event. The
experiment can instead give a resolved continuous counterexample to a claim
that every smooth neck crossing automatically supplies a momentum jump.

## Full equivariant system, analytic predictions

Einstein frame, unit lapse, zero shift:
gE=-dt^2+A^2 ds^2+B^2 dOmega^2; mixed extrinsic curvature diag(k,l,l),
K=k+2l, with dot A=-Ak and dot B=-Bl. Write the quartet as
phi=(u,v n), with independent real u(t,s),v(t,s), momenta P=dot u,Q=dot v.
Both amplitudes evolve, so the norm is NOT frozen. This Cartesian form is
regular at zeros of either amplitude. Let f=1-(u^2+v^2)/6>0,
G=I/f+phi phi^T/(6f^2), U=Lambda/f^2, Lambda=3/2.

The target Christoffel prediction is
Gamma^A_BC=(delta^A_B phi_C+delta^A_C phi_B)/(6f).
Derive/verify it independently and check equivariant closure before runs.
For D=B''-(A'/A)B', spatial Ricci eigenvalues are
Rs=-2D/(A^2 B), Ro=[1-(B'/A)^2-B D/A^2]/B^2.
Let Z=(u'^2+v'^2)/f+(u u'+v v')^2/(6f^2).

    dot k = Rs+K k-Z/A^2-U
    dot l = Ro+K l-v^2/(f B^2)-U
    dot P = K P+[u''+(2B'/B-A'/A)u']/A^2
            -(uP+vQ)P/(3f)+(u u'+v v')u'/(3f A^2)-2Lambda u/(3f)
    dot Q = K Q+[v''+(2B'/B-A'/A)v']/A^2-2v/B^2
            -(uP+vQ)Q/(3f)+(u u'+v v')v'/(3f A^2)-2Lambda v/(3f)

These are evolution equations, not replacements for the constraints.
Check unused Hamiltonian H=Rs+2Ro+K^2-k^2-2l^2-2rho and momentum
M=-2l'+2(B'/B)(k-l)+current, where
current=G(Phi_dot,Phi'), rho=(G(Phi_dot,Phi_dot)+Z/A^2+2v^2/(f B^2))/2+U.
Normalize H by max(1,|Rs|+2|Ro|+K^2+k^2+2l^2+2|rho|), M by
max(1,2|l'|+2|B'/B|(|k|+|l|)+|current|). Also record unnormalized values.

## Initial data, topology, numerics

Hash-verify and replay #307's parent archives and retain its 8/8 verdict.
Use only archived L=5.5, eta=0 and .3 finest reconstructed solutions.
A=B=psi^2; k=2a psi^-6,l=-a psi^-6;
(u,v)=q(sin theta,cos theta), (P,Q)=psi^-6 p q(cos theta,-sin theta),
q=sqrt(3)/2. On [-L,L), u and P are antiperiodic, while all other
evolved variables are periodic. This implements the chosen sign bundle
and antipodal angular map; it must not flip radial particle momentum.

Fourth-order centered spatial derivatives, classical RK4, no filtering,
constraint projection or damping. Fixed N=512,1024,2048 uniform grids;
fixed t in [0,.01], respectively 200,400,800 steps. Record diagnostics
at 201 common times and final fields. Stop and report failure on nonfinite
data, f<=.1, A<=.005 or B<=.005, or failed evolution validation. Do not
extend time or tune resolution if the registered gates fail.

Require at all stored times finest normalized H,M <1e-3 and positive
metric/f. Require medium-to-fine errors in A,B,u,v and neck/bulk Jordan
radius ratio <1e-3 relative to max(1,field magnitude). Coarse/medium versus
medium/fine max differences must have ratio [4,32] if the coarser difference
exceeds 1e-8, otherwise record small-error limitation. Constraints must
decrease across resolutions (unless already <1e-8). These new evolution
tolerances do not change the old initial-data certificate.

## Background separation and validation

Verify the equations against exact round Jordan ESU breathing data:
phi=q cos(2tau) x on unit S3, f=1-cos^2(2tau)/8,
dt=sqrt(f) dtau, A=B=sqrt(f) sech(s), k=l=-f_tau/(2 f^(3/2)).
This predicts initial norm acceleration -4q/f=-3.95897327444315 in
Einstein proper time. Do not mistake Einstein-frame scale change for a
Jordan-frame instability. The exact Jordan scale remains one.

Evolve independent homogeneous Jordan conformal-time controls
a''=-a+a^3, q_c''=-4q_c, dt=a sqrt(f) deta,
f=1-q_c^2/(6a^2), with a0=1,1+1e-4,1-1e-4,
a0'=(a0^2-1)/sqrt(2), q_c0=q, q_c0'=0. DOP853 rtol=1e-11,
atol=1e-13. Check the unused Friedmann constraint to 1e-9. Retain the
homogeneous growing mode; no mean subtraction enters the equations.

Report Jordan bulk and neck areal radii B/sqrt(f), their ratio, bulk
fractional changes and homogeneous controls separately. The short interval
is chosen before measurement (much less than a breathing period), not
selected afterward to hide instability. Do not infer long-time stability
or traversability. The background growth factor exp(sqrt(2) pi) from #296
uses Jordan/conformal ESU time, not the Einstein clock used here.

## Test worldlines and momentum windows

Evolve radial unit-mass Jordan geodesics with unwrapped s and covariant p_s.
H_particle=sqrt(1/f+p_s^2/A^2),
dot s=p_s/(A^2 H_particle),
dot p_s=[f'/(2f^2)+p_s^2 A'/A^3]/H_particle.
Start at s=L-.05 with local speeds .25,.5,.75, plus reflected partners
s=-L+.05 with opposite speeds. Use periodic cubic interpolation of the
metric and its spatial derivatives, RK4 at the field steps; record both
covariant p_s and local orthonormal p_hat=p_s sqrt(f)/A.

Locate s=+/-L by interpolation, never by an imposed kick. Require all
registered curves to cross within the interval, medium/fine crossing times
to differ by <1e-4 and p_hat at crossing by <1e-3. Verify timelikeness,
reflected trajectories and continuous bundle transport. A minimum of the
Jordan areal radius must persist at the seam (checked on adjacent grid
points); otherwise label it a seam crossing, not an evolved neck crossing.

At each crossing measure signed p_s and p_hat changes on symmetric time
windows of half-width .01/16,.01/32,.01/64. Compare all resolutions.
A candidate finite impulse needs |finest jump|>1e-6, both consecutive
absolute-jump ratios in [.8,1.2], and fine/medium agreement <10% of the
candidate jump. Vanishing windows with ratios [1.5,2.5] and a finite smooth
force instead support continuous transfer. Neither a coordinate sign change
nor a finite integral over a chosen window establishes a quantum.

## Matter momentum balance and controls

Use X=partial_s, explicitly not a Killing vector, and the one-sided tube
[L-.2,L]. With W=A B^2, P_tube=4pi integral W T^t_s ds,
T^t_s=-current. Let Lm=(G(Phi_dot,Phi_dot)-Z/A^2-2v^2/(fB^2))/2-U,
Ss=Z/A^2+Lm, So=v^2/(fB^2)+Lm. The required balance is

    Delta P_tube = -4pi integral [W Ss]_edges dt
                   +4pi integral integral W[(A'/A)Ss+2(B'/B)So] ds dt.

The geometric source term is mandatory. Independently integrate the RHS
as an auxiliary RK4 variable using spatial Simpson quadrature at 129 fixed
tube points. Compare to momentum integrated from evolved fields at all
201 output times. Finest defect <1e-5 times max(1,|P_tube|,|integrated RHS|),
and decreasing with refinement unless already <1e-8. Retain omission of
the geometry term as a negative control; demand its peak defect >10 times
the correct peak defect and >1e-8. This is a named matter momentum balance
with gravitational exchange, not an ADM momentum at infinity or a separate
mouth-particle four-momentum.

## Evidence and verdicts

Archive the freeze, input/source hashes, schedules, final fields, all
diagnostics, worldlines and signed window measurements. Replay must verify
inputs and recompute results; malformed or tampered evidence clears verdicts.
Tests cover topology, norm evolution, exact round breathing, geodesic sign,
missing geometry term, constraints and failed-output withdrawal.

Report EVOLVED_TEST_WORLDLINE_CROSSING, COVARIANT_MATTER_MOMENTUM_BALANCE,
FINITE_CROSSING_IMPULSE, and DISCRETE_RECIPROCAL_MOMENTUM_EXCHANGE separately.
The last requires a converged finite impulse AND an independently evolved
receiver with an equal-and-opposite transfer ledger; test particles alone
cannot satisfy it. If no finite impulse appears, report that result without
adding a detector. Do not infer action quantization or quantum statistics.

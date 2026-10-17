# Two moving antipodal mouth pairs: pre-run protocol

This is the next effective-mechanics/visualization step after PRs #338–339.
The two handles of #338 have FOUR mouths. Each pair remains antipodal and
has one independently evolving axis X_h and orbital angular momentum L_h.
The paired constraint is imposed, not derived. L_1 = -L_2 is NOT imposed.

## Equations and scope

Use #339's real spectral scalar on the complete round unit S²,
A = -Delta + 1/4, omega_l = l+1/2. Let x_hs = s X_h, s = +1,-1,
and w_hs = exp[-b² l(l+1)/2] conjugate(Y_lm(x_hs)). Each mouth has its
own deformation coordinate q_hs and canonical momentum p_hs. Set

H = (||p||² + ||omega a||²)/2
  + sum_h |L_h|²/(2 M)
  + sum_hs [p_hs²/(2 mu) + mu Omega² q_hs²/2
            + kappa (q_hs - Re(w_hs* a))²/2].

The torque on each pair is the SUM of its two mouth torques:
T_h = kappa sum_s (q_hs-ubar_hs) Re(w_hs* D a).
Evolve Xdot_h = (L_h/M) cross X_h, Ldot_h = T_h,
pdot = -omega² a + kappa sum_hs (q_hs-ubar_hs) w_hs,
qdot_hs = p_hs/mu, pdot_hs = -mu Omega² q_hs - kappa(q_hs-ubar_hs).
D_i are the three rotation derivatives from #339.
Total J = sum_h L_h - Re(p* D a) is conserved. The constraint supplies
internal pair forces, not another independent momentum reservoir.
M is the inertia of an entire pair; equal endpoint effective masses are M/2.
For R=1, P_hs = (M/2) (L_h/M cross (s X_h)); these velocities are
kinematically opposite within one pair. This is not an inter-handle exchange
law, and ambient Cartesian translation is not a symmetry of S².

Use exact free-field, free-shape and sphere-geodesic flows and reciprocal
interaction kicks (Strang splitting). No focusing switch, damping, prescribed
trajectory, active drive, or corrective kick is permitted. No retuning after
the first coupled run.

This Hamiltonian does NOT evolve #338's handled-mesh wave operator. The
round-sphere caps remain in the field domain; the four finite footprints are
effective couplings. The connecting necks carry no solved scalar or inertia
beyond M. No propagation through the neck, MTY clock shift, GR metric,
viscous-fluid dynamics, or quantum action rule is computed. The antipodal
constraint and smoothing tails preclude claiming strict local causality.

## Initial data and fixed budget

Initial axes: +z and +x, matching #338. Each has L=0 and both mouths q=pq=0.
Define f_x and f_z as #339's smooth target profile centered at +x and +z.
Target a = f_x + 0.65 f_z; target p = -eta (D_z f_x + 0.65 D_x f_z).
Normalize their JOINT field energy to 0.1, then freely backpropagate by pi.
This is a prepared incoming field, not spontaneous mouth emission. It also
interacts with the negative mouths at launch. Therefore launch recoil must
be reported separately from the later positive-mouth focusing event.
Initial coupling energy is additional and must be included in total H.

Parameters: cutoff 32; sigma .20; footprint b=.12; M=.10 per pair;
mu=.05 per mouth; Omega=5; kappa=1. M is the exploratory review movie's
inertia choice, not the original #339 M=.30. eta=.7 directional, 0 symmetric,
-.7 mirrored. The mirror is reflection y -> -y: X transforms as a polar
vector, L as an axial vector diag(-1,+1,-1). This fixes both initial axes.
Additional controls: kappa=0 directional and zero field/zero mouth state.
Symmetric superposed packets need not give zero individual pair recoil:
the other packet breaks axial symmetry. Do not turn that into a null gate.

Horizon pi+1.2; dt<=.002; exact samples at pi-.6, pi, pi+.6.
Refine directional at cutoff 40 and at dt<=.001; output interval .04.
Report each mouth's integrated torque, pair momenta, field momentum, shape,
kinetic and coupling energies, displacement, pair separation and budget error.

## Failure criteria, fixed before running

Numerical admissibility: every case has maximum relative energy error <1e-4,
total-J error / max(H_initial, |J_initial|, 1e-12) <1e-8, unit-axis and
tangency errors <1e-10. The cutoff/time refined FINAL vector of both pair
momenta differs by <max(2 percent of its norm, 1e-7). Mirror full trajectories
must agree within 1e-8; zero field and uncoupled controls remain at rest to
1e-13. A failure here yields UNRESOLVED_NUMERICALLY.

To support TWO INDEPENDENTLY MOVING PAIRS, BOTH directional final pair
momenta must exceed 1e-6 and the angle between the two positive axes must
change by >1e-5 rad at some time. Otherwise report NOT_SUPPORTED for this
family; a spinning camera or rigid global rotation does not qualify.
The minimum distance between any mouths of different handles must stay
above 2*.3+.05 rad; crossing invalidates this finite-mouth picture.
No gate claims discrete action, irreversible absorption, or stronger total
exchange at focus. Report focus-window torque and launch-window torque
without silently attributing the latter to the antipode.

## Drawing contract

Use #338's two Klein handle meshes, with the actual gluing permutation and
shared rim vertices. Rotate each complete drawn neck with its solved axis.
Warp the surrounding sphere mesh smoothly to match the moving rims while
keeping distant coordinate landmarks fixed. This is a mesh-coordinate map,
not a solved fluid velocity or a new scalar evolution operator. Sample the
spectral field at the MOVED sphere locations; moving mesh labels must not
advect a cached field. Keep the neck neutral because no neck field is solved.
Colored longitudinal fibers are mesh tracers, not an extra dynamical field.
Mark initial rims, show actual motion without magnification, and include a
separate clearly labeled magnified angular track and momentum ledger.
Check rim continuity, unchanged non-orientable connectivity, no inverted
sphere triangles, and disjoint mouth collars over archived frames.

Freeze this document in a local commit before the first coupled run. Publish
the hash with results; this is not external preregistration.

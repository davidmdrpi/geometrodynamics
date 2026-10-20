# Coupled mouth–throat candidate: protocol before running

**History and disclosure (2026-10-20):** The original protocol was committed
locally before pilots; its GitHub publication at `5c6e971` occurred after those
pilots. It is a locally frozen protocol, not independent external preregistration.
After the pilots, the horizon was extended from 6π to **8π for every production
case**: three complete acceleration windows had left no later interval in which
to test the subsequent collapse. Equations, parameters, controls and the first
three-window acceleration criterion were unchanged. The initial 6π wording is
retained below to make this amendment visible. The numerical window thresholds
used in the original experiment and replay are Δseparation < −0.001 rad and
Δcombined speed > 0.001; those tolerances were not specified in the original
published prose and are additional implementation choices, not preregistered
thresholds. All baseline window separation changes are positive, so this
distinction does not affect its rejection.

The unpublished implementation was lost in a session reset. The reconstruction
uses this action and the preserved binary initial states as recovery checks;
it does not reselect the parameters. See the results document for the recovery
provenance and full replay comparisons. No claim of an independent replication
is made by reconstructing the same equations against their earlier output.

This is an effective classical **port model**, not a moving excised-cap mesh,
GR solution, or QED derivation. It adds independent meridional mouth positions
and a positive-radius neck degree of freedom to one coupled sphere/interval
scalar field. The sphere remains round and unholed. Finite aperture springs
connect its moving samples to the ends of a resolved interval field. All
forces, including radius forces, are derivatives of one autonomous action.
This is a constitutive choice, not a result derived from PR #345's rim flux.

We test whether a finite packet drives repeated approach, throat shrinkage,
and a regularity event without repeated launches or assigned kicks. The same
field persists throughout the run; empty-sphere antipodal refocusing comes from
omega_l=l+1/2. The interval is a finite transmission path, not an instantaneous
force. No local amplitude gain is applied at either mouth.

## Action and clocks

For field coefficients u, canonical momenta p, mouth coordinates q_A,q_B,
neck radius r and mechanical momenta P_A,P_B,P_r:

    L = 1/2 u_dot^T M(q,r) u_dot + u_dot^T B(q,r) u
        - 1/2 u^T K(r) u - V_port(u,q,r)
        + m/2 (q_A_dot²+q_B_dot²) + mu/2 r_dot² - V_neck(r).
    H = 1/2 (p-Bu)^T M^-1 (p-Bu) + 1/2 u^T K u + V_port
        + (P_A²+P_B²)/(2m) + P_r²/(2mu) + V_neck.

M, B, K are block diagonal between the sphere spectral field and a piecewise
linear interval of length L=.565. The ports couple their fields reciprocally:
V_port=g(r)/2 sum_i [z_i - <u>_footprint(q_i)]², with g(r)=g0(r/r0)².
The selected restoring potential is V_neck=k(r-r0)²/2. Its collapse barrier is
reported, not secretly reduced after seeing the results. Neck weight is
w(r)=(r²+eps²)/(r0²+eps²), with eps=.1r0. This positive floor is an explicit
regularity law; r=0 is never evolved. The event threshold is r_cut=.25r0.

Use external time t=tau+Theta on the sphere, Theta=kappa(z-z_source). On the
interval the smooth clock difference is A=Theta(q_B)-Theta(q_A)-Delta.
A jump Delta at the internal cut supplies the clock holonomy. On the sphere
M_ij=int (1-|grad Theta|²)Y_iY_j, B_ij=int Y_i gradTheta.gradY_j,
and K_ij=(l+1/2)² delta_ij. On interval coordinate xi in [0,1],
M=w L(1-(A/L)²) M_unit, B=w A B_unit, K=w/L K_unit.
Finite ports identify values on the same global tau slice. Their finite
aperture is specified in that chart. This is not exact finite-aperture MTY
matching on moving mouth worldtubes. The mechanical kinetic terms also use
tau; this is a nonrelativistic constitutive model. Clock support is a fixed
background, not a solved differential-aging process.

Initially kappa=-Delta/[cos(q_A)-cos(q_B)] so A=0. The clock profile is then held
fixed, not retuned to rescue evolution. Stop if max(|kappa|,|A/L|)=.9. This
means this chart loses its safety margin; it is NOT proof of a chronology
horizon or collapse. Mouth overlap (separation=.3) also stops the run.

The radius and mouth forces include derivatives of M, B, K and the port
potential. Setting radius or mouths fixed is a static constraint, not an
unreported sequence of corrective kicks. Total H is conserved after launch.
For zero clock offset there is also a rotational Noether check; the fixed
nonzero clock background can exchange angular momentum, so autonomous energy
conservation alone is not claimed as a complete isolated angular ledger.

## Frozen experiment and acceptance gates

Baseline: l_max=10, neck nodes=17, sigma=.30, footprint=.20, L=.565,
r0=sin(.3)², m=.3, mu=.2, k=10, g0=1, packet energy=.08,
q_A=.35, q_B=pi-.35, Delta=-1.2. Initially mechanical velocities and neck field
are zero. One smooth spectral packet starts at q_source=1.1 with directional
velocity toward A. Normalize the complete field plus initial port energy to
.08. Neck cutoff costs .5*k*(.75*r0)² (~.02145), so energy alone does not forbid
it; neither does energy sufficiency predict it. Horizon=6pi, about three
empty-sphere full return times. No packet is relaunched at pi or 2pi.

Controls: no wave; zero clock offset; mistimed offset Delta=-.6; fixed mouths;
fixed radius; temporally dephased packet. Also rerun baseline at tighter ODE
tolerances and higher sphere/neck resolution. Do not tune g0, k, pulse energy,
or r_cut after inspecting outcomes. Report every case, including chart stops.

Success requires all of the following:
1. Mechanical speed increases AND separation falls across three consecutive
   return windows of width 2pi (not three numerical samples).
2. Radius falls by at least 25%, and reaches the positive cutoff if claiming
   a collapse candidate.
3. Energy relative residual <1e-5 and tolerance-refined diagnostics agree;
   wave-off stays at rest; results survive a resolution check.
4. No chart or overlap stop precedes the claimed event. A cutoff crossing alone
   is not evidence of sustained self-acceleration or an MTY mechanism.

No cutoff crossing means no collapse in this candidate over this horizon.
Failure of any gate rejects the full proposed sequence for these parameters.
The same outcome in the no-offset control does not establish advanced-time
causation. A result exclusive to one tuned arrival phase is reported as such.

## Cutoff and release

At r_cut stop and retain the complete field and mechanical state and energy
ledger. Do not erase the interval, duplicate its field on the sphere, or add a
pulse. A positive cutoff is only a regularized collapse **candidate**. An empty
sphere and outgoing burst require a separate topology-transfer law preserving
energy and accounting for momentum. That law will be implemented only after a
candidate meets the acceleration/shrinkage gates; no gamma-ray spectrum or
quantization follows from this scalar experiment.

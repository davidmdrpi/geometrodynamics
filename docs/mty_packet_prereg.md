# Prospective finite-packet MTY response experiment

Date: 2026-10-08 UTC. Base: main `0412d4e053838bea7179a6eeb817dfda4f44fde5`.
Publish this specification AND both new implementation/producer modules before
the validation schedule is run. Preliminary pilots are disclosed below.

## Question and mandatory scope

Does a finite, independently specified source pulse admit a bounded,
self-consistent nonlinear response of two independent localized coordinates
on a time-shifted wave graph, with conservative local scattering, accounted
outgoing flux, a measurable advance, and response changes when the handle is
disconnected? No confirmation wave, equal/opposite impulse, or phase closure
is inserted. No corrective kicks or outcome-dependent timing adjustments.

This is an explicitly conditional REDUCED SCALAR MODEL. It is not a solution
of the Einstein equations, the S3 focusing problem, a solved finite-mouth
metric junction, or a model of gravitational mouth COM recoil. A traversable
channel and the support that keeps it open are assumed. The two coordinates
are dynamical local response modes with canonical momenta, not two oriented
faces re-labelled as independently gravitating particles. The scalar action,
material coefficients and twisted handle boundary condition are assumptions.
They are not derived from vacuum GR, topology alone, or the earlier neck data.
The graph retains only an antipodal channel's travel time, in units where
pi*a/c=1. Each node also has a lead accounting for source injection and outgoing
wave radiation; this is not a complete closed-universe field simulation.

Always report GR traversable support, gravitating-mouth recoil, unique history
and action discreteness as NOT_ESTABLISHED, regardless of the reduced verdict.

## Action and measured quantities

For each wave edge use S_edge = integral dt dx [(phi_t)^2-(phi_x)^2]/2.
At node j identify every incident edge's endpoint field with q_j, with action

    S_node = integral dt [M_j*qdot_j^2/2 - K_j*q_j^2/2 - lambda_j*q_j^4/4].

M=(.1,.15), K=(1,1.4), lambda=(10,15). Each node has three ports: bulk, handle,
lead. With incoming characteristic derivative a_jp and outgoing b_jp,
continuity and variation of this action give

    b_jp = v_j - a_jp
    M_j*vdot_j = 2*sum_p(a_jp) - 3*v_j - K_j*q_j - lambda_j*q_j^3.

Thus dE_j/dt=sum_p(a_jp^2-b_jp^2). The generalized force from each edge is
a_jp-b_jp; the on-site potential transfers generalized impulse to a fixed
anchor. That force and p_j=M_j*v_j are NOT longitudinal radiation momentum
or the GR four-momentum of a mouth. The anchor potential energy is included;
a finite-mass anchor's recoil is outside this reduced action.

Use the unweighted Euclidean canonical state (q_A,q_B,p_A,p_B) for response
comparisons. Its scale and units belong to this declared reduced action.
No conservation of an unspecified "wavefront norm" substitutes for energy.

## Clock preparation and route

Reuse NetworkMouth's actual clock maps. A: tau_A=t_A. B: tau_B=t_B-Delta.
Both rates are one during scattering. Delta is calculated from a prescribed
ideal flat-clock out-and-back excursion of B before the packet experiment:
speed sqrt(3)/2, two inertial segments of duration D/2, then rest. The
acceleration/deceleration events are idealized instantaneous operations.
Proper duration is D/2, hence Delta=D/2. The worldline displacement sums to
zero. Per unit preparation rest mass, the nonzero excursion requires positive
actuator work 2 and recovers 2; net work zero does not mean free preparation.
This is a kinematic clock history, not a moving-throat GR solution. It is not
used to claim that the required stresses or actuator are available in BAM.

Bulk transit L=1, proper handle transit h=.125. Arrival-node delays are

    bulk: A<-B and B<-A both 1
    handle: A<-B h-Delta; B<-A h+Delta.

The B-to-A passage has local duration h>0 even when its exterior delay is
negative. The bulk A->B plus handle B->A cycle has exterior duration
1.125-Delta. D is swept prospectively; it is never solved for packet closure.
The handle multiplies a real scalar by -1, an explicitly twisted line-bundle
matching rule; spatial orientation alone does not force a scalar sign.
A change of basis at B changes BOTH edge signs and the local B coordinate.
Measured energies and canonical states transformed back to the physical
basis must be unchanged.

## Finite packet, self-consistency and controls

The source lead at A carries

    g(t)=amplitude*cos(pi*t)^4*cos(4*pi*t), |t|<.5; zero otherwise.

The B lead has no independent source. Its outgoing response is measured.
All four internal incoming histories and the two nonlinear node-force
histories are solved simultaneously. Each iteration solves retarded node
response to trial incident histories, then imposes the delayed/advanced edge
matching. The complete history, not a time-forward guessed confirmation,
must satisfy the residual gate. Shift arrays with zero padding, NEVER a
periodic FFT wrap. No incoming free field is supplied outside the window.
Finite-window boundary assumptions are tested by tails and extension.

The fixed schedule in mty_packet_probe.schedule has 32 runs:

- D=(0,2.25,2.75,3,3.25), each at amplitudes 1 and 2, each dt=1/128 and 1/256.
  D=0 is a causal control; 2.25 is the zero-cycle-time case; the final three
  test advances over an interval of offsets, not one adjusted phase point.
- D=3 with handle disconnected, both amplitudes and both steps. The unused
  handle ports become outgoing reservoirs and their energy MUST be counted.
- D=3 at each amplitude, dt=1/256: extended interval, alternative initial
  solver seed, midpoint nonlinear force, and B-basis reversal, separately.

Main interval [-16,64]; extended interval [-24,96]. Initial node states at
the left boundary are zero. These intervals were chosen after the disclosed
pilot showed that a short [-4,8] window loses appreciable flux at its edges.
The extent comparison must cover the entire main interval.

Use a trapezoidal linear oscillator solve implemented as a digital filter,
verified against direct two-state solves. For the main method, the quartic
force is its discrete gradient, so local work equals the change of potential
energy. The alternative uses lambda*((q_n+q_n+1)/2)^3. This is an independent
nonlinear discretization, not an independent GR solver. Time refinement and
its energy error are measured separately.

Anderson history solve: alpha=.5, M=10, w0=.01, maxiter=1500, Armijo line
search, target max residual 2e-9*amplitude. Independent acceptance is residual
<=1e-8*amplitude, even if the iteration budget is exhausted; record that
condition separately. Abort nonfinite iterates or |q|>100*amplitude as
NUMERICALLY_UNRESOLVED. Two seeds can test sensitivity, not prove uniqueness.

## Conservation and acceptance gates

Reconstruct all diagnostics from saved incident/force histories and q,v.
Reject altered, missing, reordered or nonfinite evidence; never trust flags.

- Relative network/nonlinear matching residual <=1e-8.
- Integrated scalar-field matching defect / amplitude <=1e-5, with zero
  past integration constant, checked separately from characteristic-derivative
  matching so an advanced boundary cannot hide an additive field mismatch.
- Maximum cumulative local energy error / incident packet energy <=1e-6
  (alternative midpoint <=2e-3).
- Canonical momentum balance defect / source amplitude <=1e-7, including
  the anchor impulse; discrete kinematic defect / amplitude <=1e-8.
- Whole-history input minus outgoing-lead energy minus node energy change,
  divided by input energy, <=1e-4. Include disconnected-handle outgoing flux.
- Incoming+outgoing energy in the first/last two units of each window <=1e-4
  of input; initial+final node energies <=1e-4 of input.
- Relative L2 canonical state difference <=.02 for refinement, window extent
  and alternative discretization; <=1e-5 for alternate seed and basis reversal.

Positive handle energy is also measured on a common HANDLE-clock slice:
integral from tau-h to tau of [b_A(u)^2+b_B(u+Delta)^2] du. Do not add it to
bulk energy at the same exterior coordinate time and call that a global
Cauchy-slice energy: the time-shifted identifications do not provide such a
slice. Whole-history scattering balance and local flux balance are the gates.

Mechanism tests, assessed only after numerical gates:

- A-lead outgoing energy before source onset t=-.5 is <1e-9 of input for
  D=0 and disconnected controls; >1e-5 for each of the three advanced offsets
  and both amplitudes.
- At D=3, disconnecting the handle changes the canonical response by >.01 in
  relative L2 norm for each amplitude.
- At each advanced offset, the amplitude-2 state divided by 2 differs from
  the amplitude-1 state by >1e-5 relative L2: measurable dynamical nonlinearity.

All numerical and mechanism gates -> REDUCED_MTY_PACKET_SCATTERING_SUPPORTED.
Numerically valid but any mechanism gate fails -> REDUCED_MTY_PACKET_SCATTERING_FAILED.
Otherwise -> NUMERICALLY_UNRESOLVED. Every scheduled failure is retained.
No threshold or horizon changes after observing validation outcomes.
No pass licenses unique classical prediction, particle exchange, quantization,
no-signaling, or physical traversability of the repository's trapped mouths.

## Pilot disclosure and provenance

Six pilot runs used D=1 or 3.5, amplitude=.5 and dt=1/64, on [-4,8], [-8,32]
and [-16,64]. They are outside the validation schedule and not counted as
evidence for its verdict. Summaries are in mty_packet_pilots.json. The pilot
solver had maxiter=600 and target 2e-10*amplitude. At [-16,64], both pilot
histories met the proposed physical/residual gates, although the D=1 solver
hit its stricter iteration target limit. This motivated maxiter=1500 and a
2e-9 target, with the declared acceptance gate still 1e-8. Local energy
errors were near roundoff; short-window global deficits were 0.7%-1.4%.
All pilot results, including those deficits, are disclosed. No production
schedule trajectory was observed before publishing the final code/specification.

Production requires a published freeze commit, verifies its bytes with git
show, hashes the sources and records dependency versions. Each completed
case is checkpointed in lossless NPZ/base64. A pinned manifest authenticates
the complete archive inventory. Diagnostic replay recomputes decisions without
rerunning the nonlinear history solver. Source-hash checks bind what ran.

References: Morris, Thorne & Yurtsever, PRL 61, 1446 (1988),
doi:10.1103/PhysRevLett.61.1446; Friedman & Morris, CMP 186, 495 (1997),
arXiv:gr-qc/9411033. The latter's existence results do not establish existence
or uniqueness for this nonlinear graph. Earlier repository context: #216,
#217, #219 and the finite-mouth/traversable-throat audits.

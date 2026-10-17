# Antipodal kick: frozen first-transit protocol

Written before the first coupled run. This is a scalar mechanism probe, not a
GR mouth calculation, a soap-film solver, or a quantization claim.

## Question and model

Does coherent antipodal refocusing concentrate momentum transfer to a movable
deformable attachment, including when the incoming focus is not symmetric?

Use R=c=1 on a round **S² without handles**, with the explicitly chosen operator
A=-Δ+1/4. Its frequencies are ω_l=l+1/2. This is the conformally coupled
ultrastatic 2+1 scalar comparator, not PR #338's minimal scalar. The positive
constant potential is appropriate to this round sphere; it must not be copied
to the negative-Euler-characteristic handled surface. Smooth finite-width
data and a converged harmonic cutoff give finite focus amplitudes. There is
no 1/sin(theta) multiplier, threshold, reset, or focus-triggered force.

Real fields are represented by complex orthonormal spherical harmonics with
the usual a_(l,-m)=(-1)^m conjugate(a_(l,m)) reality condition. The Hamiltonian is

    H = 1/2 Σ(|p_lm|² + ω_l² |a_lm|²)
        + |L|²/(2M) + p_q²/(2μ) + μ Ω0² q²/2
        + κ [q - ubar(X)]²/2.

X is a unit vector locating the attachment, L is its tangent orbital angular
momentum (L·X=0), and q is a deformation coordinate. The smoothed field is
ubar(X)=Σ h_l a_lm Y_lm(X), h_l=exp[-b² l(l+1)/2]. This footprint is a heat
kernel, with smooth global tails, not a compact mouth or a topology change.
The same interaction energy supplies the field source, deformation force,
and translational torque. Thus the proposed coupling is explicit and
reciprocal. There is no viscosity in this first conservative probe.

The target is X=(1,0,0), source antipode (-1,0,0), focus time π. Specify the
*free* target profile f_lm=l(l+1) exp[-σ²l(l+1)/2] conjugate(Y_lm(X)).
Its monopole vanishes. At the target time use a=f and p=-η D_z f, where
D_z=∂φ is the rotation generator. η=0 is symmetric. η=±0.7 supplies opposite
directional flux through the focus, with the same target and spatial profile.
Free backward evolution determines initial data exactly:

    a(0) = -P p(π)/ω,       p(0) = P ω a(π),       P_l=(-1)^l.

This identity follows from cos(πω_l)=0, sin(πω_l)=(-1)^l. It is a
phase-space refocus, not a claim that every displacement profile repeats
unchanged after a half circuit. Normalize each launch to wave energy 0.1.
Initially X is at the target and L=q=p_q=0. Total energy additionally includes
the initial coupling energy, which must be reported rather than hidden.

## Cases and fixed parameters

- Symmetric η=0; directional η=+0.7; mirror η=-0.7.
- Dephased directional control: rotate each l oscillator's (ωa,p) pair by a
  deterministic random phase (seed 338), shared by all m at that l. This
  preserves every mode's energy and total field angular momentum, but not
  the localized source envelope. It tests coherence, not identical source
  support. Record this limitation.
- Uncoupled directional control κ=0; half-amplitude directional run (energy
  0.025), to inspect continuous scaling rather than assume discrete action.
- σ=0.20, b=0.12, M=0.30, μ=0.05, Ω0=5, κ=1.
- First-transit horizon T=π+1.2; focus window W=[π-0.6,π+0.6]. No later
  recurrence is interpreted in this study. Observe throughout [0,T].
- Harmonic cutoffs 24,32,40; production dt≤0.002; time check dt≤0.001.
  Exact free oscillator/geodesic flows and reciprocal interaction kicks
  are composed with Strang splitting. Output at ≤0.02 intervals.

## Diagnostics and pre-stated failures

Use global rotational Noether momentum J=L-Σ Re[p* D a], not a sum of
untransported tangent momenta. Check all three components, energy, sphere
constraints, and harmonic reality. Report torque impulse ΔL in W and at T,
maximum displacement, signed directional recoil, and deformation energy.
Integrate |torque| to measure how much *activity* lies in W, separately from
the signed net exchange (which can cancel). Plot both.

The numerical calculation is unresolved if normalized energy drift exceeds
1e-4, total angular-momentum drift exceeds 1e-8 (scale max(E0,|J0|)), or the
32→40 and dt/2 changes in directional final ΔL exceed 2% (absolute floor
1e-7). Symmetric transverse recoil must be below 1e-8; mirrored recoil must
reverse within 1e-7; the uncoupled attachment must remain at rest.

For this fixed family, call the **localized recoil hypothesis supported**
only if the directional case passes the numerical checks, at least 75% of
integrated torque magnitude lies in W, the coherent |ΔL_W| is at least twice
the dephased value, and directional |ΔL_W| exceeds 1e-5. Otherwise report
failure or unresolved evidence, without changing the window or parameters.
Also report post-window evolution: a transient kick that later cancels is
not a lasting transfer. The window and comparison factors are study choices,
not universal constants. This test cannot establish action quantization.

## Visualization

Show actual solved displacement on three spheres (symmetric, directional,
dephased), the free target and moving attachment, and simultaneous angular
momentum/torque traces. Use a fixed common displacement scale; label any
display exaggeration. The attachment's position comes from X(t), never from
a manually prescribed antipodal jump. Show its q response separately from
its motion. Archive numerical time series and renderer inputs.

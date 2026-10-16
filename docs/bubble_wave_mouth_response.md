# Bubble-wave interpretation: refocusing and a carried mouth

The intended picture is a closed background carrying a diffuse deformation
wave between antipodes, with mouths needed for a localized interaction or
measurement event, not for the background wave to refocus. A small attached
feature is carried, stretched and lifted by that wave, then returns. The
outgoing wave should recover its undisturbed profile in the appropriate
transparent limit. This is a useful **target to test**, distinct from the
existing handled-surface scalar calculation.

## 1. Refocusing does not require a mouth

On the round background, rays from a point converge at its antipode and
continue beyond it. In a polar-distance plot this looks like a reflection:
the distance from the original source begins decreasing again. No reflecting
wall or throat needs to be inserted at that coordinate turning point.
A smooth finite-band scalar solution need not acquire a curvature singularity
there: a ray caustic and a spacetime curvature singularity are different.

Undamped propagation, however, does not by itself guarantee an identical
packet on every circuit. For the continuum equation currently approximated by
the **bare S²** viewer,

    u_tt = c² Δ_S² u,          ω_l = (c/R) sqrt(l(l+1)).

These frequencies are not equally spaced. The field can keep oscillating
without dissipation while a generic packet changes shape through dephasing.
Mesh dispersion is an additional effect, not the only source of imperfect
recurrence. Antipodal parity Y_lm(−x) = (−1)^l Y_lm(x) is a spatial identity;
turning it into a fixed time-translation identity requires compatible mode
frequencies. Attaching actual handles also changes the spectrum and breaks
the round background's focusing property.

An explicit comparator that does have exact recurrence is the conformal
scalar on a static round **S³** background:

    u_tt = c²(Δ_S³ − 1/R²)u,
    −Δ_S³ Y_l = l(l+2)Y_l/R²,
    ω_l = c(l+1)/R.

Writing T = πR/c and P u(x) = u(−x), each oscillator gives

    (u, u_t)(t+T) = −P (u, u_t)(t),
    (u, u_t)(t+2T) = (u, u_t)(t).

This follows from cos((l+1)π) = (−1)^(l+1) and sin((l+1)π) = 0.
It holds mode by mode for a finite expansion and extends to appropriate
finite-energy solutions. It is a comparator with a specified curvature term,
not a property already proved for the current S² handle solver, minimally
coupled scalar theory, or tensor gravitational waves. The S² drawing of this
comparator would be a slice of an S³ field, not an S² wave solve renamed S³.

## 2. A feature carried by the wave is a different control model

A small attached bubble suggests applying a local material displacement map
F_t to the **whole feature**: its center and shape respond to the surrounding
surface displacement and gradient. A uniform displacement translates it;
a spatial gradient stretches or tilts it. If F_t returns to the identity,
the feature returns as well. Diffuse means its variation scale is large
compared with the feature radius, not merely that its amplitude is small.

A passive marker defined this way has no reaction on the carrier wave.
Its transparency is therefore exact **by construction**. That is a valuable
visual control, but it cannot demonstrate that a gravitating mouth is
transparent. The current viewer instead modifies the propagation domain by
adding an intrinsic neck: reflection, transmission, and altered interference
are possible even though the displayed mesh itself does not feed back.
PR #338's ambient-vector breathing is also not a solved material response of
an attached bubble or mouth.

For a physical mouth, return to equilibrium requires its inertia, restoring
law, and coupling to the scalar or metric perturbation. A conservative
coupled mouth can retain oscillatory energy after a pulse; damping can relax
it while transferring energy elsewhere. Even where a long-wavelength limit
suppresses scattering, exact restoration of the outgoing wave and exact
mouth recovery must be measured, not supplied as drawing rules.
The identification of a mouth interaction with measurement remains the BAM
postulate to be tested, not a result of this visualization.

## 3. Memory is a separate observable

A bound object returning to equilibrium is not equivalent to absence of
gravitational-wave memory. Restoring forces can erase a detector's residual
shape even when freely falling test bodies retain a displacement. Nonlinear
memory is sourced by radiated gravitational-wave energy and is not restricted
to unequal masses, recoil, or a visibly asymmetric pulse. Angular structure
of the integrated flux matters; antipodal symmetry alone does not eliminate
its even multipoles.

For the standard asymptotically flat discussion, see Marc Favata,
[The gravitational-wave memory effect](https://arxiv.org/abs/1003.3486),
especially sections 1 and 4–6. That asymptotic construction is not directly
a memory definition for a closed S³ universe. In BAM the next probe should
measure a relational residual (separation, strain, and relative velocity of
specified tracked bodies) in a preselected interval after a passage and
before the next return. Persistent oscillation, permanent displacement,
velocity change, and numerical drift must be reported separately.

## Proposed next comparison (not executed here)

1. **Bare round background:** measure antipodal overlap and full-return error
   for a fixed packet over several circuits. Keep the S² minimal and S³
   conformal operators explicitly separate. Compare both field and velocity,
   with the predicted antipodal sign, rather than a moving intensity maximum.
2. **Passive attached feature:** advect its center and shape with the sampled
   carrier deformation and gradient. Confirm exact carrier equality as an
   implementation check; do not count that equality as physical evidence.
3. **Coupled mouth:** prescribe a conservative response law and its energy
   ledger, then scan packet wavelength / mouth radius. Measure residual
   displacement, velocity, internal excitation, reflected/transmitted energy,
   and outgoing profile error relative to the bare run on a common exterior
   observation region. Report any time shift separately from shape error.

Before running case 3, set the transparency/recovery tolerance, packet family,
observation window, and refinement criteria. A residual that remains above
that tolerance under refinement, or retained mouth excitation inconsistent
with the claimed return, counts as failure of the transparent-mouth law for
that regime. No reset of the field or mouth may be inserted at the antipode.
This comparison can support the bubble intuition in a measured limit without
making it impossible for the model to fail.

## Second-review qualifications

The [follow-up audit](bubble_wave_review_audit.md) makes the recurrence
statement explicit for Maxwell and unsourced TT modes, derives the
Gauss–Bonnet obstruction to importing a conformal potential onto the
negative-Euler-characteristic 2D surface, and separates reviewer-reported
low neck uptake from demonstrated transparent passage. The next baseline
should use an independently converged radial/spectral bare solution;
near-caustic comparisons between two unrelated meshes are not yet adequate.

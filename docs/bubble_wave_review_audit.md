# Audit of the second PR #338 review

The [independent follow-up review](https://github.com/davidmdrpi/geometrodynamics/pull/338#issuecomment-6185049784)
reproduced all 36 tests at `d622909`, checked the rendered fixes, and accepted
the scope of the bubble note. Its additional physical claims need the
qualifications below. This is a review audit, not a new transmission experiment.

## Which fields repeat on the round background?

| Field and assumptions | Frequency | Consequence |
|---|---|---|
| Minimal scalar on round S² | c sqrt(l(l+1))/R | Generic packets do not exactly repeat each circuit |
| Conformal scalar on static round S³ | c(l+1)/R | Exact full return at 2πR/c; half-return is −P |
| Source-free Maxwell on the round Einstein cylinder | c n/R, n=2,3,... | Exact full return at 2πR/c, with vector/tensor transport treated correctly |
| Unsourced TT modes on static round S³, vanishing perturbed anisotropic stress | c sqrt(n²−1)/R, n≥3 | Generic packets have modal phase slip |

For the Maxwell spectrum see Lechtenfeld,
[Electromagnetic knots from de Sitter space](https://arxiv.org/abs/2102.13534),
section 4 (the Maxwell problem is solved on the conformally related round
cylinder). Integral frequencies establish full periodicity; they do not
license applying the scalar's −P rule to individual Cartesian components
without the appropriate pullback of the electromagnetic field.

For the standard TT equation, the tensor Laplacian eigenvalue is
`−(n²−3)/R²`; the curvature term adds `2/R²`, giving `ω_n²=(n²−1)c²/R²`.
This agrees with the lowest frequency sqrt(8)/R in
`geometrodynamics/waves/backreaction.py`. The perfect-fluid background
calculation is also described by Barrow et al.,
[On the Stability of the Einstein Static Universe](https://arxiv.org/abs/gr-qc/0302094),
equations (22)–(23). Relative to the integer comparator, the phase deficit
per half circuit is

    δ_n = π(n − sqrt(n²−1)) = π/(n + sqrt(n²−1)) ~ π/(2n).

This obstructs exact recurrence of a generic broadband TT packet in that
sector. It does not say that every gravitational-wave solution is aperiodic:
an individual normal mode is periodic. Nor is it a theorem for every BAM
support law. The repo's [coupled multiplet response](coupled_multiplet_response.md)
and [supported FRW tensor calculation](frw_supported_tt.md) already change
the tensor equation when the support responds. Their frequencies/Floquet
maps must be used in those sectors. Maxwell recurrence cannot be transferred
to photons modelled as tensor gravitational waves without deriving that map.

## Do not transplant the conformal term onto the handled cartoon

In an ultrastatic 2+1 scalar model the conformal coupling gives

    u_tt = Δu − (K/4)u,       A = −Δ + K/4,

where K is the Gaussian curvature. On a closed surface, the constant trial
function has quadratic potential energy

    E[1,0] = (1/8)∫K dA = πχ/4.

For either two-handle topology χ=−2, this is negative. Equivalently the
Rayleigh quotient of A on the constant test function is
`πχ/(2 Area) < 0`, so A has a negative eigenvalue and the wave equation has
an exponentially growing mode. The constant test function need not itself
be an eigenfunction when K varies. This argument also works for a
non-orientable surface, using its area density or orientable double cover.
It presumes a closed metric with the usual curvature matching; any additional
singular seam curvature must be included in Gauss–Bonnet.

Thus adding K/4 would introduce a new instability, not repair recurrence.
The current minimal scalar solver is unchanged. The S³ conformal reference
must remain a separate model. The negative-energy test alone does not settle
χ=0 or positive-χ cases.

## What the reported transmission measurements establish

The reviewer reports substantial angular selectivity: concentric scalar
arrivals couple more strongly than oblique arrivals, with oblique uptake
falling as the mouth shrinks. That is a promising hypothesis for the next
probe. The numerical values remain **reviewer-reported, not independently
reproduced or archived by this PR**. A reproducible baseline still needs the
exact source, frequency content, mesh, observation window, flux estimator,
interpolation, and common energy normalization. In the comment, concentric
fractions are normalized to the incident half; oblique fractions are
normalized to the total. They cannot be divided directly as equal-budget
transmission fractions. Mouth angle and physical rim radius must also be
identified separately.

The angular-selectivity mechanism is plausible already at the scalar level.
For the intrinsic neck metric `ds²+r(s)² dψ²`, an angular Fourier component
has the centrifugal term `m²/r²` in its radial operator. A narrow waist can
therefore discriminate between rim-averaged and higher-m components. This
observation does not reproduce the reported 52–85% or 0.4–8% transmission
numbers, and it does not identify that 2D m with a physical GW multipole l.

**Small neck uptake is not yet transparent exterior passage.** Energy can
be reflected, delayed, temporarily stored, or redistributed on the same side.
A fixed mouth also cannot demonstrate lifting and recovery of a dynamical
mouth. Keep low uptake, small exterior field error, and mouth recovery as
three separate observables.

The reviewer reports a bare-run mismatch around 0.26 at 24000 vertices near
the antipodal caustic, comparable to a handled-versus-bare difference of
0.24–0.28. These are not yet a certified error bound, but they show why this
comparison cannot support an unperturbed-reformation claim. Use a converged
axisymmetric/spectral bare reference, then sample both solutions on a common
exterior region away from the focus. Refine spatial and temporal resolution,
and interpolation, independently. Near-focus diagnostics need their own
convergence study. Do not assign a 1e-6 accuracy label until it is achieved.

## Gravitational-wave qualifications

Pure vacuum Einstein gravity in 2+1 dimensions has no local propagating
spin-2 modes. A scalar on this surface is a stand-in, not a 2+1 GW solution.
In 3+1 spherical vacuum symmetry there is no radiative monopole. But that
means no spherically symmetric **field**: a quadrupolar radiation field can
still have geometrically spherical wavefronts with angularly varying
amplitude and polarization. Uniform scalar-ring results cannot be copied
into that tensor sector.

The review's `T_l ~ (ωb)^(2l+2)` is not accepted here as a universal throat
law. Such exponents require a specified dimension, radial potential,
asymptotics, junction/support conditions, and a regime excluding resonances.
For comparison, explicit low-energy wormhole scattering derives amplitudes
from the effective potential; see Azad, Loran and Mostafazadeh,
[Transmission of low-energy scalar waves through a traversable wormhole](https://arxiv.org/abs/2010.15023).
That scalar result is not a derivation for the BAM tensor channel. In
particular, blocking throat transmission can mean strong reflection; it does
not by itself imply that an incident gravitational wave passes the object
unperturbed. The coupled-mouth protocol must derive and test its own
reflection, transmission, phase, and residual-motion observables.

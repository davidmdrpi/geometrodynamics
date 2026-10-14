# Reciprocal-map packet evolution and leading Einstein mass response

**The finite packet propagates, but the fixed reciprocal null junction fails
induced-metric matching at O(A²).** This completes the linear 2D characteristic
packet calculation and its leading spherical gravitational mass/area response.
It does not construct a nonlinear backreacted spacetime or a surface matter law.
The distinction is physical: the two otherwise consistent exterior responses
cannot be joined with the prescribed reciprocal null identification.

The [production protocol](reciprocal_packet_protocol.md) was written before
implementation and production, after the analytic compatibility issue had been
identified. #334's submitted protocol and preflight are unchanged. The numerical
thresholds and order of gravitational approximation here are stated in that
new protocol, not represented as a pass of every gate in the original submission.

## Specified characteristic problem

Set b=c=kappa5=1, with Einstein equation G_ab=kappa5 T_ab. The scalar is constant
on S3 for the gravitational calculation; it is not the unit-normalized Y_0.
Thus angular integration supplies Omega3=2 pi². A=0.01 is the amplitude of the
rescaled field u=r^(3/2) Phi. The background is vacuum 5D Tangherlini, with no
compact cosmological support field and no Gauss-Bonnet term.

On the incoming rectangle U in [-8,0], V in [.5,4], prescribe the smooth compact
packet from the protocol on U=-8 and zero data on V=.5. Solve

\[
 u_{UV}+Q_l u=0,\qquad
 Q_l={F\over4}\left[{l(l+2)\over r^2}+{3f\over4r^2}
                                      +{3b^2\over r^4}\right].
\]

The nonzero horizon coefficient from #334 is retained. At U=0 extract a(V).
On the second exterior, prescribe u_out(-c/V,0)=(-1)^l a(V) and zero data on
U_out=-2. Its second coordinate covers V_out in [0,80]. The characteristic
corner data agree. The outgoing rectangle therefore has a complete Goursat
problem; it is not a one-dimensional trace identity or an assumed free-wave
translation. The incoming horizon trace can have a scattering tail even though
the original boundary pulse is compact.

Both l=0 and l=1 are evolved. Only l=0 receives a spherical gravitational
response. A single l=1 stress tensor is not spherically symmetric, and no exact
Misner-Sharp mass is assigned to that control. We also evolve a second independent
real solution with initial profile (V-2)a to test a nonzero symplectic pairing.

The midpoint diamond discretization uses the true cell sizes on the nonuniform
outgoing U=-1/V grid. A manufactured nonzero-potential solution tests convergence
independently of Tangherlini, and a zero-potential control tests the two-leg data.
The observer field is interpolated on r=2b, UV=-(b²/3)exp(4), with its clock
origin transformed under passive Kruskal boosts. It is a finite observation
window, not an integral over all future radiation or infinity.

## The mass response is an Einstein constraint calculation

For ds²=-F dU dV+r² dOmega3², F=exp(2sigma), a canonical spherical scalar obeys

\[
 r_{UU}-2\sigma_U r_U=-{\kappa_5 r\over3}\Phi_U^2,\qquad
 r_{VV}-2\sigma_V r_V=-{\kappa_5 r\over3}\Phi_V^2,
\]
\[
 r_{UV}=-{2\over r}(r_Ur_V+F/4),\qquad
 \Phi_{UV}=-{3\over2r}(r_U\Phi_V+r_V\Phi_U).
\]

Differentiating the geometric mass definition and using these equations gives

\[
 m=r^2(1+4r_Ur_V/F),\quad M={3\Omega_3\over2\kappa_5}m,
\]
\[
 m_U=-{4\kappa_5 r^3\over3F}r_V\Phi_U^2,\qquad
 m_V=-{4\kappa_5 r^3\over3F}r_U\Phi_V^2.
\]

On Tangherlini, r_U=-FV/(2b), r_V=-FU/(2b). The computed O(A) field supplies
the complete O(A²) mass source. Corrections to the metric in this source enter
at higher order. The mass integration constant is zero perturbation at the
vacuum corner, not arbitrary on each generator. Two rectangular integration
paths give an independent conservation/integrability check. No Vaidya ansatz,
geometric-optics assumption, or zero outgoing potential is used.

The resulting delta m is a constrained leading response, not a separately
evolved nonlinear metric. In particular, no purported exact geometry is built
by substituting a varying m into the static metric. Amplitude scaling is an
identity at the computed order and is not evidence for nonlinear accuracy.

At the observer, compute the outward Killing-energy flux directly from the
scalar stress tensor:

\[
 {dE_{out}\over dt}=-\Omega_3 r_{obs}^3\Phi_t\Phi_{r_*}.
\]

Its integral agrees with minus the independently reconstructed local mass change.
These are energies normalized to the stationary Killing time. A static local
observer's energy is larger by 1/sqrt(f(r_obs)); their proper time is
sqrt(f(r_obs))t. Neither raw integral T_VV dV nor a clock-origin shift is treated
as a physical energy change.

## Why the geometric continuation fails

The horizon constraints must also agree on the induced radius r(V). A continuous
scalar and its tangential derivative give the same T_kk on both sides when the
common generator parameter is V. At zeroth order the incoming affine generator
has inaffinity 0; the reciprocal outgoing generator has inaffinity
g''/g'=-2/V. This is the pressure-bearing shell from #334.

Write r=b+delta r, Phi=Phi_1+O(A³), and
S(V)=kappa5 b (Phi_1')²/3. With vacuum data and zero expansion before the pulse,
the two constraints at O(A²) independently require

\[
 \delta r_{in}''=-S,\qquad
 \delta r_{out}''+{2\over V}\delta r_{out}'=-S.
\]

Their difference is computable without fitting a response law:

\[
 \delta r_{in}(V)-\delta r_{out}(V)
 =-{1\over V}\int_{V_0}^{V}(V-s)^2S(s)\,ds<0
\]

once a nonzero pulse has passed. This is a positive-kernel obstruction, not a
small residual to remove by increasing numerical resolution. Requiring one
shared radius and subtracting the constraints gives (2/V)delta r'=0; either
constraint then requires S=0. The packet violates that requirement.

Perturbations of the lapse or the map of O(A²) do not cure this: they multiply
the zero background expansion and cannot cancel the nonzero zeroth-order
inaffinity jump times delta r'. This rules out a regular perturbative completion
of **this fixed null reciprocal background with the continuous canonical
scalar and the stated vacuum-past data**. It is not a theorem against every
regular focus law, a surface that changes causal character, additional bulk
stress with a jump in T_kk, or a nonperturbative change of background matching.

On the two independently evolved horizon faces the leading masses obey
delta m_out(V)=-delta m_in(V). Their difference is reported as a quasi-local
mass-parameter jump. It is not an independently derived shell Hamiltonian,
and matching energy flux alone cannot repair the induced-radius mismatch.

## Numerical evidence

N denotes incoming U intervals. Incoming V and both outgoing coordinates have
2N intervals. Thus the finest incoming grid is 513 by 1025 and the outgoing
grid is 1025 by 1025, for each of two independent scalar variations.

| N | Incoming mass-path relative error | Outgoing mass-path relative error | Radius mismatch at V=4 |
|---:|---:|---:|---:|
| 64 | 4.160e-4 | 8.604e-4 | -1.24595055e-4 |
| 128 | 1.047e-4 | 2.181e-4 | -1.25228451e-4 |
| 256 | 2.619e-5 | 5.473e-5 | -1.25387478e-4 |
| 512 | 6.551e-6 | 1.370e-5 | -1.25427277e-4 |

At N=512, A=.01, V=4:

- delta r_in=-2.48870442e-4, delta r_out=-1.23443165e-4.
- delta m_in=+5.10586557e-4, delta m_out=-5.10586557e-4.
- Observer outward Killing energy over t=[.75754668,2.83698822]: .01510646277.
- Observer mass change over that window: -.01510644813; relative discrepancy
  9.69e-7. This is finite-window mass balance, not all-time transmission unity.
- l=0 symplectic relative residuals: 4.80e-7 incoming, 1.03e-6 outgoing.

The archive includes both angular controls, measured field convergence orders,
boosts .5,1,2,4 at fixed physical grid points, and A/2 and 2A checks. The shared
physical grid under boosts tests covariance, not discretization independence.
It also records numerical gate results separately from junction admission.

## Reproduction and artifacts

```
python -m experiments.throat_dynamics.reciprocal_packet_probe
pytest -q tests/test_reciprocal_characteristic_packet.py tests/test_null_horizon_continuation.py
```

[Result manifest](../experiments/throat_dynamics/runs/20261014_reciprocal_packet/result.json)
contains source hashes and package versions. The same directory holds finest-grid
observer and seam arrays in `waveforms_N512.csv` and full 2D field checkpoints
at N=64 for both l values in `fields_N64_l0.npz` and `fields_N64_l1.npz`.
The smaller checkpoints are explicitly labelled; the producer reproduces all
four resolutions rather than storing every large intermediate array.

The next model must supply a response law that satisfies both null constraints
and induced-metric continuity. A different pressure-free null soldering on a
finite generator interval is one possible test; a finite-thickness focus is
another. Neither is silently substituted for the reciprocal map tested here.

Background reference: [Birukou et al., spherical scalar collapse in arbitrary
dimension](https://arxiv.org/abs/gr-qc/0201026). The equations, matching constraint,
and mass normalization used here are explicitly derived and tested above;
the paper is not claimed to establish this antipodal continuation.

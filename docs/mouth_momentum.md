# Nontrivial vacuum momentum completion on a compact twisted handle

Public freeze: [#303](https://github.com/davidmdrpi/geometrodynamics/pull/303),
`2b1bfca7b2650db7f8dcf74ca1c54ecd8bdffa17`, based on merged #300 at
`4cd86541d3b836b35561b0c4a3a54629d28851cd`.
See [the unchanged specification](mouth_momentum_prereg.md).

**Result: all eight frozen gates pass.** The repository now has an explicit
nonspherical, nonzero vacuum momentum-constraint completion and nonlinear
Hamiltonian solve on a compact spatially antipodal handle. The solution has
a reflection-symmetric minimal section and smooth tensor matching across
the identification. No matter, surface stress, prescribed mouth force or
quantum weighting is introduced.

**Scope: this is a new initial-data benchmark.** Its manifold is the
nonorientable S2 bundle over S1 with monodromy n -> -n. It does not insert
a localized mouth into #300's four-scalar universe, retain the earlier
round-S3 exterior, or establish traversability. The section is anti-trapped
on the chosen branch. A worldline crossing has not been evolved, and no
radiative momentum transfer, discrete event, action quantum or probability
law has been calculated.

## 1. What was missing and what has changed

The [earlier physical neck](which_throat_is_physical.md) has K_ij=0. For
zero matter momentum this satisfies the momentum constraint identically.
Its spatial geometry and Hamiltonian-response work therefore do not test
nontrivial extrinsic-curvature completion at a mouth.
[The #300 construction](nonlinear_supported_tt.md) does solve nontrivial
momentum constraints, but in a homogeneous four-scalar sector without a
mouth. This experiment addresses their separation with a distinct, entirely
vacuum compact-handle construction.

The new data are not obtained by placing waves on a fixed background and
declaring their stress small. A deliberately nontransverse gravitational
seed is corrected by a solved longitudinal field, and the metric conformal
factor responds to the resulting tensor through the nonlinear Hamiltonian
equation. The full physical constraints are then checked independently.

Vacuum has j_i=0: gravitational wave data reside in gamma and K, not in an
added matter momentum density. Nonzero K and its transverse completion do
not by themselves measure a propagating wave's energy flux or a particle's
momentum. The freely chosen axial seed is called wave-like initial data
only in this limited sense; propagation requires an Einstein evolution.

## 2. The manifold and tensor completion

Take

    (s+2pi,n) ~ (s,-n), gbar=ds^2+dOmega^2, Rbar=2.

The orientable double cover has period 4pi. A cut at s=0 exposes two sphere
faces, identified antipodally. There is no physical boundary and no shell
action; the fields must descend as smooth tensors. The twist acts in space,
with the timelike normal unchanged. This does not implement time reversal
through a throat or demonstrate retrocausality.

Use gamma=psi^4 gbar and K=psi^-2 Abar, with tr K=0 and psi>0. The background
tensor has A0_ss=2C and A0_AB=-C Omega_AB. C=1/2 and Lambda=1/4 are chosen
benchmark inputs in units of the conformal sphere radius, not derived
constants. They give |A0|^2=3/2 and the exact baseline psi=1.

For the axial l=2 sphere covector V_phi=3 cos(theta) sin(theta)^2, let
S_AB=D_A V_B+D_B V_A. Direct symbolic covariant differentiation verifies

    div V=0, tr S=0, div S=-4V,
    |V|^2=9u^2(1-u^2), |S|^2=18(1-u^2)^2, u=cos(theta).

Both V and S are odd under antipodal pullback. Start with the nontransverse
seed T_sA=epsilon sin(ks) V_A. The correction W_A=epsilon w(s)V_A obeys

    w''-4w=-k cos(ks).

For k=1/2 the solution and corrected tensor are

    w=k cos(ks)/(k^2+4),
    Abar_sA=epsilon [sin(ks)+w'] V_A,
    Abar_AB=-C Omega_AB+epsilon w S_AB, Abar_ss=2C.

The half-integer radial dependence and odd angular pullback together give
a single-valued tensor on the quotient. Treating its components as scalar
periodic functions would give the wrong gluing rule.

The momentum equation was also solved independently with centered finite
differences on the double cover. Its maximum errors against the exact
solution are 2.22093e-5, 5.55689e-6 and 1.38951e-6 on 32, 64 and 128 points:
ratios 3.99671 and 3.99918. The Hamiltonian solve uses the exact completed
tensor, not the finite-difference approximation. This distinction is
intentional: the numerical ODE demonstrates completion and convergence,
while its exact solution avoids importing its truncation error into the
nonlinear metric solve.

## 3. Nonlinear Hamiltonian response and independent constraints

The nonlinear scalar equation is

    Delta_bar psi-psi/4+|Abar|^2 psi^-7/8+Lambda psi^5/4=0,

with

    |Abar|^2=6C^2+2 epsilon^2(sin(ks)+w')^2 |V|^2
                       +epsilon^2 w^2 |S|^2.

Its scalar source is even in u and 2pi-periodic in s. The solve uses Fourier
collocation in s and Gauss-Legendre collocation in u, with positive Newton
line search and a spectral linear preconditioner. Every run starts at
psi=1; no failed amplitude or alternate branch was discarded.

All 15 registered primary combinations were retained: epsilon=0,.02,.05,
.1,.2 on grids (24,12), (40,20), (64,32). The finest solution is also checked
on a separate (96,48) evaluation grid.

| Check | Observed result |
|---|---:|
| Largest collocation equation residual | 3.01e-12 |
| Largest finest-solution off-grid residual | 3.13e-12 |
| Largest medium-to-fine scalar difference | 2.89e-15 |
| Correct tensor seam mismatch, including first derivative | 5.55e-12 |
| Deliberately wrong antipodal seam mismatch | 0.3470 |
| Independent normalized Hamiltonian residual, finest coordinate step | 1.84e-7 |
| Independent normalized momentum residual, finest coordinate step | 5.14e-9 |

The physical check does not substitute the solved scalar equation back into
a conformal curvature identity. At 30 off-grid coordinate points it forms
gamma and K, differences them, builds Christoffels and their derivatives,
contracts the Ricci tensor, and evaluates the physical covariant momentum
divergence. At coordinate steps .001,.0005,.00025 the largest normalized
Hamiltonian errors are 2.84e-6,7.13e-7,1.84e-7; the momentum errors are
8.23e-8,2.06e-8,5.14e-9. The residuals decrease at the expected second order
until differencing roundoff affects the curvature calculation.

The damaged conformal-weight control K=psi^+2 Abar gives Hamiltonian and
momentum errors .02129 and .009531. Omitting or damaging the momentum
correction also fails the physical divergence test. These controls show
that the independent checker detects errors hidden by a scalar solve alone.

## 4. What the minimal section means

Reflection symmetry gives psi_s=0 at s=0. Therefore the physical section
mean curvature H=4 psi_s/psi^3 vanishes there. Its area is computed using
the induced metric, Area(s)=2pi integral psi(s,u)^4 du.

| Wave amplitude | Area(0) | Area(pi) | Future null expansions at s=0 |
|---|---:|---:|---:|
| .02 | 12.567528 | 12.569754 | .999860 to .999866 |
| .05 | 12.573604 | 12.587507 | .999127 to .999162 |
| .10 | 12.595287 | 12.650751 | .996517 to .996656 |
| .20 | 12.681741 | 12.901354 | .986225 to .986768 |

Areas at s=.05 and .1 also exceed Area(0) in every positive-amplitude case.
These are section-area diagnostics and exact minimality by reflection;
they are not a stability proof against all nonspherical surface variations.

With the frozen convention K=-L_n gamma/2,

    theta_+ = H-tr_surface K,
    theta_- = -H-tr_surface K,
    tr_surface K=-2C psi^-6.

Consequently both future expansions are positive at the section: it is
anti-trapped. Full time reversal changes their signs, giving a trapped
section. Neither result establishes a traversable wormhole. This limitation
is a measured part of the result, not a condition suppressed to call the
constraint solve a measurement mechanism.

## 5. Continuity, topology and what remains open

The scalar linearization at epsilon=0 is Delta_bar-5/4, which has no kernel
on this compact domain. Smooth dependence of the elliptic equation gives
a local continuous family of positive solutions near the baseline by the
implicit-function theorem. This is not a proof of evolution stability or a
rigorous finite-amplitude bound for every value up to .2.

The amplitude-reversal control leaves psi unchanged and reverses the axial
part of K. Full time reversal leaves psi unchanged and reverses all of K.
The untwisted S2 x S1 control, with k=1 and ordinary periodic gluing, also
solves the constraints. Thus this experiment establishes admissibility of
the chosen twist, not its dynamical selection.

The half-integer radial label is imposed by tensor descent on the chosen
topology. The amplitude remains continuously variable. Neither that label
nor the existence of a neck is action quantization. No Born weights,
quantized carrier, detector threshold, crossing rule or Phi selection enters
the computation.

The remaining physical task is substantial: embed or match a localized
mouth to the intended bulk configuration, obtain an admissible spacetime
development with the required throat behavior, and calculate a reciprocal
wave/mouth response through an actual crossing. The open four-scalar
interface identified in #302 is not resolved by switching to this vacuum
benchmark. The advance is a tested compact-handle constraint solver and
nontrivial gravitational initial data on which further evolution work can
be based.

## 6. Reproducibility and evidence

The [raw archive](../experiments/closure_ledger/runs/20260915_mouth_momentum/probe.json.gz)
retains every primary scalar solution, momentum solve, physical metric and
extrinsic-curvature sample, coordinate residual, section diagnostic and
control. The [gate report](../experiments/closure_ledger/runs/20260915_mouth_momentum/probe.md)
records all eight results and separate scoped verdicts. Its decompressed
SHA256 is `bb9ad63c3468e4795d5147a512d1e045df14089075bae66d42937ea2d5bce89b`.
An independent complete rerun reproduces that compressed archive byte for
byte on this runtime. Source hashes are retained in the run's manifest.

The targeted suite passes 27 tests. It includes missing/duplicated cases,
nonfinite and changed solutions, missing/altered K, changed curvature and
momentum evidence, altered ODE source data, missing/unknown gates, forged
seams and incorrect parameters. Actual CLI failure paths overwrite a stale
affirmative report. Gate scoring recomputes PDE residuals and crosschecks
the archived physical data against the stored solutions rather than trusting
passing flags. Additional validation checks did not change the equations,
grids, thresholds or numerical archive.

```sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.mouth_momentum_probe \
  --output-dir experiments/closure_ledger/runs/20260915_mouth_momentum
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_mouth_momentum.py
```

The construction uses the standard conformal method; see the references
in the public specification. The contribution is a repository implementation
and verified benchmark, not a novel general wormhole existence theorem.

### Python 3.10 compatibility correction (2026-09-16)

Repository-wide CI initially stopped during collection because the probe used
Python 3.11 starred-subscript syntax. Explicit tuple construction restores
the declared Python 3.10 support. A grammar regression was added; 28 targeted
tests pass. No equation, tolerance, archived solution or scientific gate
changed. The manifest retains original and corrected source hashes.

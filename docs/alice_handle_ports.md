# What crosses an Alice handle: result

Probe of the pre-registration `docs/alice_handle_ports_prereg.md`, frozen at
`36ba9c2`. Archive: `experiments/closure_ledger/runs/20261118_alice_handle_ports/`
(`probe.json`, `probe.md`). Code: `geometrodynamics/bulk/alice_handle_ports.py`,
`experiments/closure_ledger/alice_handle_ports_probe.py`. Tests:
`tests/test_alice_handle_ports.py`.

## Verdict

**`COUNTEREXAMPLE`**, on the split in P4. The frozen rules make P4's split a
structural statement, and it fails at both registered resolutions. P5's
identity tolerance is also missed at both, and the miss grows with `N`, so it
is not `UNRESOLVED_NUMERICALLY` either. Everything else holds at both
resolutions, and P7's geometry and spectrum checks pass.

| statement | N = 2000 | N = 4000 | gate |
|---|---|---|---|
| P1 `m = 0` gluings identical | pass (0) | pass (0) | ≤ 1e-10 |
| P2 torus degeneracy | pass (1.3e-14) | pass (1.9e-14) | ≤ 1e-10 |
| P2 torus `L` conserved | pass (1.7e-10) | pass (2.1e-10) | ≤ 1e-9 |
| P3 Klein = `c = 0` ∪ `c = π/m` | pass (1.2e-14) | pass (3.0e-14) | ≤ 1e-10 |
| **P4 lowest six pairs split** | **FAIL (1.01e-8)** | **FAIL (1.01e-8)** | > 1e-6 |
| P4 first sign change at `π/(2Δ)` | pass (0.32 of bound) | pass (0.32 of bound) | ≤ `Δ/ω̄` |
| **P5 identity** | **FAIL (7.0e-7)** | **FAIL (8.6e-6)** | ≤ 1e-8 |
| P5 `∫Φ_B dt` non-vacuous | pass (0.99) | pass (0.99) | ≥ 0.1 |
| P6 Klein `L_R` = −torus `L_R` | pass (7.8e-8) | pass (1.7e-8) | sign |
| P6 tolerance and size | pass (1.6e-7) | pass (3.4e-8) | ≤ 1e-3 |

The number in brackets is the worst case over all runs at that resolution: the
smallest `Δ/ω̄` for the split, the largest timing error as a fraction of its
bound, the largest identity residual relative to `|L(0)|`, the smallest
`max|∫Φ_B dt| / |L(0)|`, and the largest antisymmetry residual.

## What failed, and what it means

**P4's split.** Every Klein doublet is split, and the splitting is resolved:
the smallest one, `a = 0.3`, `m = 3`, is `1.01082e-8` at `N = 2000` and
`1.01134e-8` at `N = 4000`, so it is converged to `5e-4` relative and lies far
above eigensolver precision. But it is below the frozen threshold `1e-6`, so
the statement "the lowest six pairs are split, `Δ > 1e-6 ω̄`" is false for
this geometry. Four of the six `a = 0.3`, `m = 3` pairs fall below it, and
no other pair does.

The threshold came from the expectation, which was stated as "not a gate",
and the expectation is wrong by a growing factor:

| `a` | `m` | measured `Δ/ω̄`, lowest pair | expectation `exp(−2m asinh(cot a))` |
|---|---|---|---|
| 0.3 | 1 | 3.16e-3 | 2.3e-2 |
| 0.3 | 2 | 5.51e-6 | 5.2e-4 |
| 0.3 | 3 | 1.01e-8 | 1.2e-5 |
| 0.5 | 1 | 2.62e-2 | 6.5e-2 |
| 0.5 | 2 | 3.79e-4 | 4.3e-3 |
| 0.5 | 3 | 5.72e-6 | 2.8e-4 |

The measured splitting does fall exponentially in `m`, but faster than
expected. The log of the ratio between successive `m` is 6.35 and 6.30 at
`a = 0.3`, against an expected 3.78. At `a = 0.5` it is 4.24 and 4.19,
against 2.73. That makes the exponent 1.5 to 1.7 times the hand estimate.
Where that factor comes from is not derived here.

The reversal law that P4 is about survives where the split is this small. In
all 36 pair runs at each resolution, `L_band + L_neck` changes sign first at
`π/(2Δ)` to within the frozen bound. That includes `Δ/ω̄ = 1.0e-8`, where the
sign change comes at `ω̄t ≈ 1.6e8`. The worst timing error is 0.32 of
its bound. The verdict does not change because of this: the split statement
is false as frozen.

**P5's identity.** The bookkeeping identity
`L(t) − L(0) + 2∫₀ᵗ Φ_B dt = 0` holds to 1e-12 on the P6 packet runs, which
last a time of order 1. On the P4 runs it has to hold until `t = π/Δ`, which
is `ω̄t ≈ 3e8` for the smallest split. Its residual there grows roughly
linearly in `t` (the archive keeps the series as `identity_series`) and
reaches 7.0e-7 at `N = 2000` and 8.6e-6 at `N = 4000`. Over the same runs,
the residual `‖Kv − λMv‖/‖Kv‖` of the two eigenvectors goes from 1.0e-9 to
4.2e-9.

The residual tracks run length and eigenvector accuracy, not the gluing. The
likeliest reading is that the mode sum's eigenvectors are not exact enough to
carry the identity over `ω̄t ~ 1e8`. It is not a measured leak of angular
momentum, but it is a missed tolerance, and it is reported as missed. The
identity is not established to 1e-8 for the longest runs. The non-vacuity
half of P5 holds: in every run, `|∫Φ_B dt|` reaches 99% of `|L(0)|` or more,
so the handle does carry essentially all of `L(0)` across rim B.

## What held

* **P1.** For `m = 0` the torus and Klein spectra, and the energy that crosses
  each rim, are identical. They are the same discrete problem.
* **P2.** Torus gluing with `c = 0` and `c = 1` keeps every `m ≥ 1` level
  doubly degenerate and conserves `L_band + L_neck`.
* **P3.** The Klein spectrum in each `m ≥ 1` is the union of the `c = 0` and
  `c = π/m` torus spectra, with the cosine sector carrying the first and the
  sine sector the second.
* **P6, the reversal.** A packet sent through the handle arrives at the far
  side of the band with its angular momentum reversed relative to the torus
  run. In the window before anything can arrive round the sphere,
  `L_R(Klein) + L_R(torus)` is at most 1.6e-7 (`N = 2000`) and 3.4e-8
  (`N = 4000`) of the largest `|L_R|`, while their difference is 2.000000. The packet's peak `|L_R|` is 0.91 (`m = 1`)
  and 0.75 (`m = 2`) of `L(0)`.
* **P7.** The #338 mesh agrees with the separable solver with no fitted
  parameter. The rim azimuth difference is constant across glued torus-mesh
  vertices (spread 1e-15) and equals π. The sum is constant for the Klein
  mesh (spread 1e-15). The lowest 20 nonzero `ω²` agree to 0.33% at
  `n_sphere = 6000`, `n_ring = 48` and to 0.18% at 12000/64, for both gluings.
  The largest `m` among them is 3, below the bound of 8.

## Implementation choices the freeze left open

* **Nodes and mass.** Each loop segment, band and neck, gets uniform nodes,
  with the rims as nodes. The mass matrix is lumped, so that the discrete flux
  `Φ_B` makes the semi-discrete continuity equation exact, as P5 requires.
* **Eigenvalues and evolution.** Eigenvalues are Rayleigh quotients of the
  computed eigenvectors. Evolution is an exact mode sum with long-double
  phases.
* **P4 sign search.** P4's "first changes sign" is found by a certified
  search. Every interval is either cleared, using a fast/slow split, a
  Lipschitz bound or a curvature bound, or bisected to a sign change. The
  search raises an error rather than skip an interval it can do neither with.
* **Where P5 is checked.** P5's identity is checked on the P4 runs and also on
  the P6 Klein runs.
* **P7 sizes.** P7 counts the separable levels with `m` up to 8 at
  `N = 2000`. The torus mesh's twist was measured, not assumed: it is π.

## Corrections after review

The draft reviewed on PR #355 before the full run had two defects. Both were
fixed in `caf3912`, before these numbers were produced:

* **The classifier.** It returned `COUNTEREXAMPLE` from a single resolution,
  and `UNRESOLVED_NUMERICALLY` when a miss grew. It now needs both registered
  resolutions, needs a miss to shrink, and allows only the one further
  doubling. It returns `NO_VERDICT`, a label outside the freeze's vocabulary,
  only in cases the freeze does not name; that did not happen here.
* **The sign search.** It could step over a shallow first crossing, for
  example `1 − 1e-7 − cos(t − 0.005)`. It is now the certified search
  described above, with that case as a test.

The archive was assembled from stages cached by separate runs, each stage
computed whole by the current code. A full `N = 4000` stage is longer than the
session's two-hour limit for one command.

## Scope

As frozen. This is a scalar on a fixed 2D surface: it is not a GR solution
and not on `S³`. The reversal is classical, a consequence of the broken
`U(1)`, and is not claimed to be quantum. Because the verdict is not
`CONFIRMED`, the second gate in the pre-registration, a dynamical reflection
axis, is not opened by this result.

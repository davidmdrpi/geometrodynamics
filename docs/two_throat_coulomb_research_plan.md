# Two-throat Coulomb probe — research plan

A falsification test flagged in `docs/THESIS.md`: the finite-separation
two-throat Coulomb law on `S³` is listed as "a near-term falsification
test, not a demonstrated result." This probe runs it.

The PR #45 / PR #46 exchange-kernel thread showed that the `S³` scalar
Green function reduces in the flat-space limit to the Coulomb potential
`1/(4π d)` and Fourier-transforms to the photon propagator `1/q²`
(scattering / momentum-space picture). This probe tests the
**static / position-space** picture directly: two charged BAM throat
mouths on `S³`, interacting through the same `S³` Green function,
should reproduce

  (a) the Coulomb potential `V ∝ 1/r` in the flat-space limit, and
  (b) the inverse-square force `F ∝ 1/r²` in the flat-space limit,

and should exhibit the genuine compact-`S³` modifications at finite
separation.

## Setup

A point charge `q` on `S³` of radius `R` sources the zero-mean
potential

```
V(ψ) = q · G(ψ),    G(ψ) = ((π − ψ) cot ψ − ½) / (4π² R)
```

(from `geometrodynamics.transaction.s3_geometry.s3_green_potential`,
the same kernel as PR #45). The zero-mean subtraction is the
compact-manifold neutralizing background (Gauss's law on a closed
manifold requires total charge zero; the uniform `−q/V_S³` background
makes the Poisson equation solvable).

The interaction energy of two charges `q₁, q₂` at geodesic angle ψ:

```
U(ψ) = q₁ q₂ G(ψ)
```

The geodesic force (gradient along the proper distance `s = R·ψ`):

```
F(ψ) = −dU/ds = −(q₁ q₂ / R) · dG/dψ
     =  q₁ q₂ · |dG/dψ| / R         (magnitude; sign below)
     =  q₁ q₂ · [(π − ψ) + sin ψ cos ψ] / (4π² R² sin²ψ)
```

using `dG/dψ = −[(π − ψ) + sin ψ cos ψ] / (4π² R sin²ψ)`. The repo
field kernel `s3_green_field_kernel(ψ) = |dG/dψ|/R` provides `|E|`.

## Predictions

### P1. Coulomb potential in the flat-space limit

For small ψ with `r = R·ψ`: `V(ψ) → q / (4π r)`. The `1/r` Coulomb
potential (Heaviside–Lorentz units).

### P2. Inverse-square force in the flat-space limit

For small ψ: `F(ψ) → q₁ q₂ / (4π r²)`. The `1/r²` inverse-square law.

### P3. Exact S³ force form

```
F(ψ) = q₁ q₂ · N(ψ) / (4π² R² sin²ψ),    N(ψ) = (π − ψ) + sin ψ cos ψ
```

The leading behaviour is `1/sin²ψ` (the THESIS claim
`F ∝ 1/sin²ψ`), but the exact force is modulated by `N(ψ)`, which
runs from `N(0) = π` (Coulomb regime) to `N(π) = 0`. The probe tests
whether the THESIS `1/sin²ψ` is exact or leading-order:

  - **Confirmed as leading**: `F · sin²ψ → const` only near ψ = 0;
    the full `N(ψ)` modulation is the compact-`S³` correction.
  - The modulation is a genuine prediction, not a defect: it
    encodes the antipodal image (P4).

### P4. Antipodal equilibrium

`N(π) = 0`, so `F(π) = 0`. A charge on `S³` has its field lines
converge at the antipode; the force on a test charge placed exactly
at the antipode vanishes by symmetry. This is the compact-`S³`
signature absent in flat space — a charge "feels" its own antipodal
image.

### P5. Gauss's law on S³

The electric flux through the geodesic 2-sphere at colatitude ψ
(area `4π R² sin²ψ`) equals the enclosed charge:

```
Φ(ψ) = |E(ψ)| · 4π R² sin²ψ = q · N(ψ)/π
     = q · [(π − ψ) + sin ψ cos ψ] / π
```

This equals the **point charge plus the neutralizing background
within the cap**:

```
Q_enclosed(ψ) = q − q · (cap volume fraction)
              = q − q · (ψ − sin ψ cos ψ)/π
              = q · [(π − ψ) + sin ψ cos ψ]/π
```

so `Φ(ψ) = Q_enclosed(ψ)` exactly. Flux falls from `q` at the source
to `0` at the antipode, consistent with the uniform background. This
is a strong internal-consistency test of the whole picture.

### P6. Attraction / repulsion sign

Like charges (`q₁ q₂ > 0`) repel (force pushes them to larger ψ);
opposite charges attract. Verify the sign of `−dU/ds`.

## Two-throat-mouth realization

Each BAM particle is a throat-pair: a front mouth and a back mouth
at the `S³` antipode. The probe places two actual mouths using the
`s3_geometry` machinery (`antipode4`, `geo4`) and computes the
interaction via the Green function, verifying it matches the analytic
`V(ψ)`. The throat-pair (dipole) structure — front `+q`, back at
antipode — is tested as an extension: the antipodal-image effect (P4)
is automatically built into a single throat-pair.

## Tests

  T1. **Potential = S³ Green function**: `V(ψ) = q·G(ψ)` from the repo.
  T2. **Coulomb potential (flat limit)**: `V(ψ)·r → q/(4π)` as ψ → 0.
  T3. **Force from field kernel**: `F(ψ) = q₁q₂·|dG/dψ|/R` matches
      the analytic `q₁q₂·N(ψ)/(4π²R²sin²ψ)`.
  T4. **Inverse-square force (flat limit)**: `F(ψ)·r² → q₁q₂/(4π)`.
  T5. **Exact S³ force / 1/sin²ψ claim**: characterise `N(ψ)`;
      confirm `F·sin²ψ → q₁q₂·π/(4π²R²)` only near ψ = 0.
  T6. **Antipodal equilibrium**: `F(π) = 0` (N(π) = 0).
  T7. **Gauss's law on S³**: `Φ(ψ) = Q_enclosed(ψ)` to machine
      precision across ψ ∈ (0, π).
  T8. **Attraction/repulsion sign**: like charges repel, opposite
      attract.
  T9. **Two-mouth placement realization**: place two mouths at
      explicit `S³` points; verify `geo4`-separation interaction
      matches `V(ψ)`.

## Verdict structure

  - **COULOMB_REPRODUCED**: P1, P2 (flat-limit Coulomb + inverse
    square) hold to machine precision; P5 (Gauss's law) exact; P3/P4
    characterise the compact-`S³` modifications (`N(ψ)` modulation,
    antipodal equilibrium). The two-throat `S³` Green-function
    interaction reproduces electrostatics, with the THESIS `1/sin²ψ`
    claim confirmed as the leading behaviour and refined by the
    antipodal-image modulation.

  - **COULOMB_FALSIFIED**: the flat-limit Coulomb / inverse-square
    law fails, or Gauss's law is violated. Would falsify the
    "charge = throat mouth on `S³`" identification.

## What this leaves open

  - **Dynamical / radiative**: this is the static potential. Moving
    charges (radiation, magnetic fields) need the time-dependent
    `S³` Green function — a separate probe.
  - **Throat-pair internal structure**: whether the back mouth
    carries `−q` (dipole) or `+q` (monopole pair) affects the
    long-range field; the probe tests both but a first-principles
    determination from the Hopf charge is open.
  - **Self-energy**: the `ψ → 0` divergence of `G(ψ)` is the usual
    point-charge self-energy; regularisation by the finite throat
    radius is open.

## Cross-references

  - `geometrodynamics/transaction/s3_geometry.py` —
    `s3_green_potential`, `s3_green_field_kernel`, `antipode4`, `geo4`.
  - PR #45: `bam_exchange_kernel_probe.py` — same Green function,
    momentum-space (scattering) picture.
  - `docs/THESIS.md` — flags this as a falsification test.
  - `experiments/closure_ledger/two_throat_coulomb_probe.py` —
    this probe.

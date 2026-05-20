# Two-throat Coulomb probe

**Run:** 2026-05-20T06:17:58+00:00

Falsification test (flagged in docs/THESIS.md): two charged BAM throat mouths on S³, interacting through the S³ Green function, vs the Coulomb potential and inverse-square force law.

## Electrostatics on S³

```
potential:    V(ψ) = q·G(ψ), G(ψ) = ((π−ψ)cot ψ − ½)/(4π²R)
force (exact): F(ψ) = q₁q₂·N(ψ)/(4π²R²sin²ψ), N(ψ)=(π−ψ)+sinψcosψ
flat limits:  V → q/(4π r), F → q₁q₂/(4π r²)
Gauss law:    Φ(ψ) = q·N(ψ)/π = Q_enclosed(ψ)
```

## Test summary

| # | Test | Key finding | PASS? |
|---|---|---|---|
| T1 | `T1_potential_is_s3_green_function` | V = q·G(ψ), max diff = 2.78e-17 | **PASS** |
| T2 | `T2_coulomb_potential_flat_limit` | V·r → q/(4π), residual = 3.80e-07 | **PASS** |
| T3 | `T3_force_from_field_kernel` | F kernel = analytic, max diff = 8.88e-16 | **PASS** |
| T4 | `T4_inverse_square_force_flat_limit` | F·r² → q₁q₂/(4π), residual = 2.65e-12 | **PASS** |
| T5 | `T5_exact_s3_force_form` | 1/sin²ψ leading: True; N(ψ) modulation: True | **PASS** |
| T6 | `T6_antipodal_equilibrium` | N(π)→0 (F→0 at antipode): -1.22e-16 | **PASS** |
| T7 | `T7_gauss_law_on_s3` | Φ = Q_enclosed, max diff = 1.11e-16 | **PASS** |
| T8 | `T8_attraction_repulsion_sign` | like repels: True; opposite attracts: True | **PASS** |
| T9 | `T9_two_mouth_placement_realization` | V matches, max diff = 0.00e+00; back-mouths at π: True | **PASS** |

## T2: Coulomb potential (flat-space limit)

| ψ | r = Rψ | V | V·r | q/(4π) | residual |
|---:|---:|---:|---:|---:|---:|
| 1e-05 | 1e-05 | 7.9577e+03 | 0.079577 | 0.079577 | -3.80e-07 |
| 1e-04 | 1e-04 | 7.9574e+02 | 0.079574 | 0.079577 | -3.80e-06 |
| 1e-03 | 1e-03 | 7.9539e+01 | 0.079539 | 0.079577 | -3.80e-05 |
| 1e-02 | 1e-02 | 7.9195e+00 | 0.079195 | 0.079577 | -3.83e-04 |
| 5e-02 | 5e-02 | 1.5522e+00 | 0.077612 | 0.079577 | -1.97e-03 |
| 1e-01 | 1e-01 | 7.5521e-01 | 0.075521 | 0.079577 | -4.06e-03 |

## T4: Inverse-square force (flat-space limit)

| ψ | r = Rψ | F | F·r² | q₁q₂/(4π) | residual |
|---:|---:|---:|---:|---:|---:|
| 1e-05 | 1e-05 | 7.9577e+08 | 0.079577 | 0.079577 | +2.65e-12 |
| 1e-04 | 1e-04 | 7.9577e+06 | 0.079577 | 0.079577 | +2.65e-10 |
| 1e-03 | 1e-03 | 7.9577e+04 | 0.079577 | 0.079577 | +2.65e-08 |
| 1e-02 | 1e-02 | 7.9580e+02 | 0.079580 | 0.079577 | +2.64e-06 |
| 5e-02 | 5e-02 | 3.1857e+01 | 0.079642 | 0.079577 | +6.42e-05 |
| 1e-01 | 1e-01 | 7.9826e+00 | 0.079826 | 0.079577 | +2.49e-04 |

## T5: Exact S³ force form / 1/sin²ψ claim

| ψ | F | F·sin²ψ | N(ψ) | F·sin²ψ / Coulomb-const |
|---:|---:|---:|---:|---:|
| 0.0100 | 7.9580e+02 | 0.079577 | 3.141592 | 1.000000 |
| 0.1000 | 7.9826e+00 | 0.079561 | 3.140927 | 0.999788 |
| 0.5000 | 3.3748e-01 | 0.077570 | 3.062328 | 0.974769 |
| 1.0000 | 9.2877e-02 | 0.065764 | 2.596241 | 0.826409 |
| 1.5708 | 3.9789e-02 | 0.039789 | 1.570796 | 0.500000 |
| 2.0000 | 2.3381e-02 | 0.019332 | 0.763191 | 0.242931 |
| 2.5000 | 1.1466e-02 | 0.004107 | 0.162131 | 0.051608 |
| 3.0000 | 2.3975e-03 | 0.000048 | 0.001885 | 0.000600 |

1/sin²ψ confirmed as leading: **True**; N(ψ) modulation present (compact-S³): **True**.

## T6: Antipodal equilibrium

| ψ | π − ψ | N(ψ) | force |
|---:|---:|---:|---:|
| 2.6416 | 5e-01 | 7.9265e-02 | 8.7353e-03 |
| 3.0416 | 1e-01 | 6.6533e-04 | 1.6909e-03 |
| 3.1316 | 1e-02 | 6.6665e-07 | 1.6887e-04 |
| 3.1406 | 1e-03 | 6.6667e-10 | 1.6887e-05 |

N near antipode = `-1.22e-16` → force vanishes at ψ = π. A charge on compact S³ feels its own antipodal image; the antipode is a force-free equilibrium absent in flat space.

## T7: Gauss's law on S³

| ψ | flux = E·area | Q_enclosed (point + bg cap) | diff |
|---:|---:|---:|---:|
| 0.1000 | 0.999788 | 0.999788 | 1.11e-16 |
| 0.5000 | 0.974769 | 0.974769 | 1.11e-16 |
| 1.0000 | 0.826409 | 0.826409 | 1.11e-16 |
| 1.5708 | 0.500000 | 0.500000 | 0.00e+00 |
| 2.0000 | 0.242931 | 0.242931 | 0.00e+00 |
| 2.5000 | 0.051608 | 0.051608 | 6.25e-17 |
| 3.0000 | 0.000600 | 0.000600 | 2.02e-17 |

Flux falls from q (source) to 0 (antipode) as the neutralizing background charge is enclosed; Φ(ψ) = Q_enclosed(ψ) exactly.

## T9: Two-mouth placement realization

| χ | geodesic ψ | U interaction | U analytic | diff | back-mouth dist |
|---:|---:|---:|---:|---:|---:|
| 0.20 | 0.2000 | 0.354911 | 0.354911 | 0.00e+00 | 3.141593 |
| 0.50 | 0.5000 | 0.109817 | 0.109817 | 0.00e+00 | 3.141593 |
| 1.00 | 1.0000 | 0.022167 | 0.022167 | 0.00e+00 | 3.141593 |
| 1.50 | 1.5000 | -0.009716 | -0.009716 | 0.00e+00 | 3.141593 |
| 2.00 | 2.0000 | -0.025899 | -0.025899 | 0.00e+00 | 3.141593 |
| 2.50 | 2.5000 | -0.034420 | -0.034420 | 0.00e+00 | 3.141593 |

## Verdict

**COULOMB_REPRODUCED.** COULOMB LAW REPRODUCED. Two charged BAM throat mouths on S³, interacting through the S³ Green function, reproduce electrostatics:
  (P1) Coulomb potential V → q/(4π r) in the flat-space limit (T2);
  (P2) inverse-square force F → q₁q₂/(4π r²) in the flat-space limit (T4);
  (P5) Gauss's law Φ(ψ) = Q_enclosed(ψ) holds exactly on compact S³, with flux falling from q at the source to 0 at the antipode as the neutralizing background is enclosed (T7).
The THESIS finite-separation claim F ∝ 1/sin²ψ is CONFIRMED as the leading short-range behaviour and REFINED: the exact force F(ψ) = q₁q₂·N(ψ)/(4π²R²sin²ψ) carries the modulation N(ψ) = (π−ψ)+sinψcosψ (T5), which encodes the compact-S³ antipodal image — the force vanishes exactly at the antipode (T6), an equilibrium absent in flat space. Like charges repel, opposite charges attract (T8). The two-mouth placement realization confirms the geometric picture, with each particle's back-mouth at geodesic distance π (T9). The "charge = throat mouth on S³" identification is NOT falsified — it reproduces Coulomb electrostatics with calculable compact-manifold corrections.

## What this leaves open

- **Dynamical / radiative**: this is the static potential. Moving charges (radiation, magnetic fields) need the time-dependent S³ Green function — a separate probe.
- **Throat-pair internal structure**: whether the back mouth carries −q (dipole) or +q (monopole pair) affects the long-range field; a first-principles determination from the Hopf charge is open.
- **Self-energy**: the ψ → 0 divergence of G(ψ) is the usual point-charge self-energy; regularisation by the finite throat radius is open.

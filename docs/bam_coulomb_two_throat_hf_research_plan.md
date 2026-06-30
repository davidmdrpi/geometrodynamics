# The BAM Coulomb-photon kernel for the two-throat HF: replacing the Yukawa stand-in (PR #190)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## Replacing the stand-in with the real photon

PR #187 and #189 used a screened-photon (Yukawa) interaction as a **regulated
stand-in** for the BAM throat-fibre exchange. This probe replaces it with the
genuine BAM Coulomb-photon kernel — the **unscreened** Coulomb
`V(d) = 1/(4πd)` (real space) ⟷ `1/q²` (Fourier), the photon propagator BAM
derives from the throat-fibre exchange geometry (the #42–#44 result) — and
recomputes the two-throat direct + exchange energies.

The kernel is the **flat-space limit of the BAM S³ scalar Green function**

```
G(ψ) = ((π − ψ) cot ψ − ½) / (4π² R)   (the repo's s3_green_potential),
```

which reproduces the Euclidean `1/(4πs)` singularity near the source. On the
local weak-field patch the throats see the unscreened Coulomb, with S³
curvature corrections `O(1/R²)` carried by `G`. The interaction is regulated
**properly for an isolated system** (the Hockney zero-padded Coulomb solver),
**not** by the ad-hoc Yukawa screening.

## The kernel

- The isolated-system Coulomb of a unit point source reproduces `1/(4πd)` to
  machine precision (`Φ = {1.25: 0.0637, 2.5: 0.0318, 3.75: 0.0212} = 1/(4πd)`),
  whose 3D Fourier transform is `1/q²` — the photon propagator in Feynman gauge
  (#42–#44).
- The S³ Green function's near-source limit is the same Coulomb coefficient:
  `G·4πs = 0.957 → 1` (`s = Rψ`, just above the function's `eps = 0.08`
  regulator floor).

## The regulator

The unscreened Coulomb is solved by the **Hockney zero-padded convolution**
(the density padded to a 2× box so periodic images do not interact, convolved
with the free-space `1/(4πr)` Green function). Validated against the analytic
Gaussian Coulomb self-energy:

```
U = 0.0167  vs  exact ½·(1/4π)·1/(σ√π) = 0.0167 ,   ratio 0.9992  (~0.08%) .
```

So the kernel is the true unscreened Coulomb — the Yukawa screening parameter
`μ` is gone, replaced by a numerical open-boundary regulator with no spurious
long-range cutoff.

## The recomputed energies (on the #180 orbitals)

| R | direct `J` | exchange `K_ex` |
|---:|---:|---:|
| 0.0 | 0.0627 | 0.0627 |
| 1.0 | 0.0536 | 0.0383 |
| 2.0 | 0.0377 | 0.0096 |
| 3.0 | 0.0265 | 0.0012 |
| 4.0 | 0.0200 | 0.0001 |
| 6.0 | 0.0133 | 0.0000 |

The direct `J(R)` is now correctly **long-ranged** — `J(6) = 0.0133 ≈
1/(4π·6) = 0.01326` (ratio `1.003`), the point-charge Coulomb tail — unlike
the Yukawa stand-in's exponential decay. The exchange `K_ex(R)` stays
**short-ranged** (set by the overlap density), so far-separated throats feel
the Coulomb direct field but not the exchange — the physically correct
long-range structure the screened stand-in lacked.

## The #187 physics survives

With the overlap-normalized `E±(R) = (J ± K_ex)/(1 ± S²)`:

| R | `S` | `E₊` (boson) | `E₋` (fermion) |
|---:|---:|---:|---:|
| 1.0 | 0.79 | 0.0565 | 0.0409 |
| 2.0 | 0.41 | 0.0405 | 0.0337 |
| 3.0 | 0.15 | 0.0271 | 0.0258 |
| 4.0 | 0.05 | 0.0200 | 0.0199 |

For the **repulsive** photon the antisymmetric (Pin⁻) branch `E₋` lies below
the symmetric `E₊` at every finite separation — the **fermion-lower** result
of #187, established there with the Yukawa stand-in, holds with the real
unscreened photon. At coincidence (`R = 0`) `J = K_ex = 0.0627` (numerator
`0`) and `S → 1`, so the antisymmetric state is the **zero vector**
(Pauli-forbidden). The statistics are a property of the geometry (the Pin⁻
sign + the overlap structure), not of the interaction's screening.

## Honest scope

- The kernel is now the genuine BAM photon — the **flat Coulomb limit** of the
  S³ Green function; the `O(1/R²)` S³ curvature corrections are carried by `G`
  but not applied here (the weak-field local patch).
- The unscreened Coulomb is regulated by the Hockney **open-boundary** solver
  (validated to ~0.08% on the Gaussian self-energy), **not** physical
  screening (the Yukawa `μ` is gone).
- The orbitals are the **rigid** #180 throat-solitons; the self-consistent SCF
  with the Coulomb kernel (the #189 relaxation with the real photon) is the
  follow-up.
- Energies in code units; weak-field. The upshot: the Yukawa was a faithful
  short-range stand-in, and replacing it with the unscreened photon leaves the
  #187 statistics intact while making the direct channel correctly
  long-ranged.

## Reproduce

```bash
python -m experiments.closure_ledger.bam_coulomb_two_throat_hf_probe
# Verdict: BAM_COULOMB_PHOTON_KERNEL_REPLACES_YUKAWA_STANDIN_LONG_RANGED_DIRECT_TWO_THROAT_HF_PHYSICS_ROBUST
```

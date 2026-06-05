# Antipodal-horizon exchange kernel from the throat boundary data (PR #135)

**Naming note.** PRs #42–#44 already built a "BAM exchange kernel" for the
**gauge** sector — the photon propagator `1/q²` from the S³ scalar Green
function (`bam_exchange_kernel_probe.py`). This PR builds the complementary
**matter-sector** exchange kernel: the two-point Green's function / resolvent of
the matter cavity operator (#116) with the **antipodal horizon** boundary data
(#129).

The geometric throat arc fixed the boundary data the antipodal horizon imposes:
the l-parity condition (even-l Neumann, odd-l Dirichlet, #129), a unitary
mirror, with a real normal-mode spectrum (#130). This PR builds the object that
boundary data defines — the matter exchange kernel: the amplitude for a quantum
to be exchanged between two points through the antipodal throat, assembled as a
sum over the stable modes #130 supplies.

## The kernel = the cavity resolvent with antipodal boundary data

For each angular channel `l`, the exchange kernel is the resolvent of the matter
cavity operator `H_l = −d²/dr*² + V_l` (`V_l = f[l(l+2)/r² + 3rs²/r⁴]`, #116)
with the antipodal boundary data of #129 — Neumann `ψ'(throat)=0` for even `l`,
Dirichlet `ψ(throat)=0` for odd `l`, Dirichlet shell wall at `R_OUTER`:

```
K_l(r, r'; ω) = ⟨ r | (H_l − ω²)^{−1} | r' ⟩ .
```

The operator is exactly self-adjoint (symmetric matrix, `max|H − Hᵀ| = 0`).

## Spectral representation: exchange of the stable modes

Because the antipodal operator is self-adjoint (#129/#130), it has a complete
real eigenbasis `H_l ψ_n = ω_n² ψ_n`, and the kernel is the **mode sum**

```
K_l(r, r'; ω) = Σ_n ψ_n(r) ψ_n(r') / (ω_n² − ω²) .
```

| l | poles ω² (kernel) | mode-sum vs resolvent |
|---:|---|---:|
| 0 | 1.363, 10.585, 28.867 | 1e-14 |
| 1 | 5.269, 19.046, 41.991 | 7e-16 |
| 2 | 2.018, 11.587, 29.880 | 2e-15 |
| 3 | 6.728, 20.581, 43.523 | 8e-16 |

The **poles are the real normal-mode spectrum of #130** — so the exchange kernel
is a propagator assembled as a sum over the stable exchanged modes, with no
decaying contribution (mode sum = matrix resolvent to ~1e-14).

## Reciprocity

The antipodal boundary data makes `H_l` self-adjoint, so the kernel is
symmetric, `K_l(r, r') = K_l(r', r)` (reciprocity; `|K − Kᵀ|/|K| ~ 1e-14`).

## Unitary vs lossy exchange — the boundary data decides

The antipodal boundary data (real BC) ⟹ Hermitian `H_l` ⟹ **real poles** ⟹ an
undamped, **unitary** exchange kernel. The absorbing horizon (ingoing BC) would
give a non-Hermitian operator with **complex poles** ⟹ a decaying, **lossy**
kernel (#130). So the antipodal horizon boundary data is exactly what makes the
matter exchange kernel unitary — the propagator-level face of the unitary mirror
(#129) and the global CPT/unitarity (#64).

## Angular antipodal-parity grading

The full kernel factorises `K(x, x') = Σ_l K_l(r, r'; ω) · C_l(Ω·Ω')` with `C_l`
the S³ zonal harmonic. Under the throat ↔ antithroat (antipodal) exchange
`Ω' → AΩ'`, `C_l(−Ω·Ω') = (−1)^l C_l(Ω·Ω')` (degree-l harmonics, #129/#134), so
each `l`-channel of the exchange kernel carries the antipodal sign `(−1)^l` —
even-`l` channels symmetric, odd-`l` channels antisymmetric under the C-swap
(#63). The exchange kernel is parity-graded by the same `(−1)^l` that fixed the
boundary condition.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | build the matter exchange kernel from the antipodal data (#129) |
| T2 | kernel = resolvent | `K_l = (H_l − ω²)^{-1}`, self-adjoint |
| T3 | spectral rep | `K_l = Σ ψψ/(ω_n²−ω²)`; poles = #130 spectrum (~1e-14) |
| T4 | reciprocity | `K_l(r,r') = K_l(r',r)` (~1e-14) |
| T5 | unitary vs lossy | antipodal real poles (unitary) vs absorbing complex (lossy, #130) |
| T6 | parity grading | angular channels carry `(−1)^l` (#129/#134) |
| T7 | scope | free/one-loop kernel; interacting kernel / normalisation open |
| T8 | assessment | `ANTIPODAL_HORIZON_EXCHANGE_KERNEL_UNITARY_RECIPROCAL_PARITY_GRADED` |

## Established and open

  - **Established (BAM-native):** the matter exchange kernel is the cavity
    resolvent with the antipodal horizon boundary data (#129) — a reciprocal,
    unitary (real-pole) two-point kernel whose spectral representation is a sum
    over the stable modes (#130) and whose l-channels are antipodal-parity-
    graded `(−1)^l`.

  - **Does not / open:** this is the FREE / one-loop kernel on the fixed
    antipodal background (the S_BAM fluctuation-measure propagator, #115–#122);
    it does **not** include the interacting / multi-loop kernel (vertices,
    self-energy) or fix the absolute normalisation; the bulk-scale (#133) and
    flavor (#134) residuals stand. (The gauge-sector photon kernel `1/q²` is the
    separate PR #42–#44 probe.)

## Cross-references

  - `docs/null_throat_boundary_conditions_research_plan.md` — #129, the
    antipodal l-parity boundary data this kernel is built from.
  - `docs/antipodal_vs_absorbing_qnm_research_plan.md` — #130, the real spectrum
    that are the kernel's poles (and the absorbing complex-pole counterfactual).
  - `docs/tangherlini_fluctuation_determinant_research_plan.md` — #116, the
    cavity operator.
  - `bam_exchange_kernel_probe.py` (PRs #42–#44) — the complementary
    gauge-sector exchange kernel (photon `1/q²` from the S³ Green function).

## Run

```
python -m experiments.closure_ledger.antipodal_horizon_exchange_kernel_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_antipodal_horizon_exchange_kernel_probe/`.
Expected verdict:
`ANTIPODAL_HORIZON_EXCHANGE_KERNEL_UNITARY_RECIPROCAL_PARITY_GRADED`, 8/8 PASS.

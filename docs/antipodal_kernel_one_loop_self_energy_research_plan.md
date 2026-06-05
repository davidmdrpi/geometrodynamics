# One-loop self-energy audit for the antipodal matter kernel (PR #136)

PR #135 built the matter-sector exchange kernel — the **free** propagator
`G_0 = (H_l − ω²)^{−1}` of the cavity operator with the antipodal horizon
boundary data (#129) — and showed it is unitary, with real poles (the stable
#130 spectrum). PR #135 listed the **interacting / self-energy** kernel as the
lead open item. This PR audits the leading interacting correction: the
**one-loop self-energy** `Σ`, and asks whether it preserves the tree-level
stability and unitarity.

## The Dyson-dressed propagator

The self-energy `Σ(s)` (`s = ω²`) dresses the free propagator,

```
G(s) = G_0(s) + G_0 Σ G_0 + … = 1/(s − ω_k² − Σ(s)),
```

so the pole shifts: `Re Σ` is a mass renormalisation, `Im Σ` a width. A mode
stays a sharp, stable particle iff `Im Σ = 0` at its pole.

## The one-loop self-energy = the two-particle bubble

For a cubic self-interaction on the cavity (vertex `g_{knm} = ∫ ψ_k ψ_n ψ_m dr*`,
the triple overlap of the antipodal modes), the one-loop self-energy of mode `k`
is the bubble

```
Σ_k(s) = Σ_{n≤m} c_{nm} |g_{knm}|² / (s − (ω_n + ω_m)² + i0⁺),
```

(`c_{nm}` the symmetry factor) — the amplitude for `k → (n,m) → k` through the
two-particle intermediate state.

## Im Σ = 0 below threshold ⟹ the lightest mode is exactly stable

By the optical theorem `Im Σ_k(s)` is (minus) the two-particle phase space — it
is **nonzero only** when `s` reaches a two-particle threshold `(ω_n + ω_m)²`.
The lowest threshold is `2ω_0`. The lightest mode sits at `ω_0 < 2ω_0`, so its
pole `s = ω_0²` lies **below** `s_thr = (2ω_0)²`:

| quantity | value |
|---|---:|
| `ω_0` (lightest mode) | 1.167 |
| two-particle threshold `2ω_0` | 2.334 |
| pole `s = ω_0²` | 1.362 |
| threshold `s_thr = (2ω_0)²` | 5.448 |
| `Im Σ_0(ω_0²)` | ≈ 0 |

The lightest matter mode cannot decay (energy conservation) and stays a sharp,
real-pole, **stable** particle through one loop. (The two-particle density of
states is empty below `2ω_0`.)

## The real mass shift is finite

`Re Σ_0(ω_0²)` is a finite mass renormalisation: the cubic vertex overlaps
`g_{0nm}` decay with the mode index, so the mode sum converges —

| internal-mode cutoff | `Re Σ_0(ω_0²)` |
|---:|---:|
| 10 | −0.277 |
| 20 | −0.279 |
| 30 | −0.279 |
| 40 | −0.280 |

— stable to `~1e-3`, the residual UV piece being the same zeta/heat-kernel
regularisation as the #116 fluctuation determinant. A finite, real mass shift
(× coupling²) — no UV catastrophe on the discrete antipodal cavity.

## Unitarity survives one loop — and no horizon-absorption width

`Im Σ_k(s) ≤ 0` (a width) above the two-particle threshold and `= 0` below: the
dressed kernel respects the optical theorem (unitarity). Crucially, because the
throat is a **unitary mirror** (#129) there is **no horizon-absorption
contribution** to `Σ` — the only width source is genuine multi-particle decay,
which the lightest mode is kinematically forbidden from. This is the sharp
contrast with the absorbing horizon, which gives **every** mode a tree-level
width (#130). So the antipodal kernel's one-loop self-energy is
unitarity-preserving: it extends the tree-level stable spectrum (#130/#135) to
one loop.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | one-loop self-energy audit for the antipodal matter kernel (#135) |
| T2 | Dyson dressing | `G = 1/(s − ω_k² − Σ)`; Re Σ mass shift, Im Σ width |
| T3 | two-particle bubble | `Σ_k = Σ c\|g_{knm}\|²/(s−(ω_n+ω_m)²)`; vertex = mode triple overlap |
| T4 | lightest mode stable | `Im Σ_0(ω_0²) = 0` below threshold (`ω_0 < 2ω_0`) |
| T5 | finite mass shift | `Re Σ_0` mode sum converges (regularised, #116 scheme) |
| T6 | unitarity preserved | `Im Σ ≤ 0`, `= 0` below threshold; no horizon-absorption width (#129) |
| T7 | scope | one loop, fixed background, modelled vertex |
| T8 | assessment | `ANTIPODAL_KERNEL_ONE_LOOP_SELF_ENERGY_REAL_SHIFT_STABLE_LIGHTEST_MODE` |

## Established and open

  - **Established (BAM-native):** the one-loop self-energy of the antipodal
    matter kernel is a finite real mass shift with an imaginary (width) part
    that vanishes below the two-particle threshold — so the lightest matter mode
    stays exactly stable (`Im Σ = 0`), and unitarity survives one loop with no
    horizon-absorption width (the antipodal mirror, #129). One loop extends the
    tree-level stable spectrum (#130/#135).

  - **Does not / open:** the leading (one-loop) correction on the fixed
    background; the interaction **vertex** is **modelled** (a generic cubic
    triple-overlap), not derived from the S_BAM measure, and the coupling is an
    input — so `Re Σ` is fixed only up to the coupling. Higher loops, the
    absolute normalisation (#133), and the flavor residuals (#134) stand.

## Cross-references

  - `docs/antipodal_horizon_exchange_kernel_research_plan.md` — #135, the free
    propagator this self-energy dresses.
  - `docs/antipodal_vs_absorbing_qnm_research_plan.md` — #130, the tree-level
    stable spectrum (and the absorbing tree-level-width counterfactual).
  - `docs/null_throat_boundary_conditions_research_plan.md` — #129, the unitary
    mirror (no horizon-absorption width).
  - `docs/tangherlini_fluctuation_determinant_research_plan.md` — #116, the
    zeta/heat-kernel regularisation of the real mass shift.

## Run

```
python -m experiments.closure_ledger.antipodal_kernel_one_loop_self_energy_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_antipodal_kernel_one_loop_self_energy_probe/`.
Expected verdict:
`ANTIPODAL_KERNEL_ONE_LOOP_SELF_ENERGY_REAL_SHIFT_STABLE_LIGHTEST_MODE`, 8/8 PASS.

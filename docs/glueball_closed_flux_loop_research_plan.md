# Closed flux loops: pure-confinement benchmark + the BAM Möbius glueball tower (PR #100)

Glueballs — bound states of pure glue, i.e. closed flux loops with no
valence quarks — are the cleanest probe of the confinement geometry:
they carry NO quark masses and so are untouched by the flavor puzzle
(#97–#98). This probe benchmarks the BAM closed-flux-loop spectrum
against lattice QCD, and identifies where BAM's non-orientable topology
makes a **different** prediction than orientable-string QCD —
legitimately, because glueballs are not experimentally observed.

## The pure-confinement benchmark

A glueball is a closed flux loop of tension `σ` (the PR #99 string
tension, `√σ ≈ 0.42 GeV`); its mass is set entirely by `σ` — no quark
input — so glueball masses scale as `√σ`. The lattice values are
`0++ ≈ 4.1√σ` (1.73 GeV), `2++ ≈ 5.7√σ`, `0-+ ≈ 6.1√σ`. The BAM
closed-loop ground state (lowest closed-string level, `M² = 4πσ`) is

```
M = √(4πσ) ≈ 1.50 GeV ≈ 3.5 √σ,
```

the same `√σ` scale as the lattice 0++ (4.1√σ), agreeing to ~13% on the
`O(few)` coefficient. BAM and lattice agree because both are closed flux
loops of the SAME `σ` — a parameter-free benchmark (given `σ`) that BAM
passes.

## The closed-string Regge slope (half the meson)

A closed string has half the open-string Regge slope, so the glueball
trajectory slope is `α'_glueball = 1/(4πσ) = 0.44 GeV⁻²` (vs the meson
`1/(2πσ) = 0.88`), with an `M²` tower spacing of `2πσ = 1.13 GeV²`. The
observed pomeron (glueball) slope is ~0.25 GeV⁻², the same ballpark.

## Where BAM's topology diverges: the Möbius tower

The BAM machinery has TWO closed-loop sectors:

  - **orientable** — the glueball ring (`make_glueball_ring`, periodic,
    orientation `+1`): the standard closed string, matching lattice;
  - **non-orientable** — the Möbius tube (`make_mobius_tube`,
    antiperiodic, orientation `−1`): a half-twisted closed loop, where a
    single traversal reverses orientation.

The Möbius (antiperiodic) boundary condition shifts the mode
quantization from integer `n` to half-integer `n + ½`, so the
non-orientable glueball tower is shifted by `πσ` in `M²`:

> ⚠️ **The Möbius column below is superseded — see the correction after
> the table. The `½` and hence the `πσ` *level* shift stand; the absolute
> masses do not, because both towers were given a common intercept.**

| n | orientable (periodic) | Möbius (antiperiodic, **masses superseded**) |
|---:|---:|---:|
| 0 | 1.50 GeV | 1.68 GeV |
| 1 | 1.84 GeV | 1.99 GeV |
| 2 | 2.13 GeV | 2.26 GeV |
| 3 | 2.38 GeV | 2.49 GeV |

So BAM predicts a non-orientable Möbius glueball tower **interleaving**
the orientable one (~0.17 GeV above each), effectively **doubling** the
glueball spectrum. Orientable-string lattice QCD has no such sector.

> **CORRECTION (PR #229).** Two parts of the above have since been
> checked and one does not survive.
>
>   - **The `½` is confirmed.** This document *asserts* the
>     integer → half-integer shift; the probe's T4 returns `pass: True`
>     with no computation. PR #229 computes it three independent ways
>     (the #227 twisted ring, the repo's own `MobiusSpectrum`, a direct
>     antiperiodic solve), all giving exactly ½, and shows it is
>     topological: `O(1/M²)` convergence exactly, circumference-independent
>     to 1.4e-11. The level shift `+πσ` stands.
>   - **The common intercept does not.** The table above uses the *same*
>     `M₀²` for both towers. Antiperiodic moding also shifts the
>     zero-point (Casimir) energy — by `ΔE₀ = +π/(4C)` per transverse
>     polarization, computed two independent ways agreeing to 4.1e-07 —
>     and that scales as `1/C`, so it cannot be absorbed into a common
>     `M₀²`. **The Möbius ground state is therefore not `M₀² + πσ`, and
>     the absolute masses in the table above should not be used.** The
>     corrected intercept needs the full closed-loop dynamics this
>     document defers; #229 reports the gap rather than publishing a
>     replacement tower.
>
> The interleaving *structure* and the doubling of the spectrum are
> unaffected in kind, but their numerical placement is open.
> See `docs/mobius_identification_quantitative_research_plan.md`.

## Why this is legitimate (the non-observable point)

Glueballs are **not experimentally observed**: QCD predicts them, but
they mix heavily with ordinary qq̄ mesons of the same `J^PC` and have
never been cleanly identified. So the Möbius tower is a genuine
BAM-vs-lattice difference for a NON-observable — testable against lattice
QCD (which can isolate pure-glue states), but not contradicted by any
experiment. BAM's non-orientable topology is free to predict states that
orientable QCD does not, precisely where nature has not yet ruled.

## Tests

| # | test | finding |
|---|---|---|
| T1 | pure confinement | glueballs = closed loops, no quark input; `∝ √σ` |
| T2 | scale benchmark | orientable ground `√(4πσ) ≈ 1.50 GeV` (3.5√σ) ≈ lattice 0++ (4.1√σ), ~13% |
| T3 | Regge slope | closed-string glueball slope = half the meson |
| T4 | topology fork | BAM has orientable (periodic) + non-orientable (Möbius) loops |
| T5 | Möbius tower | half-integer modes shift `+πσ`, interleave orientable (≈2× states) |
| T6 | non-observable | glueballs unseen ⟹ Möbius tower legit vs lattice, not experiment |
| T7 | honest scope | scale benchmark passed; Möbius tower a topological extension |
| T8 | assessment | `GLUEBALL_SCALE_BENCHMARKS_LATTICE_MOBIUS_TOWER_IS_BAM_TOPOLOGY` |

## Established and open

  - **Established (BAM-native):** glueballs are the pure-confinement
    benchmark (no quark/flavor input); the BAM orientable ground state
    `√(4πσ) ≈ 1.50 GeV` benchmarks the lattice 0++ `√σ` scale (to ~13%);
    the closed-string glueball Regge slope is half the meson; and BAM's
    non-orientable Möbius sector predicts an extra glueball tower
    (half-integer modes, shifted by `πσ`) interleaving the orientable one
    — a topological prediction testable against lattice but not experiment
    (glueballs unobserved).

  - **Open:** the precise glueball `M/√σ` coefficients (full closed-loop
    dynamics: spin, mixing) — the robust statements are the `√σ` scale and
    the topological doubling; and a lattice cross-check with
    non-orientable boundary conditions (whether the Möbius tower is seen).

## Cross-references

  - `geometrodynamics/qcd/topology.py` — `make_glueball_ring` (orientable)
    and `make_mobius_tube` (non-orientable) closed loops.
  - `docs/qcd_confinement_cornell_audit_research_plan.md` — PR #99, the
    string tension `σ` and the flux-tube confinement audited here.

## Run

```
python -m experiments.closure_ledger.glueball_closed_flux_loop_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_glueball_closed_flux_loop_probe/`.
Expected verdict: `GLUEBALL_SCALE_BENCHMARKS_LATTICE_MOBIUS_TOWER_IS_BAM_TOPOLOGY`, 8/8 PASS.

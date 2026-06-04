# Non-perturbative convergence of the Z₂-graded BAM sector sum (PR #126)

PR #122 assembled the factorized sector sum

```
Z = Σ_{k odd, c₁, n_part} (−1)^k ∫ (dL/L) det^{−1/2}_matter · e^{i(π/2)(1−2a)} · e^{−S_BAM},
```

with the Z₂ grading `(−1)^k` (Möbius/odd-k orientation, PR #121) and the
continuous η-phase `e^{i(π/2)(1−2a)}` (U(1) holonomy, PR #121). It left **one
question open: does this sum converge non-perturbatively**, or is it a formal
expression that diverges term-by-term or in the sum? This PR audits that
convergence in its three independent pieces.

## The three pieces

The sector sum factorizes over three independent labels — the winding `k`, the
Hopf charge `c₁`, and the modulus `L` (loop length). Convergence of the whole
is convergence of each:

```
Z  =  ( Σ_{k odd}  winding )  ×  ( Σ_{c₁∈ℤ}  Hopf )  ×  ( ∫ dL/L  moduli ).
```

### 1. The winding sum is FINITE (a closure cap, not a tower)

The physical winding is **not** an infinite tower. The odd-k lemma (PR #121)
restricts `k` to odd values, and the available closure phase

```
Φ_avail(k) = 2π(k+1) + 50π·max(0, k−3)²
```

caps the admissible winding at `k ∈ {1, 3, 5}` — the three generations, with
the bulk dimension `k₅ = 5` the cap (PR #73). For `k = 7` the phase cost
(2563.5) is already far beyond the budget. So the winding sum is a **finite
sum of three terms**, trivially convergent — there is no winding tower to
diverge.

### 2. The Hopf-charge sum is a CONVERGENT theta

The Hopf charge `c₁ ∈ ℤ` carries an EM/Hopf action cost `~ c₁²`, so

```
Σ_{c₁∈ℤ} e^{−A c₁²} = √(π/A) · θ₃(0, e^{−π²/A}) → √(π/A),
```

a convergent Jacobi theta (verified at `A = 0.5, 1, 2` → `2.507, 1.773,
1.271`). Charge conservation `Σ c₁ = 0` (PR #58) constrains it further. The
Gaussian suppression makes the integer-charge sum absolutely convergent.

### 3. The moduli integral is FINITE at BOTH ends

The Z₂-graded moduli integral

```
I(m) = ∫₀^∞ (dt/t) [θ_per(t) − θ_anti(t)] e^{−m²t}
```

is finite at both ends, for two distinct reasons:

  - **UV (`t → 0`): the Z₂ grading cancels the divergence.** Individually
    `θ_per`, `θ_anti` ~ `1/√t` (a Weyl/heat-kernel divergence as `dt/t`). But
    the orientation difference cancels the boundary-condition-independent Weyl
    term, leaving `θ_per − θ_anti ~ e^{−π²/t} → 0` (PR #122). The integrand
    vanishes faster than any power at the UV — **no UV divergence**. (Numerically
    the integrand is `~9·10⁻¹⁴` already at `t = 0.02`.)

  - **IR (`t → ∞`): the mass gap kills the tail.** As `t → ∞`, `θ_per → 1`,
    `θ_anti → 0`, so the bracket → 1; the `e^{−m²t}` mass gap (the bounce
    saddle / the physical masses) makes the large-`t` tail integrable — **no IR
    divergence**.

The integral is finite for every mass gap: `I = 0.606` at `m = 0.3`, `0.173`
at `m = 0.5`, `0.0075` at `m = 1.0`.

## Overall

```
Z = (finite winding, 3 terms) × (convergent Hopf theta) × (finite moduli integral)
```

Each factor is finite, so the Z₂-graded sector sum **converges
non-perturbatively**. The grading is doing real work: it is precisely the
orientation signs `(−1)^k` that cancel the UV Weyl divergence in the moduli
integral — without the Z₂ grading the moduli integral would diverge at `t → 0`.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | audit non-perturbative convergence (PR #122 open item) |
| T2 | winding sum finite | `k ∈ {1,3,5}` closure cap; 3 terms, not a tower |
| T3 | Hopf sum convergent | `Σ e^{−A c₁²} = √(π/A)` convergent theta |
| T4 | moduli UV-convergent | Z₂ cancellation `θ_per − θ_anti ~ e^{−π²/t} → 0` |
| T5 | moduli IR-convergent | mass gap `e^{−m²t}` ⟹ finite tail |
| T6 | overall convergence | finite × convergent × convergent ⟹ converges |
| T7 | scope | convergence established; absolute scale open |
| T8 | assessment | `Z2_GRADED_SECTOR_SUM_CONVERGES_NONPERTURBATIVELY` |

## Established and open

  - **Established (BAM-native):** the factorized Z₂-graded sector sum of
    PR #122 converges non-perturbatively — the winding sum is a finite
    closure-capped sum (3 terms), the Hopf-charge sum is a convergent theta,
    and the moduli integral is finite at both ends (UV by Z₂ cancellation,
    IR by the mass gap). Finiteness is a genuine consequence of the grading.

  - **Does not / open:** this audits **convergence** (finiteness), not the
    **absolute normalization** — the bulk `κ₅²/Λ₅` anchor (PR #112) and the
    exact value of `Z` are not fixed here, nor is the multi-loop / interacting
    measure addressed.

## Cross-references

  - `docs/bam_factorized_sector_sum_research_plan.md` — PR #122, the factorized
    sum whose convergence is audited here.
  - `docs/bam_sector_phase_ledger_research_plan.md` — PR #121, the `(−1)^k`
    grading and continuous η-phase.
  - `docs/combined_matter_sector_aps_ledger_research_plan.md` — PR #125, the
    APS partition ledger over the same sum.

## Run

```
python -m experiments.closure_ledger.z2_graded_sector_sum_convergence_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_z2_graded_sector_sum_convergence_probe/`.
Expected verdict: `Z2_GRADED_SECTOR_SUM_CONVERGES_NONPERTURBATIVELY`, 8/8 PASS.

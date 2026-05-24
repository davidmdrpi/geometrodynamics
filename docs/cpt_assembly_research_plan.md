# Geometric CPT assembly for BAM throat histories

Assembles the three discrete symmetries — **C** (the inner/outer swap,
PR #63), **P** (parity), and **T** (`iσ_y`, B2) — into one geometric CPT
statement on throat histories (worldlines in the bulk `time × S³ ×
radial`). Each piece is already BAM-native; this probe shows their
product is the antiunitary CPT symmetry, guaranteed by the throat's
*local* Lorentz invariance (PRs #59–#60), and that it maps a throat
history to the **antithroat history run backwards** — the
Feynman–Stückelberg antiparticle, which is exactly the pair-production
structure (PR #58).

## The three operations (already BAM-native)

  - **C — charge conjugation = the inner/outer swap** (PR #63):
    `S: r ↦ 2R_MID − r`, an involution fixing the throat that flips the
    integrated Hopf curvature, `c₁ → −c₁` (throat → antithroat).
    `C² = +1`.

  - **P — parity = spatial reflection** of the `S³` angular slice:
    `x → −x`, `p → −p`; spin (axial) is `P`-even. `P² = +1`.

  - **T — time reversal = `iσ_y K`** (B2, antiunitary, `K` = complex
    conjugation): `t → −t`, `p → −p`, spin `s → −s`, energy `E → +E`
    (antiunitarity keeps `E > 0`). `T² = −I` — the fermionic signature
    (the non-trivial RP³ spin structure, B2).

## The transformation table

How the observables transform (sign = parity of the quantity):

| obs | C | P | T | **CPT** |
|---|---:|---:|---:|---:|
| charge `q` (`c₁`) | − | + | + | **−** |
| momentum `p` | + | − | − | **+** |
| position `x` | + | − | + | **−** |
| spin `s` | + | + | − | **−** |
| time `t` | + | + | − | **−** |
| energy `E` | + | + | + | **+** |

So **CPT**: `q → −q`, `p → +p`, `x → −x`, `s → −s`, `t → −t`, `E → +E`.
A particle `(q, p, s, E>0)` maps to an antiparticle `(−q, p, −s, E>0)`
with `x, t` reversed — the antiparticle running **backwards** in
spacetime.

## CPT on throat histories (Stückelberg = pair production)

A throat worldline going forward in time with charge `c₁ = +1` maps under
CPT to an antithroat (`c₁ = −1`) running backwards. This is the
Feynman–Stückelberg statement (antiparticle = particle reversed in time),
and it **is** the pair-production structure of PR #58: a throat–antithroat
pair is one worldline that turns around in time at the nucleation point
(a "V" in time), the two arms related by CPT.

## Why CPT holds (the theorem)

CPT is guaranteed for any local, Lorentz-invariant theory (Lüders–Pauli).
BAM's throat carries **local** Lorentz invariance (PRs #59–#60), so CPT
is exact at the local level. The closed `S³` breaks *global* Lorentz
invariance (a preferred frame, PR #59), so any CPT violation is suppressed
by the same `(R_MID/R_cosmo)² ~ 10⁻⁷⁸` — calculable and unobservably
small. An `O(1)` CPT violation would falsify; BAM passes.

## B4 accounting

C, P, T and CPT are **discrete geometric operations**; the charge `c₁` is
a dimensionless topological integer, the spin signature `T² = −1` a
group fact, the sign table pure parities. Nothing carries a scale —
independent of the single anchor `m_e` (B4-consistent).

## Tests

  T1. **The three operations.** C = inner/outer swap (#63); P = `S³`
      reflection; T = `iσ_y K` (B2). Defined, each acting on the throat.
  T2. **(Anti)involution signatures.** `C² = +1`, `P² = +1`, `T² = −I`
      (fermionic; verify `(iσ_y K)² = −I`).
  T3. **CPT transformation table.** The sign product `C·P·T` gives
      `q→−, p→+, x→−, s→−, t→−, E→+`.
  T4. **BAM realizations.** C ↔ the inner/outer swap (`c₁ → −c₁`, #63);
      T ↔ `iσ_y` (`T² = −I`, B2).
  T5. **Stückelberg / pair production.** CPT(throat forward) =
      antithroat backward = one worldline turning in time (PR #58).
  T6. **CPT theorem from local Lorentz invariance.** Guaranteed by local
      Lorentz (PRs #59–#60); global breaking on `S³` → suppressed
      violation `(R_MID/R_cosmo)²`.
  T7. **Falsification / B4.** `O(1)` CPT violation would falsify; BAM
      passes (suppressed). The operations are dimensionless/geometric.
  T8. **Assessment.**

## Verdict structure

  - **CPT_ASSEMBLED** (expected): C (inner/outer swap, #63), P (`S³`
    reflection), and T (`iσ_y`, B2) compose to the antiunitary CPT
    symmetry — `q→−, p→+, x→−, s→−, t→−, E→+`, with `C²=P²=+1`,
    `T²=−I` — mapping a throat history to the antithroat history run
    backwards (Stückelberg = pair production, #58). CPT is guaranteed by
    the throat's local Lorentz invariance (#59–#60); global `S³` breaking
    gives a suppressed `(R_MID/R_cosmo)²` violation. The discrete-symmetry
    sector is unified.

  - **CPT_FAILS**: the operations do not compose to CPT, the signatures
    are wrong, or CPT is not realized on throat histories.

## What this leaves open

  - **The full CPT operator on the throat Dirac spinor** from `S_BAM`
    (the explicit `Θ = CPT` matrix and `Θ²`), beyond the sign table.
  - **P and the antipodal `Z₂`.** Disentangling spatial parity from the
    antipodal deck transformation of `RP³ = S³/Z₂` (B2) precisely.
  - **Observable CPT bounds.** Mapping `(R_MID/R_cosmo)²` to specific
    CPT-violation observables.

## Cross-references

  - `docs/charge_conjugation_swap_research_plan.md` — C = inner/outer
    swap (#63).
  - `docs/topological_discrete_sector_research_plan.md` — `T = iσ_y`,
    the antipodal `Z₂` / RP³ spin structure (B2).
  - `docs/pair_production_threshold_research_plan.md` — the throat–
    antithroat pair (#58).
  - `docs/stable_moving_throat_research_plan.md` — local Lorentz
    invariance (#59).
  - `geometrodynamics/embedding/transport.py` — `T = iσ_y`.
  - `experiments/closure_ledger/cpt_assembly_probe.py` — this probe.

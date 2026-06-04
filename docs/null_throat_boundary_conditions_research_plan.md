# Null throat boundary conditions for wave transport on the 5D horizon (PR #129)

PR #128 built the horizon-regular charts for the 5D Tangherlini throat and
identified the antipodal map `(U,V,Ω) → (−U,−V,Ω_antipodal)` as BAM's
throat ↔ antithroat C-swap. PR #116 ran the matter cavity operator
`H = −d²/dr*² + V_tangherlini` on the exterior. This PR answers the next
question: **what boundary condition does the null throat (the 5D horizon)
impose on transported waves?** The answer is the BAM-native one — an
**l-parity-graded antipodal condition** that makes the throat a **unitary
mirror**, not an absorbing horizon.

All facts are checked numerically (vanishing potential, S³ harmonic parity,
flux, the cavity spectrum).

## The potential vanishes at the horizon

The massless wave equation `□Φ = 0`, separated as
`Φ = e^{−iωt} Y_l(Ω) ψ(r)/r^{3/2}`, gives the Schrödinger form
`−d²ψ/dr*² + V_l ψ = ω²ψ` with `V_l = f(r)[l(l+2)/r² + 3rs²/r⁴]` (PR #116).
Since `V_l ∝ f(r) → 0` at the throat `r = rs` (verified for `l = 0..3`), the
near-horizon equation is `−ψ'' = ω²ψ`: the modes are the pure null phases
`ψ ~ e^{±iωr*}` — the ingoing/outgoing null rays of PR #128
(`v = t + r*`, `u = t − r*`).

## Three candidate boundary conditions

| BC | form at the throat | character |
|---|---|---|
| ingoing / absorbing (standard QNM) | `ψ ~ e^{−iωr*}` | flux into the horizon; **non-unitary** (sink) |
| reflective wall | Dirichlet `ψ=0` or Neumann `ψ'=0` | real discrete spectrum (the box of PR #116) |
| **antipodal (BAM)** | `Φ(U,V,Ω) = Φ(−U,−V,Ω_antipodal)` | the PR #128 identification |

## The antipodal map fixes the BC by l-parity

The scalar harmonics on the horizon `S³` carry antipodal parity
`Y_l(−x) = (−1)^l Y_l(x)` — they are degree-`l` harmonic polynomials on `ℝ⁴`
(verified for `l = 1,2,3,4`). Single-valuedness of the untwisted scalar under
the antipodal identification, `Φ(U,V,Ω) = Φ(−U,−V,Ω_antipodal)`, then forces
the radial function to carry the compensating parity `(−1)^l` across the
throat:

| l-parity | angular factor `(−1)^l` | radial parity | throat BC |
|---|---:|---|---|
| **even** | `+1` | even | **Neumann** `ψ'(throat) = 0` (antinode) |
| **odd** | `−1` | odd | **Dirichlet** `ψ(throat) = 0` (node) |

So the null throat BC is not a single choice — it is the **l-parity-graded
antipodal condition**. A twisted/Möbius field flips even ↔ odd, connecting to
the Z₂ orientation grading (PRs #67/#121).

## The throat is a unitary mirror

Both even-l Neumann and odd-l Dirichlet are **real** boundary conditions, so
the Klein–Gordon/Wronskian flux `j ∝ Im(ψ*ψ')` through the throat **vanishes**:
the throat is a perfect, unitary mirror — no net flux is lost into it. This is
the sharp contrast with the ingoing/absorbing horizon BC `ψ ~ e^{−iωr*}`, whose
flux `j = −ω|A|² ≠ 0` carries probability into the hole (verified: the ingoing
mode gives `j = −ω = −1.3`, the real standing wave gives `j = 0`). The
antipodal throat conserves flux (unitarity) — consistent with the global
CPT / unitarity of the throat histories (PR #64): what falls toward the throat
on one sheet re-emerges on the antipodal sheet, nothing is destroyed.

## The cavity spectrum is real, discrete, and l-parity-graded

On the exterior `[R_MID+ε, R_OUTER]` with the outer shell wall (Dirichlet at
`R_OUTER`) and the antipodal BC at the throat, the operator `−d²/dr*² + V_l`
has a real, positive, discrete spectrum (a unitary bound cavity), with even-l
(Neumann-at-throat) and odd-l (Dirichlet-at-throat) families giving **distinct**
spectra (e.g. lowest `ω²`: `l=0` Neumann `1.37`, `l=1` Dirichlet `5.27`,
`l=2` Neumann `2.03`, `l=3` Dirichlet `6.73`). This is the wave-transport face
of BAM's even-k/odd-k structure (#67 even-k absence, #121 odd-k lemma).

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | derive the null throat BC for wave transport on the 5D horizon |
| T2 | vanishing potential | `V_l ∝ f → 0` at the throat ⟹ null modes `e^{±iωr*}` |
| T3 | candidate BCs | ingoing/absorbing, reflective wall, antipodal |
| T4 | l-parity BC | `Y_l(−x) = (−1)^l Y_l` ⟹ even-l Neumann, odd-l Dirichlet |
| T5 | unitary mirror | real BC ⟹ zero flux; ingoing ⟹ flux `−ω` (absorbing) |
| T6 | real discrete spectrum | even-l (N) vs odd-l (D) distinct; even/odd Z₂ face |
| T7 | scope | kinematic BC established; QNM / nucleation open |
| T8 | assessment | `NULL_THROAT_BC_ANTIPODAL_L_PARITY_UNITARY_MIRROR` |

## Established and open

  - **Established (BAM-native):** the null throat imposes the antipodal,
    l-parity-graded boundary condition — even-l Neumann, odd-l Dirichlet, from
    the S³ harmonic parity `Y_l(−x) = (−1)^l Y_l(x)` under the PR #128 antipodal
    identification — making the throat a **unitary, flux-conserving mirror**
    rather than an absorbing horizon; the exterior cavity spectrum is real,
    discrete, and even/odd-graded.

  - **Does not / open:** this is the *kinematic* BC structure of classical wave
    transport. It does **not** solve the full quasinormal-mode spectrum
    (complex `ω`, ringdown), nor compute the dynamical throat ↔ antithroat
    *nucleation* amplitude (the bounce rate, PR #58/#88). The
    absorbing-vs-antipodal distinction is fixed by the BAM antipodal postulate
    (PR #128); this probe shows that postulate yields a self-consistent,
    unitary, l-graded transport condition.

## Cross-references

  - `docs/five_d_tangherlini_throat_horizon_lift_research_plan.md` — PR #128,
    the antipodal Kruskal structure realised as the BC here.
  - `docs/tangherlini_fluctuation_determinant_research_plan.md` — PR #116, the
    cavity operator whose throat BC is derived here.
  - `docs/cpt_assembly_research_plan.md` — PR #64, the global CPT / unitarity
    consistent with the unitary-mirror throat.

## Run

```
python -m experiments.closure_ledger.null_throat_boundary_conditions_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_null_throat_boundary_conditions_probe/`.
Expected verdict: `NULL_THROAT_BC_ANTIPODAL_L_PARITY_UNITARY_MIRROR`, 8/8 PASS.

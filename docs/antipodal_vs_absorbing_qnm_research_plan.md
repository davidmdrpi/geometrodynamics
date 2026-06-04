# Antipodal vs absorbing throat quasinormal spectrum (PR #130)

PR #129 derived the boundary condition the null throat imposes on the matter
waves: the BAM-native **antipodal** condition (real, l-parity-graded — even-l
Neumann, odd-l Dirichlet) makes the throat a **unitary mirror**, whereas a
standard black-hole horizon would impose the **absorbing** (ingoing)
condition. PR #129 showed the antipodal cavity has a real, discrete spectrum.
This PR computes the **full complex frequency spectrum** for both boundary
conditions and contrasts them — the spectral fingerprint distinguishing BAM's
antipodal throat from an ordinary absorbing horizon.

## The eigenproblem

The separated radial wave equation `−d²ψ/dr*² + V_l ψ = ω²ψ` runs on the BAM
cavity `[R_MID+ε, R_OUTER]`, with the shell wall (Dirichlet) at `R_OUTER` and
one of two throat conditions:

| throat BC | form | operator |
|---|---|---|
| **antipodal** (BAM, #129) | Neumann `ψ'=0` (even l), Dirichlet `ψ=0` (odd l) | **self-adjoint** |
| **absorbing** (ordinary horizon) | ingoing `ψ'(throat) = −iω ψ(throat)` | **non-self-adjoint** |

The absorbing case is a **quadratic eigenvalue problem** in `ω` (ω in the BC,
ω² in the bulk), solved by the companion linearisation
`(K₀ + ωK₁ + ω²K₂)ψ = 0 → A z = ω B z`.

## The spectral contrast

| l | antipodal ω (BC) | absorbing ω (ingoing) |
|---:|---|---|
| 0 | 1.181, 3.293, 5.439 (N) | 1.893−1.241i, 3.917−1.166i, 6.009−1.044i |
| 1 | 2.295, 4.363, 6.478 (D) | 1.995−1.126i, 3.973−1.140i, 6.044−1.035i |
| 2 | 1.436, 3.443, 5.532 (N) | 2.168−1.000i, 4.063−1.100i, 6.102−1.021i |
| 3 | 2.594, 4.536, 6.595 (D) | 2.389−0.896i, 4.185−1.050i, 6.182−1.002i |

  - **Antipodal ⟹ real ω (undamped normal modes).** The self-adjoint operator
    has a real spectrum (`max|Im ω| = 0` to numerical precision). With the
    `e^{−iωt}` convention these are infinitely-lived normal modes — sharp
    spectral lines, quality factor `Q = ∞`, zero width — l-parity graded.
  - **Absorbing ⟹ complex ω (damped ringdown).** The ingoing BC makes the
    operator non-self-adjoint; the eigenfrequencies are `ω = ω_R − i|ω_I|`
    with `Im(ω) < 0` — damped quasinormal modes, lifetime `τ = 1/|ω_I|`.

## Quality factor: sharp normal modes vs leaky ringdown

| l | Q (antipodal) | Q (absorbing) | τ (absorbing) |
|---:|---|---|---|
| 0 | ∞ | 0.76 | 0.81 |
| 1 | ∞ | 0.89 | 0.89 |
| 2 | ∞ | 1.08 | 1.0 |

`Q = ω_R/(2|ω_I|)` is infinite for the antipodal throat (a sharp δ-line) but
`O(1)` for the absorbing horizon (a Lorentzian of width `Γ = 2|ω_I|`, finite
lifetime `τ = 1/|ω_I|`) — the thin cavity leaks fast into the throat.

## The physical consequence: stable matter needs the unitary throat

A matter state is a sharp mass (a stable or long-lived particle) only if its
cavity mode has a **real** frequency. The absorbing throat gives every mode a
width / complex mass (a decaying resonance); only the antipodal, unitary throat
(PR #129) yields the real, stable spectrum the BAM matter sectors (the
lepton/quark bound states) require. The undamped-vs-ringdown distinction is the
spectral face of the program's global CPT / unitarity (PR #64): **BAM matter is
stable precisely because the throat reflects antipodally rather than
absorbing.**

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | compute & contrast antipodal vs absorbing throat QNM spectrum |
| T2 | setup | cavity operator, shell wall, antipodal vs ingoing throat BC (QEP) |
| T3 | antipodal real | `Im(ω) ≈ 0`, undamped normal modes (`Q = ∞`), l-graded |
| T4 | absorbing complex | `Im(ω) < 0`, damped quasinormal ringdown |
| T5 | quality factor | antipodal `Q = ∞`, absorbing `Q ~ O(1)`, `τ = 1/|ω_I|` |
| T6 | stable matter | sharp real masses need the unitary antipodal throat (#64/#129) |
| T7 | scope | finite-cavity spectrum; idealised horizon QNM / GW open |
| T8 | assessment | `ANTIPODAL_THROAT_REAL_UNDAMPED_VS_ABSORBING_COMPLEX_RINGDOWN` |

## Established and open

  - **Established (BAM-native):** the antipodal throat (PR #129) gives a real,
    undamped normal-mode spectrum (`Im ω = 0`, `Q = ∞` — sharp stable lines),
    while an absorbing horizon gives complex quasinormal frequencies
    (`Im ω < 0`, damped ringdown, `Q ~ O(1)`). Stable matter — the BAM
    lepton/quark bound states — requires the unitary antipodal throat.

  - **Does not / open:** this is the spectrum of the *finite* BAM cavity (the
    physically appropriate region). It does **not** compute the idealised
    `r* → −∞` horizon QNMs (the full ringdown tower), the coupling to
    gravitational radiation, or the absolute mode normalisation. The absorbing
    case is the counterfactual — BAM selects the antipodal BC (PR #129); the
    contrast shows what is at stake. The `Im(ω) < 0` sign is the `e^{−iωt}`
    decay convention.

## Cross-references

  - `docs/null_throat_boundary_conditions_research_plan.md` — PR #129, the
    antipodal l-parity BC whose spectrum is computed here.
  - `docs/tangherlini_fluctuation_determinant_research_plan.md` — PR #116, the
    cavity operator.
  - `docs/cpt_assembly_research_plan.md` — PR #64, the global CPT / unitarity
    realised here as the undamped (stable-matter) spectrum.

## Run

```
python -m experiments.closure_ledger.antipodal_vs_absorbing_qnm_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_antipodal_vs_absorbing_qnm_probe/`.
Expected verdict:
`ANTIPODAL_THROAT_REAL_UNDAMPED_VS_ABSORBING_COMPLEX_RINGDOWN`, 8/8 PASS.
(Runs in ~40s: eight QEP solves.)

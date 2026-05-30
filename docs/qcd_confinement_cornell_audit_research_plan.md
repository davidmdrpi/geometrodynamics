# QCD confinement geometry: Cornell / flux-tube string-tension audit (PR #99)

The quark MASS sector terminated honestly at the flavor puzzle (#97–#98)
— the Yukawa magnitudes are not geometric. This probe pivots to the QCD
**confinement** sector, the part of QCD that IS geometric in BAM (flux
tubes = wormhole bridges), and audits the Cornell potential and the
flux-tube string tension: which content is BAM-native geometry, and which
is the single QCD scale anchor.

## The Cornell potential

The BAM QCD machinery (`geometrodynamics/qcd/bridge.py`) uses the Cornell
static energy

```
V(L) = σ·L − A·ℏc/L,    σ = 0.18 GeV²,  A = 0.30.
```

  - **Linear `σ·L` — the flux tube = a wormhole bridge of constant
    tension.** Confinement is a 1D throat-network bridge connecting the
    quark–antiquark; its energy per unit length is constant — the
    defining property of a confining string.
  - **Coulomb `−A·ℏc/L` — short-distance throat/gluon exchange.** The
    short-range piece is one-gluon exchange, the QCD analogue of the
    lepton Coulomb law BAM derived from eigenmode throat flux.

## String breaking = Schwinger pair nucleation = PR #58 (eE → σ)

The flux tube does not stretch forever: at large `L` it breaks by
producing a quark–antiquark pair. The BAM bridge nucleates with the
Schwinger rate

```
Γ ∝ exp(−π m_q² / (σ L)),
```

the QED Schwinger formula `exp(−π m_e²/(eE))` with the electric field
replaced by the string tension, **eE → σ**. This is precisely the PR #58
throat-pair-production mechanism (`e E_S · R_MID = m_e c²`) transported to
QCD: the confining string is a tense brane, and when its work `σ·L`
reaches the pair threshold `≈ 2 m_q` the throat-pair nucleates and the
string snaps. QCD string breaking and lepton pair production are the
**same BAM nucleation physics with `eE ↔ σ`**.

## Consistency checks the BAM σ passes

| check | BAM | observed |
|---|---|---|
| confinement scale | `√σ = 424 MeV` | lattice ~440 MeV |
| Regge slope `α' = 1/(2πσ)` | `0.884 GeV⁻²` | ~0.88–0.93 GeV⁻² |
| string-breaking length (`σ·L = 2 m_q`) | `L ≈ 1.4–1.6 fm` | lattice `L_break ≈ 1.35 fm` |

## The one QCD anchor (B4 analogue)

Just as the lepton sector rides on the single dimensionful anchor
`m_e = ℏc/R_MID` (B4), the confinement sector rides on a single scale —
`√σ ≈ 0.42 GeV`, the Λ_QCD scale. The Cornell FORM (linear + Coulomb),
the string-breaking = Schwinger mechanism, and the Regge slope
`α' = 1/(2πσ)` are all geometric / dimensionless-derived; the absolute
value of `σ` is calibrated to the QCD scale, not derived — exactly the B4
pattern (form geometric, one scale anchored).

## Tests

| # | test | finding |
|---|---|---|
| T1 | Cornell form | `V(L)=σL − A·ℏc/L` (Coulomb small-L, linear large-L) |
| T2 | flux tube | linear `σL` = wormhole bridge of constant tension |
| T3 | Coulomb | `−A·ℏc/L` = short-distance throat/gluon exchange |
| T4 | string breaking | Schwinger `exp(−πm_q²/(σL))` = PR #58 (`eE→σ`); `σL≈2m_q` |
| T5 | Regge slope | `α'=1/(2πσ)=0.884 GeV⁻²` vs observed ~0.88–0.93 |
| T6 | QCD anchor | `√σ ≈ 424 MeV` = the one QCD scale (B4 analogue) |
| T7 | honest scope | structure BAM-native; scale `√σ` anchored |
| T8 | assessment | `CONFINEMENT_GEOMETRIC_STRING_BREAK_IS_SCHWINGER_SCALE_IS_QCD_ANCHOR` |

## Established and open

  - **Established (BAM-native):** the Cornell linear term is the
    flux-tube wormhole-bridge of constant tension; string breaking is the
    PR #58 Schwinger throat-pair nucleation with `eE → σ`; the BAM `σ`
    reproduces the Regge slope `α' = 1/(2πσ) ≈ 0.88 GeV⁻²` and the
    string-breaking length; `√σ ≈ 0.42 GeV` is the single QCD scale anchor
    (B4 analogue).

  - **Open:** a first-principles value of `σ` — it is the Λ_QCD scale
    anchor, calibrated to lattice, like `m_e` for leptons.

## Cross-references

  - `geometrodynamics/qcd/bridge.py` — the Cornell static energy and
    flux-tube bridge.
  - `geometrodynamics/qcd/diagnostics.py` — the lattice σ extraction from
    the Schwinger bridge-nucleation rate.
  - `docs/pair_production_threshold_research_plan.md` — PR #58, the
    throat-pair nucleation / Schwinger threshold (`e E_S R_MID = m_e c²`),
    of which QCD string breaking is the `eE→σ` analogue.

## Run

```
python -m experiments.closure_ledger.qcd_confinement_cornell_audit_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_qcd_confinement_cornell_audit_probe/`.
Expected verdict: `CONFINEMENT_GEOMETRIC_STRING_BREAK_IS_SCHWINGER_SCALE_IS_QCD_ANCHOR`, 8/8 PASS.

# n_part = 233 revisited: the quark hierarchy is the program's one dynamical sector (PR #97)

Returns to the hardest open piece of the quark sector — the origin of
`n_part = 233` (the quark closure integer `N_q = 2·n_part = 466`,
`β_quark = 233π`). PR #76 classified it as a phenomenological compensator
(only parity §8-stable; `n_part` drifts 216–255) and named the right
route as "a quantitative QCD-shell model." This probe revisits it with
the machinery built since — and sharpens the diagnosis using the
now-complete lepton and neutrino sectors.

## The fresh angle: a huge hierarchy can be geometric

When PR #76 ran, the implicit worry was that the ~9-order quark mass²
hierarchy was simply too large for the geometric closure machinery. The
neutrino arc (#86–#90) overturned that intuition: it derived a hierarchy
of comparable size — the keV→TeV seesaw `M_R = m_D·e^{S}`, ~10⁶ in mass —
as a **clean geometric exponential** (the non-orientable tortoise bounce,
an O(15) action). So a huge hierarchy does NOT force a non-geometric
origin. The BAM mass program now has two geometric hierarchy *types*:

  - **charged leptons**: a closure-ledger ladder with the clean,
    §8-stable closure integer `4·k_5² = 100` (PR #71) — power-law-like;
  - **neutrinos**: `m_ν = m_D·e^{−S}`, the tortoise-bounce exponential
    (PR #88–90).

The question is then sharper: not "is the quark hierarchy too big to be
geometric" but "does it have the **regularity** of a geometric hierarchy
(power-law or exponential)?"

## The quark hierarchy is irregular

It does not. The consecutive up-type mass ratios are `m_c/m_u ≈ 588` and
`m_t/m_c ≈ 136` — not constant, so NOT a clean exponential (geometric
progression); and not the `k²`-style pattern of a power law. The
down-type ratios (`m_s/m_d ≈ 20`, `m_b/m_s ≈ 45`) differ from the
up-type, so the two Z₂ partitions are asymmetric. The quark hierarchy is
irregular — matching neither geometric type.

## The geometric shell cannot carry it

The quark shell basis (PR #77/#83: `k=0`, overtones `n=3,4,5`) has
`ω²(l=1, n) = 14.6, 22.7, 32.5` — a mass² span of only **×2.2** — whereas
the observed quark mass² span (`t/u`) is **×6.4×10⁹**. The geometric
shell under-produces by ~2.9×10⁹; `n_part` is the empirical price.

## The diagnosis: a dynamical (QCD-RG) hierarchy

Irregularity across a wide scale range is the signature of
renormalisation-group running: `α_s` runs logarithmically, so the quark
masses are QCD-dressed differently at each scale — an intrinsically
**dynamical** hierarchy. Leptons and neutrinos do not feel QCD, so their
hierarchies are geometric; the quarks do, so theirs is dynamical. This is
why the quark closure integer is the ONLY one that drifts under §8
(216–255): it is absorbing dynamical content that no geometric closure
quantity encodes. The lepton↔quark closure gap
`N_q − N_lepton = 466 − 100 = 366` quanta is precisely that dynamical
(QCD) excess.

| sector | hierarchy type | closure quantity | §8 |
|---|---|---|---|
| charged leptons | power-law (ledger) | `4·k_5² = 100` | **stable** |
| neutrinos | exponential `e^{−S}` (bounce) | `S ≈ 16` (geometric) | stable form |
| quarks | **irregular (QCD-RG)** | `N_q = 466` (`n_part=233`) | **drifts 216–255** |

So `n_part` is not merely fit on the "wrong basis" (PR #76) — the quark
hierarchy is the BAM mass program's **one dynamical sector**, and a
geometric closure integer can only compensate it, never derive it. The
right route is sharpened: a QCD-shell model **with `α_s` running** (the
missing ingredient is the RG dynamics, identified now by contrast with
the geometric lepton and neutrino sectors).

## Tests

| # | test | finding |
|---|---|---|
| T1 | recap | `n_part=233` = v3 compensator (parity-only §8 inv.; drift 216–255) |
| T2 | geometric sectors | leptons (integer 100) + neutrinos (`e^{−S}`) are geometric |
| T3 | huge ≠ non-geometric | the neutrino exponential is huge and geometric |
| T4 | irregular | quark `c/u≈588` vs `t/c≈136` (not exp), up/down asym (not power) |
| T5 | shell under-carries | shell span ×2.2 vs observed ×6.4×10⁹ |
| T6 | dynamical diagnosis | irregular ⟹ QCD-RG; gap `366` = dynamical excess |
| T7 | honest scope | #76 upheld + sharpened; route = QCD-shell + `α_s` running |
| T8 | assessment | `QUARK_HIERARCHY_DYNAMICAL_N_PART_COMPENSATES` |

## Established and open

  - **Established:** the PR #76 compensator verdict is upheld and
    sharpened. Size is not the obstruction (the neutrino exponential is
    geometric and huge); the obstruction is that the quark hierarchy is
    **irregular** (neither power-law nor exponential), the signature of
    QCD-RG dynamics. The geometric shell carries only ×2.2 of the
    ×6.4×10⁹ span; `n_part` (and the 366-quantum lepton↔quark gap) is the
    dynamical excess. The quark sector is the program's one dynamical
    (non-geometric) hierarchy.

  - **Open:** a first-principles `n_part = 233` — none exists in the
    geometric machinery, and by the diagnosis none should. The derivation
    route (a QCD-shell model with `α_s` running) is a substantial program
    outside the closure-ledger machinery.

## Cross-references

  - `docs/quark_npart_origin_research_plan.md` — PR #76, the prior
    classification (compensator; wrong-basis reading).
  - `docs/quark_beta_status.md` — the five-probe sequence and the parity
    invariance; the "366-quanta lepton/quark gap".
  - `docs/beta_lepton_derivation_research_plan.md` — PR #71, the lepton
    closure integer `4·k_5² = 100`.
  - `docs/majorana_bounce_action_research_plan.md`,
    `docs/boundary_compliance_bulk_geometry_research_plan.md` — PR #88–90,
    the neutrino exponential hierarchy `m_ν = m_D·e^{−S}`.
  - `docs/throat_shell_mass_operator_unification_research_plan.md` —
    PR #83, the unified operator (lepton throat `k≠0` vs quark shell `k=0`).

## Run

```
python -m experiments.closure_ledger.npart_dynamical_hierarchy_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_npart_dynamical_hierarchy_probe/`.
Expected verdict: `QUARK_HIERARCHY_DYNAMICAL_N_PART_COMPENSATES`, 8/8 PASS.

# Gauge–matter coupling from the antipodal throat boundary (PR #141)

> **Framing (to avoid a category error).** This is quantum field theory on the
> *fixed classical* throat geometry — geometry → fields, **not** quantum
> gravity. "Gauge" is the U(1)_Hopf field (the photon of PRs #42–#44); "matter"
> is the antipodal cavity modes (#116, #129–#140); the coupling is read off the
> *classical* antipodal throat boundary. The metric is never quantised.

The matter arc (#129–#140) built the BAM matter sector on the antipodal throat;
the gauge arc (#42–#44) built the photon exchange kernel (`1/q²` from the S³
Green function). This PR **joins them**: how the U(1)_Hopf gauge field couples to
the antipodal matter modes **at the throat**, consistently with the antipodal
identification. The answer ties the two sectors through the C-swap (#63): the
throat is the surface where charge conjugation acts, the gauge–matter vertex
inherits the **same antipodal Z₂ selection rule** as the cubic vertex (#137/#140),
and U(1) charge is conserved by the unitary mirror (#129). Only the coupling
**strength** — `α`, the EM coupling (#105) — stays the universal input residual.

## Minimal coupling

Matter charged under U(1)_Hopf with charge `c₁` (`|c₁| = 1`, #58/#74) couples
minimally through the gauge-covariant derivative

```
D_μ = ∂_μ − i c₁ A_μ,
```

giving the gauge–matter vertex `c₁ ∫ A_μ j^μ`, with `j^μ = ψ*∂^μψ − (∂^μψ*)ψ` the
matter current. This is the standard gauge-invariant coupling of the #42–#44
photon to the matter current.

## The C-swap is spatial inversion × charge conjugation

The antipodal map `A : x → −x` (the throat ↔ antithroat C-swap, #63) acts on
both sectors at once:

  - on the matter harmonics, `Y_l → (−1)^l Y_l` (the parity that fixed the BC
    #129 and graded the vertices #137/#140);
  - on the Hopf charge, `c₁ → −c₁` (charge conjugation, #63).

So the C-swap is **one operation with two effects** — a spatial inversion and a
charge conjugation — and the throat is the **particle ↔ antiparticle (C)
surface** (#63/#64). This is exactly why the gauge field can couple to the matter
there.

## The gauge–matter vertex inherits the antipodal Z₂ selection rule

The vertex couples a photon (angular content `l_γ`) to two matter legs
(`l₁, l₂`): its angular part is the triple overlap `∫_{S³} Y_{l_γ} Y_{l₁} Y_{l₂}
dΩ` — the **same structure as the cubic matter vertex** (#137). S_BAM's antipodal
invariance (the #140 Ward identity) forces `Σl = l_γ + l₁ + l₂` even:

| vertex (photon · matter·matter) | Σl | even? | ∫YYY | allowed? |
|---|---:|:---:|---:|:---:|
| γ(1) · matter(1,0) | 2 | ✓ | 0.25 | ✓ |
| γ(1) · matter(1,2) | 4 | ✓ | 0.0417 | ✓ |
| γ(0) · matter(1,1) | 2 | ✓ | 0.25 | ✓ |
| γ(1) · matter(1,1) | 3 | ✗ | 0 | ✗ |

The gauge–matter vertex obeys the same antipodal Z₂ selection rule, now with the
gauge leg (verified exactly via the S³ monomial integral — odd-Σl forbidden).

## U(1) charge is conserved at the throat

The antipodal throat is a unitary mirror for the matter (#129): zero net matter
flux through it. The same boundary conserves the gauge charge flux. Combined with
the C-swap charge flip (`c₁ → −c₁`), outgoing charge re-emerges as the conjugate
on the antipodal sheet, so charge is conserved — `Σc₁ = 0` (#58) — and the throat
balances particle against antiparticle. Charge conservation at the throat is the
**gauge face of the unitary mirror**.

## The coupling strength is α (input)

The minimal-coupling **structure** — the covariant derivative, the triple-overlap
vertex with the Σl-even selection rule, the charge conservation, the throat as
the C-surface — is derived from the antipodal geometry. The coupling **strength**
is the EM coupling `α` (the "137 problem", #105), a universal residual not fixed
by the geometry. So the structure is BAM-native; the magnitude `α` is input.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | gauge–matter coupling from the antipodal throat boundary |
| T2 | minimal coupling | `D_μ = ∂_μ − i c₁ A_μ`; vertex `c₁ ∫ A_μ j^μ` |
| T3 | C-swap | inversion (`Y_l→(−1)^l`) × charge conjugation (`c₁→−c₁`); C-surface |
| T4 | gauge vertex | Σl-even triple overlap (antipodal Z₂, #137/#140) |
| T5 | charge conservation | unitary mirror (#129) + C-swap flip ⟹ `Σc₁ = 0` (#58) |
| T6 | coupling strength | `α` (#105) — input residual |
| T7 | ledger / scope | structure derived; `α` / normalisation open |
| T8 | assessment | `GAUGE_MATTER_COUPLING_MINIMAL_Z2_SELECTED_CHARGE_CONSERVED_ALPHA_INPUT` |

## Established and open

  - **Established (BAM-native):** the U(1)_Hopf gauge field couples minimally to
    the antipodal matter modes at the throat; the C-swap is inversion × charge
    conjugation (the throat is the C-surface, #63/#64), the gauge–matter vertex
    inherits the antipodal Z₂ selection rule (`Σl` even, #137/#140), and U(1)
    charge is conserved by the unitary mirror (`Σc₁ = 0`, #58/#129). The
    coupling **structure** is derived.

  - **Does not / open:** the coupling **strength** `α` (the 137 problem,
    #105/#108) is the universal input residual, not fixed by the geometry; the EM
    normalisation, higher gauge vertices, and the running of `α` are not
    addressed. The `α` (#105/#108), bulk-scale (#133), and flavor (#134)
    residuals stand.

## Cross-references

  - `bam_exchange_kernel_probe.py` (PRs #42–#44) — the gauge-sector photon
    exchange kernel (`1/q²` from the S³ Green function).
  - `docs/antipodal_matter_interaction_synthesis_research_plan.md` /
    `docs/s_bam_vertex_generation_research_plan.md` — #139/#140, the matter
    sector and the antipodal Z₂ Ward identity the gauge vertex shares.
  - `docs/charge_conjugation_swap_research_plan.md` /
    `docs/cpt_assembly_research_plan.md` — #63/#64, the C-swap and CPT (the
    throat as the C-surface, `c₁ → −c₁`).
  - `docs/alpha_G_ledger_classification_research_plan.md` — #105, `α` as the
    universal residual (the 137 problem).

## Run

```
python -m experiments.closure_ledger.gauge_matter_coupling_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_gauge_matter_coupling_probe/`.
Expected verdict:
`GAUGE_MATTER_COUPLING_MINIMAL_Z2_SELECTED_CHARGE_CONSERVED_ALPHA_INPUT`, 8/8 PASS.

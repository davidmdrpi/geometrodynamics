# Gauge Ward identity and current conservation audit (PR #142)

> **Framing (to avoid a category error).** QFT on the *fixed classical* throat
> geometry — geometry → fields, **not** quantum gravity. "Gauge" is the
> U(1)_Hopf field (the photon of #42–#44), "matter" is the antipodal cavity
> modes (#129–#140); the Ward identity is read off the *classical* antipodal
> throat. The metric is never quantised.

PR #141 built the minimal gauge–matter coupling at the antipodal throat. This PR
audits its consistency: does the matter current stay conserved, does the
Ward–Takahashi identity hold, and is the photon protected from a mass? The
finding ties the gauge sector to the matter stability thread: **current
conservation, the Ward identity, and photon masslessness all follow from the
unitary antipodal throat (#129)** — the *same* postulate that gives stable
matter (#130/#135/#136/#138). An absorbing throat would leak charge and break
gauge invariance. So **unitarity ⟹ gauge invariance: one postulate, both.**

## The conserved Noether current

The global U(1)_Hopf phase symmetry gives the Noether current
`j^μ = i(ψ* ∂^μ ψ − (∂^μ ψ*) ψ)`, conserved on shell `∂_μ j^μ = 0`. For a
stationary cavity mode `ψ_n(r) e^{−iω_n t}` the charge density `ρ = 2ω_n |ψ_n|²`
is time-independent, so conservation reduces to the vanishing radial current.

## Current conservation at the antipodal throat — real modes carry no charge flux

The antipodal cavity modes are **real** (#135): a self-adjoint, unitary-mirror
operator has a real eigenbasis. The radial charge current is therefore

```
j^r ∝ Im(ψ_n* ∂_r ψ_n) = 0   (ψ_n real),
```

exactly zero (verified). No charge flows through the throat: the charge is static
and conserved — a stable charged particle. **Current conservation in the gauge
sector IS the zero-flux unitary-mirror property** (#129).

| mode | ω | max\|j^r\| | charge flux |
|---|---|---:|---|
| antipodal n=0 (real) | 1.166 | 0 | none (conserved) |
| antipodal n=1 (real) | 3.250 | 0 | none (conserved) |
| antipodal n=2 (real) | 5.367 | 0 | none (conserved) |
| absorbing fundamental (complex) | 1.893 − 1.12i | ≈ 0.014 | leaks into horizon |

## An absorbing throat would break it

An absorbing (ingoing) horizon gives **complex** quasinormal modes (#130); their
radial current `j^r = Im(ψ* ∂_r ψ) ≠ 0` (verified) carries charge **into** the
horizon. Current conservation fails, and with it gauge invariance — a charged
black-hole-style throat is not gauge-consistent. So **gauge invariance REQUIRES
the antipodal (unitary) throat** — the same requirement as stable matter.

## The Ward–Takahashi identity and photon masslessness

Current conservation implies the Ward–Takahashi identity

```
q_μ Γ^μ(p, p') = S⁻¹(p') − S⁻¹(p),
```

tying the gauge vertex (#141) to the matter inverse propagator (#135): the gauge
coupling is fixed by the matter dynamics, the hallmark of gauge invariance (its
differential form `Γ^μ = ∂S⁻¹/∂p_μ` normalises the vertex to the charge). It also
makes the vacuum polarisation transverse, `q_μ Π^μν = 0`, so the gauge correction
generates **no photon mass**: the `1/q²` photon (#42–#44) is protected. The Ward
identity protects masslessness.

## One postulate, both

Current conservation, the Ward identity, and photon masslessness all follow from
the unitary antipodal throat (#129) — the same real, self-adjoint, zero-flux
structure that gave the stable spectrum (#130), the unitary reciprocal propagator
(#135), the stable self-energy (#136), and the bounded vacuum (#138). **Gauge
invariance is the gauge face of the unitary mirror, not an extra assumption.**
Only the coupling strength `α` (#105) stays input.

## Tests

| # | test | finding |
|---|---|---|
| T1 | goal | gauge Ward identity and current conservation audit |
| T2 | Noether current | `j^μ` conserved; stationary mode ⟹ `ρ` static |
| T3 | conservation at the throat | real modes ⟹ `j^r = 0` (no charge flux; #129) |
| T4 | absorbing breaks it | complex modes ⟹ `j^r ≠ 0` ⟹ requires the antipodal throat |
| T5 | Ward–Takahashi | `q_μ Γ^μ = S⁻¹(p') − S⁻¹(p)` (#141 ↔ #135) |
| T6 | transversality | `q_μ Π^μν = 0` ⟹ photon massless (`1/q²` protected) |
| T7 | one postulate | gauge invariance from the unitary mirror (#129) — same as stable matter |
| T8 | assessment | `GAUGE_WARD_IDENTITY_CURRENT_CONSERVED_PHOTON_MASSLESS_FROM_ANTIPODAL_MIRROR` |

## Established and open

  - **Established (BAM-native):** the matter current is conserved at the
    antipodal throat (real modes carry no charge flux, `j^r = 0` — the unitary
    mirror #129), the Ward–Takahashi identity ties the gauge vertex (#141) to the
    matter propagator (#135), and transversality keeps the photon massless
    (`1/q²` protected, #42–#44). All follow from the unitary antipodal throat —
    the same postulate as stable matter; an absorbing throat would leak charge
    and break gauge invariance.

  - **Does not / open:** the coupling strength `α` (#105/#108) stays input;
    higher-order Ward identities and the running of `α` are not addressed. The
    `α` (#105/#108), bulk-scale (#133), and flavor (#134) residuals stand.

## Cross-references

  - `docs/gauge_matter_coupling_research_plan.md` — #141, the minimal coupling
    audited here.
  - `docs/null_throat_boundary_conditions_research_plan.md` /
    `docs/antipodal_vs_absorbing_qnm_research_plan.md` — #129/#130, the unitary
    mirror (real modes, zero flux) vs the absorbing (complex, leaking) throat.
  - `docs/antipodal_horizon_exchange_kernel_research_plan.md` — #135, the matter
    propagator the Ward identity ties the vertex to.
  - `bam_exchange_kernel_probe.py` (PRs #42–#44) — the photon `1/q²` the Ward
    identity protects.

## Run

```
python -m experiments.closure_ledger.gauge_ward_identity_probe
```

Writes `probe.json` + `probe.md` under
`experiments/closure_ledger/runs/<UTC timestamp>_gauge_ward_identity_probe/`.
Expected verdict:
`GAUGE_WARD_IDENTITY_CURRENT_CONSERVED_PHOTON_MASSLESS_FROM_ANTIPODAL_MIRROR`, 8/8 PASS.

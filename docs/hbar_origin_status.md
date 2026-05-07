# ℏ-origin status — summary of the closure-cycle action probe sequence

Closes the ℏ-origin research thread opened in
`docs/hbar_origin_research_plan.md`. A four-probe sequence under
`experiments/closure_ledger/` localized the structural content of the
closure-cycle action quantum and identified what's open.

## Probe sequence

Four probes, each with JSON + markdown archives under
`experiments/closure_ledger/runs/<timestamp>_<probe>/`:

| # | probe                                  | sub-target | result |
|---|----------------------------------------|------------|--------|
| 1 | `closure_cycle_action_probe`           | P1 (Layer 1 integer) | **PASS**: Layer-1 ledger sums per species are exactly `(2, 4, 106) · 2π`. |
| 2 | `closed_orbit_radial_action_probe`     | P1 (Layer 2 integer) | **PASS at high n**: closed-orbit radial action `S_full(l, n) / 2π → (n + 1)`; WKB-to-exact convergence at high n explains the prior C1/D1 failures as ground-state WKB error. |
| 3 | `hard_wall_boundary_verification`      | P1 physical justification | **PASS**: Dirichlet at both endpoints verified numerically AND forced by `T² = −I` (T-fixed-point argument gives ψ = 0 at the throat). DD wins decisively over DN/ND/NN/soft+soft. |
| 4 | `aharonov_bohm_hopf_fibre_probe`       | Sub-target #2 (angular) | **PASS**: Hopf holonomy `π·cos(χ)` matches numerical integration to machine precision; spinor double-cover integer-quantized at exactly the polar fibres χ ∈ {0, π/2, π}; Hopf-throat partnership at χ = 0 contributes one full closure quantum. |
| 5 | `tangherlini_omega_m_e_probe`          | Sub-target #3 (dimensional bridge) | **PARTIAL**: Q1 holds at 5 % (suggestive Compton match for the electron only); Q2 fails decisively (lepton mass ratios not inside the Tangherlini ω-spectrum). |

## Headline finding — closure cycle is structurally complete

The closure cycle is integer-quantized in units of 2π for every
species. Combining all probes:

```
N_total(species)  =  k                          [antipodal closure]
                  +  1                          [Hopf-throat partnership at χ = 0]
                  +  100 · [k = 5]              [τ-uplift closure quantum]
                  +  Σ (n_i + 1)                [radial Layer-2 per coupled mode]
```

Concrete per-species integers (under the B2_radial_ladder coupling
that gives the cleanest WKB convergence):

| species  | k | antipodal | Hopf+throat | uplift | radial | **N_total** |
|----------|--:|----------:|------------:|-------:|-------:|-------:|
| electron | 1 | 1         | +1          | 0      | +1 (n=0) | **3** |
| muon     | 3 | 3         | +1          | 0      | +2 (n=1) | **6** |
| tau      | 5 | 5         | +1          | +100   | +3 (n=2) | **109** |

Each integer is a structural count of closure quanta (factors of 2π).
The closure-phase ledger at Layer 1 + Layer 2 closes in units of 2π
for any sensible species → (l, n) coupling under the **exact-quantum
reading** (hard-wall BS, throat-T² = −I, Hopf-AB at canonical χ).

## What's principled-bounded vs derived

- **Principled-bounded**: yes. The closure cycle is integer-valued
  in units of 2π for every species, with the four constituent
  channels (antipodal, Hopf, throat, β-uplift, radial BS) each
  integer-quantized individually. The geometric structure of the
  ledger fixes the integer counts up to a per-species mode coupling.

- **Derived**: no, in the strict sense. The closure cycle's
  integer pattern is structural; the **conversion factor to ℏ in
  physical units** is not. Predicting ℏ in SI requires R_MID
  (the throat radius) to be geometrically determined from a
  self-consistency condition — open as THESIS.md "Self-consistent
  throat radius" target.

## Structural sub-target #3 verdict

Probe 5 (`tangherlini_omega_m_e_probe`) checks two natural
identifications:

- **Q1**: `ω(l=1, n=0) = 1` in canonical R_MID = ℏ/(m_e c) units?
  Observed: ω = 1.0547 — 5.47 % deviation from 1. Suggestive
  Compton-frequency match at the order-of-magnitude level, but not
  predictive.
- **Q2**: Lepton mass ratios match Tangherlini ω-ratios? Observed:
  largest ω-ratio in the catalog is 3.87, while m_μ/m_e = 207 and
  m_τ/m_e = 3477. The lepton mass ladder is **NOT in the Tangherlini
  eigenfrequency spectrum** — it lives in the locked surrogate's
  closure-quantum uplift β·max(0, k−3)² and the pinhole γ structure
  identified in the prior γ-offset probe.

The dimensional bridge is therefore **suggestive at the order-of-
magnitude level for the electron only**. BAM remains
**dimensional-ratio-complete** (predicts m_μ/m_e and m_τ/m_e at
sub-percent through the locked surrogate) and
**dimensional-scale-incomplete** (the absolute MeV scale is anchored
at m_e, not derived).

## How this refines `docs/THESIS.md` "Where does ℏ enter?"

THESIS.md flagged two natural candidates: Aharonov-Bohm flux
quantization around the Hopf fibre, and discreteness of the
Tangherlini eigenvalue spectrum. The probe sequence has now
addressed both:

- The **Hopf AB form is verified**: `π·cos(χ)` matches to machine
  precision; spinor double-cover gives 2π·cos(χ); polar fibres
  χ ∈ {0, π/2, π} are integer-quantized. The Hopf-throat partnership
  contributes exactly one closure quantum at χ = 0. **Confirmed.**
- The **Tangherlini eigenvalue identification** with m_e c² holds
  approximately (5 %) for the electron only and fails to reproduce
  the species mass ratios. **Partial only.**

The closure-cycle action quantum picture is structurally complete
in geometric units; the SI-conversion gap is precisely identified.

## What's next (post-merge)

The natural follow-up is **NOT** another closure-ledger probe — the
closure-cycle integer-quantization is now established. The remaining
ℏ-origin sub-target is:

- **Sub-target #4: R_MID self-consistency.** Determine R_MID as
  the equilibrium throat radius for the locked mass spectrum, rather
  than imposing R_MID = 1 by convention. With R_MID geometrically
  determined, the conversion `ℏ = m_e R_MID c` becomes a prediction
  in physical units. This is **outside the closure-ledger
  experiment's scope**; it requires throat dynamics machinery not
  present in the current codebase, and is flagged in `docs/THESIS.md`
  "Open problems" as the deeper geometric-self-consistency question.

## Cross-references

- `experiments/closure_ledger/closure_cycle_action_probe.py` —
  Layer-1 integer counts.
- `experiments/closure_ledger/closed_orbit_radial_action_probe.py` —
  Layer-2 integer counts under hard-wall BS.
- `experiments/closure_ledger/hard_wall_boundary_verification.py` —
  T² = −I forces Dirichlet at the throat.
- `experiments/closure_ledger/aharonov_bohm_hopf_fibre_probe.py` —
  angular Hopf-AB channel verified.
- `experiments/closure_ledger/tangherlini_omega_m_e_probe.py` —
  dimensional bridge is partial; ℏ remains anchored.
- `docs/odd_k_closure_lemma.md` — the Layer-1 universality lemma
  that the integer counts refine.
- `docs/quark_beta_status.md` — companion summary for the quark
  β-derivation thread.
- `docs/hbar_origin_research_plan.md` — original research plan.

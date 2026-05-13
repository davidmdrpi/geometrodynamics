# From closed geometry to a single scale gap

*A reading of the BAM closure-ledger experiment as a four-step
argument about where ℏ enters.*

This note distils the closure-ledger sequence under
`experiments/closure_ledger/` — eight ℏ-origin probes plus the
follow-up transport/resistance derivation — into a single continuous
argument. The probe-level technical records remain in
`docs/hbar_origin_status.md`,
`docs/transport_resistance_research_plan.md`, and the timestamped
archives; here we tell the story.

The arc is four steps:

  §2. **Integer closure ledger.** The closure cycle is exactly
      integer-quantised in units of 2π for every species.

  §3. **R_OUTER fixed point.** The lepton spectrum cross-selects a
      single outer radius R* ≈ 1.262 to 0.008 %.

  §4. **γ, transport, resistance geometrized.** Every phenomenological
      parameter of the locked surrogate is a closure-quantum invariant
      or a Tangherlini grid quantity at R*.

  §5. **Remaining m_e / 1.054 scale gap.** The single residual handle
      between geometric units and the absolute MeV scale.

---

## §1. The setting

BAM provides three structural inputs and no others:

1. **A closed spatial S³** with great-circle action `2π` — the
   antipodal cavity that quantises closure.
2. **A non-orientable throat** with transport map `T = iσ_y`,
   satisfying `T² = −I` — the geometric source of the Z₂ partition
   class and the spinor double cover.
3. **A 5D Tangherlini bulk** in the radial direction with metric
   factor `f(r) = 1 − r_s² / r²`, giving discrete radial bound
   states at angular indices `l` and radial labels `n`.

These three pieces are connected by a single ledger: the total
phase accumulated along a closed worldline. The Layer-1 ledger
sums four channels — antipodal closure, Hopf holonomy, throat
phase, and the β-uplift that locks the heaviest shell — and the
prior closure-ledger experiment showed it closes mod 2π
universally for odd k (the odd-k lemma, `docs/odd_k_closure_lemma.md`).
What's new here is treating the closure cycle as an action and
asking whether that action is integer-valued in units of 2π for
each species individually, not just the sum mod 2π.

## §2. Closed geometry yields integer 2π action counts

The first three probes establish that BAM's closure cycle is
*structurally integer-quantised* — its phase per species is a
specific integer multiple of 2π, not just zero mod 2π.

**Layer 1 — closed-form channels** (probe 1,
`closure_cycle_action_probe`). Summing the four wired channels at
the locked baseline `(action_base = 2π, χ = 0, T², β = 50π)`:

```
Φ_total(k)  =  2π·k  +  π·cos(0)  +  π  +  50π · max(0, k − 3)²
            =  2π·(k + 1)  +  50π · max(0, k − 3)²
```

Per species: electron Φ = 4π, muon Φ = 8π, tau Φ = 212π. Dividing
by 2π gives the integer counts **(N_e, N_μ, N_τ) = (2, 4, 106)**.
The Hopf holonomy π and the throat phase π *partner* at χ = 0 to
form a single closure quantum 2π; the τ-uplift contributes 100·2π
from the integer-winding lock `4β = 100·(2π)`. Every species lands
exactly on an integer multiple of 2π — Layer 1 is closure-cycle
integer-quantised by construction.

**Layer 2 — the radial bulk** (probes 2 and 3,
`closed_orbit_radial_action_probe` and
`hard_wall_boundary_verification`). The closed-orbit Bohr-Sommerfeld
action for the (l, n) bound mode satisfies

```
S_full(l, n) / (2π)  →  (n + 1)   as n → ∞
```

— hard-wall BS with the eigensolver's 0-indexed `n` mapping to
the standard BS quantum number `N = n + 1`. The hard-wall identification
is not a numerical artefact: the throat condition `ψ = T ψ` combined
with `T² = −I` forces `ψ = 0` at the throat by the T-fixed-point
argument. Numerically, the eigenfunctions vanish at both grid
endpoints to machine precision; the alternative boundary conditions
(DN, ND, NN, soft+soft) all fail by more than 25 % at high n,
while DD matches to better than 0.001 % at n ≥ 2.

The earlier closure-ledger Layer-2 failures (C1 at 0.326 rad, D1
at 0.577 rad) are now isolated to *WKB-at-ground-state* error: the
exact action is `(n + 1)·2π` per bound mode; the WKB approximation
at n = 0 underestimates it by 10–25 %. Probes that used WKB ground-
state actions saw the underestimate as a closure-spread residue.

**Angular channel — Aharonov-Bohm** (probe 4,
`aharonov_bohm_hopf_fibre_probe`). The Hopf connection `A = ½ cos(χ) dφ`
produces a holonomy `π·cos(χ)` per fibre loop, verified numerically
to machine precision. Under the spinor double cover (two traversals),
this becomes `2π·cos(χ)` — integer-quantised at exactly the three
polar fibres χ ∈ {0, π/2, π}. At χ = 0 (the canonical fibre, where
the BAM lock sits), the Hopf-throat partnership contributes a
single closure quantum.

**Result.** Per species the closure cycle decomposes as

```
N_total(species)  =  k                          [antipodal: k·2π]
                  +  1                          [Hopf-throat partnership at χ = 0]
                  +  100 · [k = 5]              [τ-uplift closure quantum]
                  +  Σ_{(l, n)} (n + 1)         [radial Layer-2]
```

with the cleanest coupling B2_radial_ladder giving:

| species | k | coupled mode | N_total |
|---|---:|---|---:|
| electron | 1 | (l=1, n=0) | **3** |
| muon | 3 | (l=1, n=1) | **6** |
| tau | 5 | (l=1, n=2) | **109** |

The closure cycle is integer-valued in units of 2π for every
species. The τ-uplift quantum 100 reappears in 109 − 3 − 6 = 100 —
the same closure-quantum integer that anchored the lepton sector
all along.

## §3. The geometry self-selects its outer radius

The closure cycle is integer-quantised in *geometric* units. To
convert to physical units (cm, ℏ, s) we need a length scale. The
canonical identification

```
R_MID · m_e c² / ℏ  =  ω(1, 0)
```

ties the throat radius to the electron Compton scale through the
lowest Tangherlini eigenfrequency. If `ω(1, 0) = 1` exactly, the
throat radius equals the reduced Compton wavelength `λ̄_C = ℏ/(m_e c)`
and the dimensional bridge closes: `ℏ = m_e R_MID c` with no fit
parameters.

This is the *Compton bridge condition*. The question is whether it
is geometrically realisable.

**Two natural conditions on R_OUTER** (probe 6,
`tangherlini_electron_scale_bridge`). On the canonical R_MID = 1
grid:

| condition | R_OUTER | ω(1, 0) | Σ V_max[0..5] |
|---|---:|---:|---:|
| Compton bridge `ω = 1` | 1.449005 | 1.000000 | 23.6308 |
| γ-lock `Σ V_max = 22.5` | 1.262266 | 1.053694 | 22.5000 |

These differ by 14.79 % and cannot both hold — the Tangherlini
spectrum and the barrier-sum γ are *both* monotonic in R_OUTER, but
they pass through the natural target values at different geometries.

**The Compton bridge is physically vetoed** (probe 7,
`compton_bridge_feasibility_probe`). Re-running the locked lepton
surrogate at the Compton-bridge γ = 23.63 breaks both μ and τ by
~46 %. Sweeping β over `{30π, …, 200π}` finds no value recovering
both species at sub-percent. The best joint fit has 42 % error;
the best μ-only fit leaves τ at 42 %; the best τ-only fit leaves
μ at 63 %. **The Compton geometry is mathematically definable but
incompatible with the locked lepton spectrum.**

**The γ-lock is selected by cross-species self-consistency**
(probe 8, `R_outer_self_consistency_probe`). Bisecting independently
on each species' mass error:

| anchor | R* | residual on the other species |
|---|---:|---:|
| muon (zero of err_μ) | **1.262239** | err_τ = +0.16 % |
| tau (zero of err_τ) | **1.262338** | err_μ = −0.16 % |

The two species select the *same* R_OUTER to within `|ΔR*|/R* =
0.008 %`. At the fixed point, Σ V_max[0..5] = 22.4994 — recovering
the locked γ_lepton value to better than 0.01. A single geometric
R_OUTER simultaneously fits both lepton mass ratios; the γ used
in the surrogate is not fit but is the radial barrier-sum on the
same grid.

This is not trivial. γ is one parameter, and m_μ/m_e and m_τ/m_e
are two independent anchors; the cross-species agreement at 0.008 %
is the geometry telling us *which* R_OUTER it wants.

The fixed point is structurally selected: phase_per_pass perturbations
shift R* by less than 10⁻⁴ % (decoupled), `resistance_scale` by up to
3 % at ±5 % perturbation, and `transport_strength` by up to 7 % at
±1 % perturbation. The closure-phase channel is completely decoupled;
the remaining sensitivities live in the cross-shell mixing parameters
`(transport, resistance)` of the locked surrogate. §4 closes those
channels.

## §4. γ, transport, and resistance geometrized

The R_OUTER fixed point of §3 uses the locked surrogate with γ, T, ρ
playing the role of phenomenological inputs. The pinhole-origin chain
(`pinhole_origin_probe`) already established the structural reading
for γ:

```
γ_lepton  ≈  Σ_{l = 1..5} V_max(l)  =  22.0    (canonical Chebyshev grid)
```

— the same radial barrier-sum operator that the QCD residual sector
locks at γ_quark = 22.25. The pinhole is structurally the **barrier-
sum on a Tangherlini grid**, not a free parameter.

The transport / resistance derivation thread
(`docs/transport_resistance_research_plan.md`) closes the remaining
two channels. Two probes:

**Origin probe (`transport_resistance_origin_probe`).** Scans six
categories of candidates for each parameter:

| parameter | locked value | best reading | %Δ |
|---|---:|---|---:|
| `transport_strength` | 25.1 | `8π = 4·(2π)` | +0.13 % |
| `resistance_scale`   | 0.2179 | `7π / 100` *or* `4·(ω(1,0) − 1)` | +0.94 % / +0.48 % |

The transport identification is decisive — `8π` is the 4th closure
quantum, structurally identical in form to the antipodal k·2π, the
Hopf-throat 1·2π, and the τ-uplift 100·2π. Resistance has two
within-1 % candidates that the origin probe alone cannot tell apart.

**Disambiguation probe (`resistance_disambiguation_probe`).**
Re-bisects R_OUTER under each resistance reading, paired with
`transport = 8π`:

| reading                          | R*_μ      | cross-species | R*-match to locked | γ at R* |
|----------------------------------|----------:|--------------:|-------------------:|--------:|
| locked baseline (control)        | 1.262239  | 0.0078 %      | (ref)              | 22.499  |
| closure-quantum (`8π`, `7π/100`) | 1.262636  | **0.0021 %**  | **0.031 %**        | 22.508  |
| eigenfrequency (`8π`, `4·(ω−1)`) | 1.258316  | 0.0025 %      | 0.311 %            | 22.417  |

Both readings produce TIGHTER cross-species agreement than the locked
baseline — they are cleaner mathematical objects than the fitted
values. The closure-quantum reading wins decisively on the
discriminating criterion (R*-match to the locked baseline, γ-match to
canonical 22.5): an order of magnitude tighter on each axis. The
eigenfrequency reading is ruled out: the numerical coincidence
`0.218 ≈ 4·(1.054 − 1)` is not preserved when R_OUTER varies through
the bisection.

The full closure-quantum ledger of the locked lepton surrogate is now:

| parameter             | locked value | structural identification |
|-----------------------|-------------:|---------------------------|
| `action_base`         | 2π           | S³ great-circle action     |
| `transport_strength`  | 25.1         | 8π = 4·(2π)               |
| `resistance_scale`    | 0.2179       | 7π / 100                  |
| `pinhole γ`           | 22.5         | Σ V_max[1..5] ≈ 22.0 (≈ 7π) |
| `β` (τ-uplift)        | 50π          | locked closure quantum     |
| `4β` (τ-uplift integer) | 100·2π     | τ-uplift quantum           |

Every parameter is a closure-quantum invariant or a Tangherlini grid
quantity at R*. With this reading the R_OUTER self-consistency loop
closes on principled inputs alone — no fitted transport or resistance
constants enter. R_OUTER is **structurally selected by the BAM
closure-quantum scaffolding**.

## §5. The remaining m_e / 1.054 scale gap

After §2, §3, and §4, the BAM framework predicts:

* **Dimensionless ratios** — `m_μ/m_e ≈ 207, m_τ/m_e ≈ 3477` at
  sub-percent through the locked surrogate, with R_OUTER selected
  by the self-consistency loop, not fit, and every surrogate
  parameter reduced to a closure-quantum invariant.
* **Integer closure quanta** — `(N_e, N_μ, N_τ) = (3, 6, 109)`
  under the B2_radial_ladder coupling, with the τ-uplift integer
  100 emerging as the residual.

What it does *not* yet predict is the absolute MeV scale. The
electron mass is anchored externally; the conversion factor from
geometric units to physical ℏ involves a single residual number:

```
ℏ · ω(1, 0)  =  1.054 · m_e c²    at the self-selected R_OUTER ≈ 1.262
```

The 1.054 is **structural**, not fitted — it is ω(l = 1, n = 0)
evaluated at the cross-species fixed point R_OUTER ≈ 1.262, where
the lepton mass ratios both fall into place. The Compton-bridge
alternative ω = 1 exactly was tested and physically vetoed.

The 1.054 is the single remaining dimensional handle. If it has a
closed form in `(k_5, π, barrier-spectrum invariants)`, the
dimensional bridge closes:

```
ℏ  =  m_e c² / ω(1, 0)  ·  (factor 1.054)
   =  m_e c² · (1 / 1.054) · ...
```

would predict ℏ in physical units from `m_e` alone. If it doesn't —
if 1.054 is an irreducible numerical artefact of the Tangherlini
spectrum at R_OUTER ≈ 1.262 with no closed form — then BAM is
*dimensional-ratio-complete and dimensional-scale-incomplete*: it
predicts every ratio at sub-percent but anchors the absolute MeV
scale at `m_e`.

This is the next concrete target, and it is a sharp one. The
1.054 is computable to arbitrary precision from the canonical
Tangherlini grid; closed-form candidates are constructible by
exhaustive search over small `(k_5, π, integer)` combinations.
A targeted probe in this branch begins that search.

---

## What this leaves open

After the §4 closure-quantum reading of `(transport, resistance)`,
two structural questions remain:

* **Closed form for the 1.054 factor.** The factor-1054 search probe
  enumerated small `(k_5, π, integer)` combinations and returned a
  clean negative result. Whether a larger ansatz space (e.g.
  barrier-spectrum invariants combined with closure-quantum integers)
  contains a closed-form match is open. If 1.054 has no closed form,
  BAM is *dimensional-ratio-complete and dimensional-scale-
  incomplete* with the m_e anchor as the unique external input.

* **The deeper R_MID self-consistency.** R_MID = 1 by convention.
  Determining R_MID dynamically from a self-consistency condition
  (equilibrium throat radius for the locked mass spectrum) would
  lift the m_e anchor itself. This is THESIS.md "self-consistent
  throat radius" — outside the closure-ledger scope.

---

## Cross-references

- `docs/hbar_origin_status.md` — eight-probe ℏ-origin summary table
  with archive pointers.
- `docs/hbar_origin_research_plan.md` — the research-plan document
  that opened the ℏ-origin thread.
- `docs/transport_resistance_research_plan.md` — §4 derivation
  thread (closure-quantum readings of transport and resistance).
- `docs/odd_k_closure_lemma.md` — the Layer-1 universality lemma
  that the integer counts refine.
- `docs/quark_beta_status.md` — companion summary for the quark
  β-derivation thread.
- `experiments/closure_ledger/runs/` — per-probe JSON + markdown
  archives.

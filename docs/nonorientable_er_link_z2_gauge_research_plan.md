# The Z₂ link twist: non-orientable ER links as a gauge field over the network, and what it does to the antipodal wave interaction (PR #227 — research plan)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity. Nothing
> below quantises the metric; the twist class is a discrete label of the
> classical background, and every claim is a statement about matter
> fields propagating on it.

## 0. The gap this closes

Two mature arcs in this repo have never been joined.

**Arc A — the antipodal wave interaction.** The antipodal identification
(#128/#129) fixes the throat boundary data (even-`l` Neumann, odd-`l`
Dirichlet), grades the matter exchange kernel by `(−1)^l` (#135), selects
the cubic and quartic vertices by `Σl` even (#137/#138), and drives
nucleation through the antipodal caustic (#166, #175). The synthesis
(#139) calls this **Thread A**, "the antipodal Z₂ `(−1)^l`" — and treats
it throughout as a **single global parity** attached to one throat.

**Arc B — the bulk ER-link network.** Bridges pair mouths (#206), surgery
composes them (#207), tri-mouth junctions carry the GHZ valence (#208),
the two-port throat projects to the advanced kernel (#216/#217), and the
two-throat network carries a genuine mouth-to-mouth two-level beat
(#223/#224). Orientation appears here only as an **inert decoration**:
`NetworkMouth.orientation` is a ±1 field multiplied into `U_BA` and
`traverse_throat`, never varied, never given a consequence
(`geometrodynamics/transaction/network.py:64,130,282,404`).

Meanwhile a third arc — the non-orientable hadron sector (#100–#114) —
uses a Z₂ orientation label with real dynamical teeth: the Möbius flux
tube's antiperiodic boundary condition shifts integer → half-integer
modes and lifts the glueball tower by `πσ` in `M²`.

**The gap.** Thread A's Z₂ is a *global parity of one throat*. Arc B is a
*network of many links*. A Z₂ that lives on a network is not a parity —
it is a **gauge field**, with link variables, gauge orbits, and Wilson
loops, and only its **holonomies** are physical. Nobody has written that
down. This plan does, and then computes what it does to Arc A.

## 1. The object: the link twist ε ∈ H¹(network; Z₂)

Assign to each ER link `b` of the network a **twist class** `ε_b ∈ {±1}`,
and to each mouth `a` a frame sign `s_a ∈ {±1}`. The transit amplitude
through link `b` carries `ε_b`; a mouth reframing acts as

```
ε_b  ⟶  s_{a(b)} · ε_b · s_{a'(b)}        (Z₂ gauge transformation)
```

so individual `ε_b` are pure convention and only the **loop holonomy**

```
W(γ)  =  ∏_{b ∈ γ} ε_b  ∈ {±1}
```

is physical. This is Z₂ lattice gauge theory whose Wilson loops are the
orientation holonomies of the network — `w₁` promoted from a per-loop
label (#121's sector-phase ledger already names it `w₁ ∈ H¹(loop;Z₂)`)
to a **connection over the ER graph**.

### 1.1 Precision note — where the non-orientability actually lives

This must be stated carefully, because the loose version is false. The
repo's own #183 result records the two deck determinants: the bulk
quotient `S³/antipodal = RP³` has deck det **+1** (orientable), the brane
angular slice `S²/antipodal = RP²` has deck det **−1** (non-orientable).
**RP³ is an orientable 3-manifold; `w₁(RP³) = 0`.** So the twist is *not*
ambient non-orientability of the exterior. It is:

  - **`H¹(RP³; Z₂) = Z₂`** — the two **flat Z₂ line bundles** over the
    quotient (equivalently its two spin structures, the choice #B2 already
    makes when it selects `T² = −I`); and
  - **`w₁(RP²) ≠ 0` at the mouth** — the genuinely non-orientable RP²
    cross-cap that carries Pin⁻ (#169/#170), whose transport is the
    orientation-reversing Hopf map `σ(z₁,z₂) = (z̄₂, −z̄₁)` ⟹ `T = iσ_y`.

Both deliver the *same* ±1 holonomy on a transit loop, which is why one
Z₂ suffices; but the plan says "twist class", not "non-orientable
spacetime", and any write-up that conflates them should be corrected.

## 2. The five claims — and their pre-flight validation

Each claim below was checked in a reduced model before this plan was
written (`experiments/closure_ledger/nonorientable_er_link_z2_preflight.py`,
runs in <5 s). These are **not** the probe: they are minimal models that
establish the claims are not vacuous and fix the expected signs and
magnitudes. The probe (§3) must reproduce them on the repo's real
geometry.

### Claim 1 — the two flat bundles ARE the #129 boundary dichotomy (derivation, not postulate)

Under the deck map, `Y_l(−x) = (−1)^l Y_l(x)`. Sections of the **trivial**
flat bundle over the quotient are the `ε = +1`, **even-`l`** modes;
sections of the **twisted** bundle are the `ε = −1`, **odd-`l`** modes.
So #129's "even-`l` Neumann, odd-`l` Dirichlet" is not two boundary
conditions imposed side by side — it is **the two flat Z₂ bundles of
`H¹ = Z₂`, one per twist sector**. Thread A's `(−1)^l` is the `ε = −1`
Wilson line.

*Status if it survives:* one item moves in the #139 epistemic ledger from
**postulated** ("the antipodal identification itself") to **classified**
(the identification is a choice of flat bundle, and `H¹(RP³;Z₂) = Z₂`
says there are exactly two — the dichotomy is forced, only the branch is
chosen). It does **not** make the antipodal postulate "forced"; it makes
its *structure* forced.

### Claim 2 — a twisted link freezes mouth exchange **exactly** (the sharp new result)

#224's two-throat network is a **ring**: mouths A and B are joined by two
exterior brane arcs, so the doublet coupling is the sum of two paths,
`J = J₁ + J₂`, and the splitting is `Δω = 2|J₁ + J₂|`. A twist on one arc
makes the loop holonomy `W = −1` and the paths **subtract**:

| arcs | `W` | splitting `Δω` |
|---|---:|---:|
| `J₁ = J₂ = 1` | +1 | **4.000000** |
| `J₁ = J₂ = 1` | −1 | **0.000000** (machine zero) |
| `J₁ = 1, J₂ = 0.7` | +1 | 3.400 |
| `J₁ = 1, J₂ = 0.7` | −1 | 0.600 |

On the **symmetric** network a single twisted link sends the doublet
splitting to **exactly zero**: `T_exchange = π/Δω → ∞`. Mouth-to-mouth
transfer is **topologically frozen** — not suppressed, not slow, but
degenerate by destructive interference of the two exterior arcs.

This is a genuinely new mechanism and it is the *complement* of #224's
two: #224 froze exchange **dynamically** (the `α⁻⁴` coupling law at the
primordial anchor) and **by asymmetry** (the Rabi/Lorentzian localization,
a 1.5% bias pinning survival above 98.6%). Both are approximate and both
degrade as the network becomes symmetric — #224 explicitly names the
"exact-degeneracy loophole" (complete transfer in a perfectly symmetric
network) and closes it only by appealing to no two throats being
identical. **The twist closes it in the opposite limit**: the more
symmetric the network, the more exact the topological freeze. The two
mechanisms cover complementary regimes, which is the structural payoff.

### Claim 3 — Z₂ gauge invariance is exact, and a twisted loop is an antiperiodic ring

Pre-flight, `N = 6` ring with random bond strengths: four random gauge
transformations `ε_b → s_a ε_b s_{a'}` leave the spectrum invariant to
**0.00e+00** (exactly — sign flips are similarity transforms), while a
genuine `W = −1` twist moves it by **0.55**. Only the holonomy is
physical, exactly as a gauge structure requires.

And the twisted ring is the **antiperiodic** ring. Pre-flight, `M = 240`
sites, mode index `kM/2π`:

```
untwisted:  0.000, 1.000, 1.000, 2.000, 2.000, 3.000   (integer)
twisted  :  0.500, 0.500, 1.500, 1.500, 2.500, 2.500   (half-integer)
```

**This is the same integer → half-integer shift that gives the Möbius
glueball tower its `+πσ` lift in `M²` (#100) and the Möbius baryons their
`2√σ` gap (#103).** The unification claim: the hadron sector's Möbius
tower and the ER network's frustrated loops are **one Z₂ holonomy in two
realizations**. If the probe confirms the shift is quantitatively the
same object, the non-orientable experimental note (#110/#114 — π₁(1600),
η₁(1855), Λ_c 3135, Λ_b 6469, the 849 MeV dipion endpoint) becomes the
**experimental face of the network twist**, which is the first
observational handle any part of Arc B has had.

### Claim 4 — on the quotient, the antipodal revival period halves, and the twist sets its sign

The strongest new result on Arc A's own turf. #166 computes the conformal
zonal wave on **S³**: exact mirror-refocus `ψ(χ, πR) = −ψ(π−χ, 0)` at
`t = πR`, full revival at `t = 2πR`. On the **quotient** the mode set is
restricted to one twist sector, and the tower `ω_l = (l+1)/R` then
collapses the revival. Pre-flight (`l ≤ 160`, Gaussian packet, errors are
relative L² norms):

| sector | `t = πR` | `t = 2πR` |
|---|---|---|
| full S³ tower | `f = −P f₀` (2e-14) — #166 reproduced | `f = +f₀` (0) |
| **untwisted (even-`l`, ε=+1)** | **`f = −f₀`** (0, exact) | `f = +f₀` (0) |
| **twisted (odd-`l`, ε=−1)** | **`f = +f₀`** (0, exact) | `f = +f₀` (0) |

Two things follow, both new:

1. **The revival period halves on the quotient**: `πR`, not `2πR`. The
   antipodal focus of #166 does not land on a far-away antipodal point —
   on the quotient the antipode **is** the emitter, so the caustic returns
   *to its own source*, and it does so a full factor of two sooner than
   the S³ picture suggests.
2. **The twist class sets the revival sign**: `ε = +1` returns `−f₀`
   (destructive against the source), `ε = −1` returns `+f₀`
   (constructive). Self-focusing is **constructive precisely in the
   twisted, odd-`l` sector** — which is the sector #202's Pin parity
   already selects for the electron (`k = 1`, odd at the neck).

The consequence to test is the nucleation threshold. #166/#58 set pair
creation at `2 m_e c²` from a caustic gain computed on S³. If the
returning focus is coherent with the source in the twisted sector, the
amplitude at the focus adds rather than passing through, and the energy
density at the caustic goes as the *square* of that sum. **The
conjecture to test is a factor-of-4 energy advantage for the twisted
sector, and exact cancellation at the focus for the untwisted one** — a
topological selection rule on which sector can nucleate at all. This is
the claim in this plan most likely to break, and §4 says how it breaks.

### Claim 5 — the vertex selection rule generalises to `Σl ≡ n_twist (mod 2)`

#137/#138 derive `Σl` even for the cubic and quartic vertices. That is the
`n_twist = 0` case. With legs carrying twist classes `ε_i`, deck
invariance of the integrand requires `(−1)^{Σl} = ∏_i ε_i`, i.e.

```
Σl  ≡  n_twist   (mod 2)          [n_twist = #{legs with ε_i = −1}]
```

Pre-flight enumeration over `l_i ≤ 3`: exactly 32 of 64 triples survive
for every `n_twist`, with `Σl` parity `{0}` for even `n_twist` and `{1}`
for odd. **A vertex with an odd number of twisted legs selects `Σl`
odd — the exact complement of #137's rule.** The Z₂ never disappears; it
moves. Immediate corollaries to check:

  - **Tri-mouth junctions (#208).** #208 derived a charged GHZ no-go by
    flux conservation (`k_A + k_B + k_C = 0` has zero solutions for a ±1
    doublet) and concluded GHZ can live only in the transported-frame
    spin label. The twist adds a *second, independent* junction rule:
    `ε_A ε_B ε_C = +1`. Junctions with an odd number of twisted legs are
    frustrated. This is a new selection rule in the sector #208 left open.
  - **Bridge surgery (#207).** The composition law `φ₁₄ = φ_a + φ_b + φ_c`
    should acquire a Z₂ term `+ π·[n_twist mod 2]` — i.e. surgery across
    an odd number of twisted bridges swaps the Bell outcome to its
    orthogonal partner. Machine-checkable against the existing 16-dim QM
    swapping calculation.

## 3. The probe (T1–T9)

`experiments/closure_ledger/nonorientable_er_link_z2_probe.py`, estimated
2–4 min. Each test must run on the **repo's real machinery**, not the
reduced pre-flight models.

| # | test | uses | pass criterion |
|---|---|---|---|
| T1 | goal + the object: link variables, gauge orbit, Wilson loops | — | statement |
| T2 | flat-bundle classification ⟹ the #129 dichotomy | `null_throat_boundary_conditions_probe` BC data | even-`l`/odd-`l` branches reproduced as the two bundles |
| T3 | Z₂ gauge invariance on the network | `transaction/network.py` (`orientation` promoted to a link variable) | spectrum invariant under `ε_b → s ε_b s'` to ~1e-14; only `W(γ)` moves it |
| T4 | **twisted-link exchange freeze** | the #224 two-throat ring | `Δω(W=−1)/Δω(W=+1) < 1e-6` on the symmetric network; `T_exchange → ∞` |
| T5 | twisted loop = antiperiodic ⟹ half-integer comb | ring solver + `qcd/topology.make_mobius_tube` | mode index shifted by exactly ½; matches the `+πσ` Möbius lift of #100 |
| T6 | **quotient revival at `πR`, sign = twist class** | the #166 conformal zonal solver | `f(πR) = −f₀` (ε=+1) and `+f₀` (ε=−1), both to ~1e-13 |
| T7 | caustic gain and the nucleation threshold per sector | #166 caustic + #58 threshold | measured gain ratio between sectors, reported honestly (see §4) |
| T8 | vertex rule `Σl ≡ n_twist (mod 2)`; junction rule `∏ε = +1`; surgery `+π` term | `cubic_vertex_ledger_probe`, `bridge_surgery_swapping_probe` | selection rules reproduced; surgery law gains the Z₂ term to ~1e-12 |
| T9 | verdict + honest ledger | — | see below |

**Verdict classes.**

  - `ER_LINK_TWIST_IS_A_Z2_GAUGE_FIELD_FREEZING_EXCHANGE_AND_GRADING_THE_ANTIPODAL_KERNEL`
    — T2–T6, T8 pass; T7 reported either way. The full claim.
  - `..._GAUGE_STRUCTURE_ONLY` — T3–T5 and T8 pass but T6/T7 fail: the
    network Z₂ is real, the Arc A consequences are not. Still publishable
    as the gauge-structure half.
  - `..._REFUTED` — T3 fails (the twist is gaugeable away on the real
    network, i.e. the ER graph is a tree with no independent loops). This
    is a real possibility and §4 treats it as the primary falsifier.

## 4. Honest scope — what is conjecture and how each piece dies

  - **Pre-flight ≠ probe.** Every number in §2 comes from a reduced model
    (two-level ring, tight-binding chain, zonal mode sums). They fix
    signs and magnitudes and prove the claims non-vacuous. They do **not**
    establish the results on the Tangherlini bulk, the MTY bridge, or the
    real two-throat network. The probe must redo all of it.
  - **Claim 2 is the safest.** It is two-path interference on a ring; the
    only way it fails is if the #224 network turns out not to be a ring
    with two independent arcs, in which case there is no loop to twist and
    the result is empty rather than wrong.
  - **Claim 4 is the riskiest** and should be expected to be trimmed. The
    mode-sector revival identities are exact and will reproduce. The
    *threshold* consequence (a factor-of-4 energy advantage, exact
    cancellation in the untwisted sector) rests on the returning focus
    being coherent with the source over the packet's whole support, which
    the reduced zonal model cannot settle — the caustic is regularized by
    the `ℓ_max ~ R/R_MID` cutoff, and the cutoff may wash out the phase
    relation before it does any work. If the measured gain ratio is not
    4, **report the measured number**; the revival identities stand
    regardless, and Claim 4 degrades to "period halves, sign is
    topological" without the threshold corollary.
  - **The primary falsifier is topological, not numerical.** A Z₂ gauge
    field has content only if the ER graph carries independent loops
    (`b₁ > 0`). #207's bridges pair mouths into a **perfect matching** —
    and a perfect matching is a *forest*, `b₁ = 0`, every `ε_b` gaugeable
    to +1. **The twist is physical only where the network has cycles**:
    the two-arc exterior ring of #223/#224, the tri-mouth junction of
    #208, and the surgered chains of #207. T3 must measure `b₁` on each
    and report it. If `b₁ = 0` everywhere in the repo's actual networks,
    this whole plan collapses to a relabeling, and that is the honest
    outcome to publish.
  - **What this cannot do.** It does not derive the antipodal postulate
    (it classifies it), does not touch the couplings `λ₃, λ₄`, does not
    reduce any of the four standing residuals (`ε`, `n_part`, `α`,
    `√σ/m_e`), and adds no dimensionful input — the twist is an integer
    mod 2.
  - **`ℏ` enters nowhere**, and "probability" remains the classical
    normalized field-norm fraction, as everywhere in the program.

## 5. Staging

The claims are independent enough to land separately, in increasing risk:

1. **§2 Claims 1, 3, 5** — classification, gauge invariance, selection
   rules. Algebraic, low risk, and they alone justify the object.
2. **Claim 2** — the exchange freeze on the #224 network. The headline.
3. **Claim 4** — the quotient revival and the threshold question. Highest
   risk, highest payoff; safe to defer to a successor PR if T6 passes and
   T7 is inconclusive.

## 6. Cross-references

  - `docs/antipodal_matter_interaction_synthesis_research_plan.md` — #139,
    Thread A (the global `(−1)^l` this plan promotes to a gauge field)
  - `docs/null_throat_boundary_conditions_research_plan.md` — #129, the BC
    dichotomy Claim 1 reclassifies
  - `docs/antipodal_horizon_exchange_kernel_research_plan.md` — #135, the
    parity-graded kernel
  - `docs/antipodal_focusing_threshold_research_plan.md` — #166, the S³
    refocus Claim 4 lifts to the quotient
  - `docs/nonlinear_antipodal_focusing_pde_research_plan.md` — #175, the
    nucleation gate at the antipodal node
  - `docs/mouth_exchange_dynamics.md` — #224, the two-level beat Claim 2
    freezes; the "exact-degeneracy loophole" it names
  - `docs/bridge_dressing_network.md` — #223, the network ring and comb
  - `docs/multi_mouth_bridge_ghz.md` — #208, the junction rule Claim 5 adds to
  - `docs/bridge_surgery_entanglement_swapping.md` — #207, the composition
    law and the perfect-matching (`b₁ = 0`) falsifier
  - `docs/glueball_closed_flux_loop_research_plan.md` — #100, the Möbius
    `+πσ` tower Claim 3 identifies with the twisted loop
  - `docs/pin_rp2_fermi_statistics_research_plan.md` — #170, Pin⁻ on the
    RP² mouth; `docs/geon_statistics_adjudication.md` — #196
  - `docs/bam_sector_phase_ledger_research_plan.md` — #121, which already
    names `w₁ ∈ H¹(loop;Z₂)` as the discrete Z₂ and proves it does not
    double-count against the U(1) holonomy
  - `geometrodynamics/transaction/network.py` — `NetworkMouth.orientation`,
    the inert decoration this plan makes dynamical
  - `geometrodynamics/embedding/transport.py` — the orientation-reversing
    Hopf map giving `T = iσ_y`

## 7. Reproduce (pre-flight)

```bash
python -m experiments.closure_ledger.nonorientable_er_link_z2_preflight
```

Prints the five checks of §2. The probe of §3 is not yet written; this
document is the plan for it.

# The Z₂ link twist: non-orientable ER links as a gauge field over the network, and what it does to the antipodal wave interaction (PR #227)

> **STATUS — EXECUTED.** The probe
> (`experiments/closure_ledger/nonorientable_er_link_z2_probe.py`, 8/8 PASS,
> ~1 min) has been run on the repo's real machinery. Verdict:
> `ER_LINK_TWIST_IS_A_Z2_GAUGE_FIELD_THAT_FREEZES_MOUTH_EXCHANGE_AND_HALVES_THE_QUOTIENT_REVIVAL_BUT_LEAVES_THE_NUCLEATION_THRESHOLD_UNTOUCHED`.
>
> **Three claims confirmed and strengthened** (Claims 1–3: the gauge
> structure, the exchange freeze — which turned out to be a *theorem*, not
> a two-level effect — and the Möbius identification, which turned out to
> be *already in the repo* under another name).
>
> **Two claims corrected by measurement, not defended** (§2.4 and §2.5).
> The conjectured factor-of-4 energy advantage at the focus is **wrong** —
> linear evolution conserves energy exactly, and the nucleation threshold
> is untouched. The advertised "new vertex selection rule" is a
> **tautology**. Both sections below now carry the measured refutation
> inline. §4 anticipated exactly this for Claim 4 and said to report the
> measured number; that is what happened.

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

### Claim 2 — a twisted link freezes mouth exchange **exactly** (the sharp new result) — **CONFIRMED, and stronger than claimed**

> **Measured on the real #224 network** (T4). Not the reduced model below —
> `build_two_throat` itself, whose ring operator closes with
> `lap[0,-1] = lap[-1,0]`, a genuine cycle:
>
> | `r_s` | `Δω` periodic | `Δω` twisted | ratio | `T_exchange` periodic → twisted |
> |---:|---:|---:|---:|---|
> | 0.25 | 9.464e-04 | 1.9e-13 | 2.1e-10 | 3319 → 1.6e+13 |
> | 0.30 | 1.714e-03 | 8.3e-14 | 4.9e-11 | 1833 → 3.8e+13 |
> | 0.40 | 4.747e-03 | 2.4e-13 | 5.2e-11 | 662 → 1.3e+13 |
> | 0.50 | 1.079e-02 | 2.4e-13 | 2.2e-11 | 291 → 1.3e+13 |
>
> **And it is a theorem, not a two-level accident.** The half-ring
> translation `S` (which swaps mouth A with mouth B) squares to the ring
> holonomy: `S² = W`. For `W = −1` the eigenvalues of `S` are `±i`, and
> since the ring operator is **real**, every level must pair into an
> exactly degenerate doublet. Measured across the low spectrum: max
> intra-pair gap **7.1e-12** twisted vs **2.1e-03** periodic. The freeze
> applies to the *whole spectrum*, not just the interior doublet — which
> the reduced two-level model could not have shown.

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

### Claim 3 — Z₂ gauge invariance is exact, and a twisted loop is an antiperiodic ring — **CONFIRMED; the Möbius sector was already in the repo**

> **Measured** (T3, T5). Gauge invariance on the real 1114-site #224 ring:
> moving the twist from the closing bond to bonds 0, 37, 400, 900 changes
> the spectrum by **0.00e+00** (it is a similarity transform by
> `diag(±1)`), while flipping the holonomy moves it by 1.7e-03. Only the
> Wilson loop is physical, exactly.
>
> **The identification is not an analogy — it is already implemented
> twice.** `geometrodynamics/qcd/topology.py` builds
> `make_glueball_ring` with `topology_kind='periodic'`, `orientation=+1`
> and `make_mobius_tube` with `topology_kind='mobius'`, `orientation=−1`.
> Those are the two holonomy sectors of one cycle. The ER-network link
> twist and the QCD Möbius label are **one Z₂ carried by two independent
> implementations that have never been connected.**

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

### Claim 4 — on the quotient, the antipodal revival period halves, and the twist sets its sign — **CONFIRMED; its threshold corollary REFUTED**

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

Both facts were **confirmed exactly** on the #166 solver (T6): the full
S³ tower reproduces #166's mirror refocus (1.7e-15), the untwisted sector
returns `−f₀` (0.00e+00) and the twisted sector `+f₀` (0.00e+00), all at
`t = πR`.

#### The threshold corollary — **REFUTED** (T7)

The plan conjectured that if the returning focus is coherent with the
source in the twisted sector, the energy density at the caustic goes as
the square of the sum — **a factor-of-4 energy advantage for the twisted
sector, and exact cancellation for the untwisted one**, a topological
selection rule on which sector can nucleate. It is wrong, on two counts:

  - **Energy is conserved exactly.** Measured ratio `E(πR)/E(0) =
    1.000000` in every sector. Linear evolution cannot amplify; a revival
    is a revival. The conjecture confused a *phase* relation with an
    *energy* one.
  - **The sectors do not differ.** Their antipodal caustic amplitude
    gains differ by **exactly `2/l_max`** (measured residual against that
    law: `0.0` at every cutoff — 1.050, 1.025, 1.0125, 1.00625 at
    `l_max = 40, 80, 160, 320`), a pure finite-cutoff artefact of the odd
    sector's mean `(l+1)` sitting one higher. It vanishes in the
    continuum, and at the physical cutoff `l_max ~ R/R_MID` it is
    negligible.

**The `2 m_e c²` threshold of #58/#166 is untouched by the twist.** What
survives from Claim 4 is real but smaller than conjectured: the
recurrence period halves and the returning sign is topological. One
genuine effect did show up and is worth recording — restricting to
*either* sector removes the `(−1)^l` alternation from the antipodal sum,
so a single-sector flat band focuses far more sharply than the full tower
(gain 2024 vs 8.9 at `l_max = 320`). That is a property of **the
quotient**, not of the twist.

### Claim 5 — the vertex rule — **the "new rule" is a TAUTOLOGY** (T8)

The plan advertised: #137/#138 derive `Σl` even; with legs carrying twist
classes `ε_i`, deck invariance requires `(−1)^{Σl} = ∏_i ε_i`, i.e.
`Σl ≡ n_twist (mod 2)`, so a vertex with an odd number of twisted legs
would select `Σl` **odd** — the complement of #137's rule. Pre-flight
enumeration "confirmed" it (32 of 64 triples for every `n_twist`).

**The pre-flight confirmed an identity.** A leg carrying twist `ε` is a
section of that flat bundle, and sections of the twisted bundle on RP³
are *exactly the odd-`l` harmonics* — so `ε_i = (−1)^{l_i}` is **forced,
not an independent label**, and `∏ε = (−1)^{Σl}` identically. The "new"
rule is #137's rule written twice. Verified against the real S³ triple
integral (`angular_triple`): every triple has `Σl mod 2 = n_twist mod 2`,
necessarily.

The lesson generalises: a Z₂ enumeration will always "pass" if the two
labels being compared are secretly the same label. What survives is
genuine but different:

  - **(a) A reinterpretation.** #137's `Σl`-even rule **is** the
    bundle-descent condition — the vertex face of the same classification
    Claim 1 gives the #129 boundary data. Not a new rule; a better
    account of the existing one.
  - **(b) A real *network* junction rule**, where the link twists genuinely
    are independent Z₂ gauge variables because they live on the ER graph,
    not on `l`. See below.

Corollaries (the second is the one that survives):

  - **Tri-mouth junctions (#208) — the rule that survives.** #208 derived
    a charged GHZ no-go by flux conservation (`k_A + k_B + k_C = 0` has
    zero solutions for a ±1 doublet) and concluded GHZ can live only in
    the transported-frame spin label. The twist adds a *second,
    independent* junction rule: **`ε_A ε_B ε_C = +1`**. Of the four
    inequivalent sign assignments to a Y-junction, **2 are admissible and
    2 are frustrated** (T8). Unlike Claim 5's vertex rule this is genuine,
    because the `ε_b` are link variables on the ER graph and carry no `l`.
  - **Bridge surgery (#207) — not computed.** The composition law
    `φ₁₄ = φ_a + φ_b + φ_c` should acquire a Z₂ term `+ π·[n_twist mod 2]`
    — surgery across an odd number of twisted bridges swapping the Bell
    outcome to its orthogonal partner. Machine-checkable against the
    existing 16-dim QM swapping calculation; **left open by this probe**.

## 3. The probe — as run (T1–T9, 8/8 PASS, ~1 min)

`experiments/closure_ledger/nonorientable_er_link_z2_probe.py`. Every test
runs on the **repo's real machinery**, not the reduced pre-flight models.

| # | test | uses | outcome |
|---|---|---|---|
| T1 | the object: link variables, gauge orbit, Wilson loops; the RP³-is-orientable precision note | — | stated |
| T2 | **the `b₁` audit — the primary falsifier** | graph cycle rank of #223/#224/#207/#208 | **does not fire** (see §4) |
| T3 | Z₂ gauge invariance | the real 1114-site #224 ring | drift **0.00e+00** across gauge copies; holonomy moves it 1.7e-03 |
| T4 | **twisted-link exchange freeze** | `build_two_throat` (#224) | ratio **2e-11 … 5e-11** over `r_s = 0.25–0.5`; degeneracy theorem `S²=W` verified (7.1e-12) |
| T5 | twisted loop = antiperiodic = the repo's Möbius sector | ring solver + `qcd/topology` | comb `0.5, 0.5, 1.5, …`; `make_glueball_ring` (+1) / `make_mobius_tube` (−1) confirmed |
| T6 | **quotient revival at `πR`, sign = twist class** | the #166 conformal zonal solver | `−f₀` (ε=+1) and `+f₀` (ε=−1), both **0.00e+00**; #166 control 1.7e-15 |
| T7 | caustic gain / threshold per sector | #166 caustic + #58 threshold | **CORRECTION**: energy ratio 1.000000; sector difference exactly `2/l_max` |
| T8 | vertex rule; junction rule | `cubic_vertex_ledger_probe.angular_triple` | **CORRECTION**: vertex rule tautological; junction rule 2 of 4 admissible |
| T9 | verdict + honest ledger | — | 8/8 |

**Verdict reached:**
`ER_LINK_TWIST_IS_A_Z2_GAUGE_FIELD_THAT_FREEZES_MOUTH_EXCHANGE_AND_HALVES_THE_QUOTIENT_REVIVAL_BUT_LEAVES_THE_NUCLEATION_THRESHOLD_UNTOUCHED`

— between the two classes anticipated below: the gauge structure and the
exchange freeze are established outright, one Arc A consequence (the
revival) holds exactly, and the other (the threshold) is refuted.

## 4. Honest scope — what held, what died

  - **Pre-flight ≠ probe, and the difference mattered.** Every §2 number
    was first fixed in a reduced model (two-level ring, tight-binding
    chain, zonal mode sums). Two of those pre-flights were misleading:
    the two-level model *understated* Claim 2 (the real effect is a
    spectrum-wide degeneracy theorem, not a doublet effect), and the
    enumeration for Claim 5 *overstated* it (it confirmed a tautology).
    A pre-flight can only show a claim is non-vacuous; it cannot tell you
    whether the labels you are comparing are independent.
  - **Claim 2 was the safest, and held** — with more than was claimed.
  - **Claim 4 was flagged as riskiest, and was trimmed exactly as
    predicted.** The plan said: "If the measured gain ratio is not 4,
    report the measured number; the revival identities stand regardless,
    and Claim 4 degrades to 'period halves, sign is topological'." That is
    precisely what happened. The threshold corollary is refuted; the
    revival identities are exact.
  - **The primary falsifier did not fire — but only on the right graph.**
    A Z₂ gauge field has content only if the ER graph carries independent
    loops (`b₁ > 0`). Measured (T2):

    | topology | `V` | `E` | `b₁` | reading |
    |---|---:|---:|---:|---|
    | #224 two-throat ring | 2 | 2 | **1** | as built |
    | #223 single bridge + exterior arc | 2 | 2 | **1** | as built |
    | #207 perfect matching (6 mouths) | 6 | 3 | **0** | bridges only — *the naive falsifier* |
    | #207 matching + shared S³ exterior | 6 | 9 | **4** | bridges + exterior — the physical graph |
    | #208 Y-junction, bare | 4 | 3 | **0** | a tree |
    | #208 Y-junction + shared S³ exterior | 4 | 6 | **3** | the physical graph |

    So the falsifier is real and would have fired on a bridges-only
    reading: #207's perfect matching **is** a forest, and there every
    `ε_b` gauges away. It does not fire because the physical graph
    includes the shared S³ exterior arcs as edges — which is not a
    convenient choice but exactly how #223/#224 build their ring, whose
    operator closes with `lap[0,-1] = lap[-1,0]`. **The twist is physical
    only where the network has cycles**, and stating which reading of the
    graph one is using is not optional.
  - **What this cannot do.** It does not derive the antipodal postulate
    (it classifies it), does not touch `λ₃, λ₄`, does not reduce any of
    the four standing residuals (`ε`, `n_part`, `α`, `√σ/m_e`), and adds
    no dimensionful input — the twist is an integer mod 2.
  - **`ℏ` enters nowhere**, and "probability" remains the classical
    normalized field-norm fraction, as everywhere in the program.

### Superseded text (kept for the record)

The pre-execution version of this plan asserted a factor-of-4 energy
advantage for the twisted sector at the caustic and a "new" vertex
selection rule `Σl ≡ n_twist (mod 2)` complementing #137. **Both are
withdrawn**, on the measurements in §2.4 and §2.5 respectively. The
original falsifier text follows, since it is what the audit was run
against:

  - *"The primary falsifier is topological, not numerical. A Z₂ gauge
    field has content only if the ER graph carries independent loops
    (`b₁ > 0`). #207's bridges pair mouths into a perfect matching —
    and a perfect matching is a forest, `b₁ = 0`, every `ε_b` gaugeable
    to +1. The twist is physical only where the network has cycles*:
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

## 7. Reproduce

```bash
# the probe (T1-T9, 8/8 PASS, ~1 min) -- writes probe.json + probe.md under
# experiments/closure_ledger/runs/<UTC timestamp>_nonorientable_er_link_z2_probe/
python -m experiments.closure_ledger.nonorientable_er_link_z2_probe

# the reduced pre-flight models the plan was written against (<5 s)
python -m experiments.closure_ledger.nonorientable_er_link_z2_preflight
```

Expected verdict:
`ER_LINK_TWIST_IS_A_Z2_GAUGE_FIELD_THAT_FREEZES_MOUTH_EXCHANGE_AND_HALVES_THE_QUOTIENT_REVIVAL_BUT_LEAVES_THE_NUCLEATION_THRESHOLD_UNTOUCHED`, 8/8 PASS.

## 8. What comes next

  1. ~~**The #207 surgery Z₂ term**~~ — **DONE, and the answer is no.**
     See §9.
  2. ~~**Is the Möbius identification quantitative?**~~ — **DONE (PR #229):
     partially.** The `+πσ` factorizes as `(1/2) × (2πσ)`. The `1/2` is
     derived, shared across three implementations, and topological
     (`O(1/M²)` exactly, circumference-independent to 1.4e-11) — #100 had
     only *asserted* it. The `2πσ` is an imported closed-string relation
     the network does not supply. And #100's common intercept across both
     towers is unjustified: antiperiodic moding shifts the zero point by
     `+π/(4C)` per polarization, so its Möbius ground state is not
     `M₀² + πσ`. Scope: the glueball tower only — #103/#109/#114's search
     table is built from `Δ = 2√σ` and stands as published. See
     `docs/mobius_identification_quantitative_research_plan.md`.
  3. ~~**Where does the twist come from?**~~ — **DONE (PR #230): energetics
     select, topology freezes.** The sectors differ by `π/(4C)` per degree
     of freedom, with the sign set by the **statistics of the link field**
     (bosonic ⟹ untwisted favoured; fermionic ⟹ twisted favoured). But the
     preference never drives a relaxation: changing `W` requires the link
     amplitude to pass through zero, severing the cycle (`b₁`: 1 → 0) —
     the #175 gate again — and `W` is exactly deformation-invariant
     otherwise. So the sector is chosen at nucleation and frozen. One
     Casimir sign gets **both** of the program's Möbius sectors right.
     See `docs/twist_sector_selection_research_plan.md`.

**With §8.1–§8.3 all closed, the #227 arc has no open roadmap items.**
What remains is named in the individual plans: the corrected Möbius `M²`
intercept (#229), whether the network can supply the `2πσ` spacing at all
(#229), whether any Bell-type observable is twist-sensitive (#228), and
absolute sector populations (#230).

## 9. §8.1 resolved — the surgery composition law has **no** Z₂ term

`experiments/closure_ledger/er_link_twist_surgery_probe.py` (T1–T5, 5/5
PASS, ~4 min) computes the conjecture on #207's own 4-mouth surgery
lattice, reusing its Hamiltonian, mouths, holonomies (`s₁₂ = s₃₄ = 2`)
and phase extraction unchanged. **The conjecture is refuted.**

Sweeping all 8 sign assignments `(ε₁₂, ε₃₄, ε_j)` at every junction
holonomy `s_j = 0..3`:

  - twisting **either outer bridge** moves the pair phase by at most
    **2.3e-13** — identically zero, not `π` — and leaves `|c_pm|`
    unchanged to 1e-14 as well;
  - twisting the **junction bridge** moves it by **4.68e-02 rad** —
    neither 0 nor `π`;
  - every odd-`n_twist` row misses the conjectured `π` by at least
    **3.09 rad**.

**The corrected law is the original law:** `φ₁₄ = φ_a + φ_b + φ_c`.
Surgery across twisted bridges does **not** swap the Bell outcome.

### Why — and why this is not a null result awaiting a better lattice

The reason is algebraic. #207 extracts the pair phase as
`arg(c_mp / c_pm)` — a **ratio of two branches that traverse the same
links**. A link-uniform sign multiplies both branches equally and cancels
identically; `|c_pm|` is a modulus, so it cancels there too. No choice of
lattice parameters can change this.

Crucially this is *not* because the twist is gaugeable away here — it
moves the spectrum by ~2.2e-02 to 2.8e-02, so the Z₂ is physical on this
lattice. It simply does not enter this observable.

### Where the twist does act — and it confirms T2 dynamically

The junction's 4.68e-02 residual is **bridge-vs-exterior interference
around the M2–M3 cycle**: those mouths sit only 8 lattice sites apart, so
the short exterior hop competes with the junction bridge and the relative
sign between the two paths is physical. It behaves like an interference
amplitude, not a topological phase:

| M2–M3 gap | `Δφ` junction twist | `Δφ` outer twist |
|---:|---:|---:|
| 4 | 2.877e-01 | 1.221e-13 |
| 8 | 4.678e-02 | 2.283e-13 |
| 16 | 2.139e-04 | 2.984e-13 |
| 32 | **2.975e-14** | 1.776e-13 |
| 48 | **1.763e-13** | 1.337e-13 |

| exterior hopping `t_x` | `Δφ` junction twist |
|---:|---:|
| 0.25 | 2.486e-05 |
| 0.50 | 1.154e-03 |
| 1.00 | 4.678e-02 |
| 2.00 | 1.770e-01 |

It **dies with the cycle length** (machine zero once the competing path
is long enough) and **scales continuously with the exterior hopping that
carries the competing path**. A topological phase does neither. The
outer-bridge twists, whose competing exterior path is 60 sites long, stay
at ~1e-13 for every gap and every hopping — they never do anything.

So the twist acts **only through cycles**, exactly as T2's `b₁` audit
requires — now confirmed *dynamically* rather than combinatorially. In
#207's geometry that cycle carries a tunable interference amplitude, not
a Bell-outcome flip.

### What this leaves open

Whether *any* Bell-type observable is twist-sensitive. It would have to
compare two paths traversing **different link sets** — which #207's
geometry does not naturally provide, since both branches of its pair
phase run through the same bridges. That is a constraint on what to look
for, not a gap in this calculation.

```bash
python -m experiments.closure_ledger.er_link_twist_surgery_probe
```

Expected verdict:
`THE_SURGERY_COMPOSITION_LAW_HAS_NO_Z2_TERM_THE_PAIR_PHASE_IS_A_RATIO_AND_THE_TWIST_ACTS_ONLY_AS_CYCLE_INTERFERENCE`, 5/5 PASS.

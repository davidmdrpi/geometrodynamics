# The RP³ spin-structure selection arc: Probes A–E (PR #230)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## 0. Why this arc exists

The twist-selection investigation had reached an honest impasse: two
candidate explanations for why matter sits in the twisted sector had been
tried and both failed — the statistics sign (wrong moding) and the
zero-mode index (not an index). Both failures traced to the same root
cause: the constructions were **1D reductions using the SUSY square root**,
which delivers `|spec(D)|` only, discarding the sign of the spectrum.

This arc replaces those constructions with the genuine objects.

| probe | file | question |
|---|---|---|
| **A** | `rp3_spin_lift_probe.py` | what are the lifts, and are B2 and `W` independent? |
| **B** | `rp3_dirac_eta_probe.py` | the genuine first-order RP³ Dirac spectrum, `η`, `h`, `det` |
| **C** | `bulk_aps_spin_structure_probe.py` | does inflow select? does the Pin⁻ mouth? |
| **D** | `twist_sector_einstein_dirac_probe.py` | does back-reaction select? |
| **E** | `network_cycle_pi1_class_probe.py` | does the #224 cycle generate `π₁(RP³)`? (decides whether B–D apply) |
| **F** | `rp3_spin_structure_on_223_ring_probe.py` | (PR #231) push the machinery through the **#223** ring, where E showed the labels *do* coincide |

## 1. Probe A — the spin-lift ledger (6/6)

`J: S³ → S³`, `x ↦ −x`, is `−I₄` with `det = +1`: it lies in **SO(4)**, which
is why RP³ is orientable (the repo's #183 deck determinant, recovered).

In `Spin(4) = SU(2) × SU(2)` acting by `x ↦ q_L x q_R†`, solving
`q_L x q_R† = −x` forces exactly two lifts:

| lift | projects to | `J̃²` |
|---|---|---|
| `J̃₊ = (−1, 1)` | `−I₄` | `(1,1) = **+1**` |
| `J̃₋ = (1, −1)` | `−I₄` | `(1,1) = **+1**` |

differing by the covering kernel `(−1,−1)`.

**Both square to +1.** This refutes a common shortcut: the two spin
structures are *not* distinguished by `J̃² = ±1` — they are distinguished by
**which SU(2) factor carries the −1**, a distinction invisible to any
argument phrased in terms of "the square of the lift", including the
`η = ±1` bookkeeping the earlier work used.

**Where BAM sits.** `σ(z₁,z₂) = (z̄₂, −z̄₁)` has `det = +1` — it is
**orientation-preserving** on S³ (it reverses the Hopf *fibre*, not S³) —
and `σ² = −I₄ = J` exactly: σ is a **square root of the antipodal map**. Its
spinor lift is BAM's `T = iσ_y`, with `T² = −I` in the *left* SU(2), so
`T² = (−1,1) = J̃₊`. **B2 selects a specific lift.**

### 1.1 The consequential finding: B2 and `W` are independent

| datum | lives in | detected by |
|---|---|---|
| B2 spin structure | `H¹(RP³;Z₂) = Z₂` | the deck loop generating `π₁(RP³)` |
| network holonomy `W` | `H¹(Γ;Z₂)`, Γ the ER graph | Wilson loops on graph cycles |

Different cohomology of different spaces. They coincide **only if** the
network cycle carrying `W` is itself the `π₁(RP³)` generator — a statement
about how the cycle embeds, not an identity.

**Consequence:** the composition rule `θ = 0 iff W·η = +1`, on which the
moding swap of `bam_spinor_spectrum_effective_action_probe` rests, is **not
justified in general**. It is conditional on the network cycle generating
`π₁(RP³)`, which this arc does **not** settle.

## 2. Probe B — the genuine RP³ Dirac spectrum (8/8)

Analytic input, no finite differencing anywhere (so no fermion doubling to
manage): on the round unit S³, `λ = ±(3/2 + k)` with multiplicity
`(k+1)(k+2)`. Writing S³ = SU(2), the antipodal map is left multiplication
by the central `−1`, acting on left-spin `j_L` by `(−1)^{2j_L}`. The two
branches sit at **different** left-spin:

```
+ branch, level k :  j_L = k/2      →  parity (−1)^k
− branch, level k :  j_L = (k+1)/2  →  parity (−1)^{k+1}
```

**Opposite parity at the same `|λ|`** — the entire origin of the spectral
asymmetry, and exactly what a `|spec|` construction cannot see.

*Consistency check built in:* S³ has `2(k+1)(k+2)` modes at each level and
the projection keeps exactly `(k+1)(k+2)` — precisely half, as
`vol(RP³) = vol(S³)/2` requires, for both lifts at every level.

### The invariants

| spin structure | `η(0)` | `h` | `\|det D\|` | `arg det D` |
|---|---:|---:|---:|---:|
| `ε = +1` | **+1/4** | **0** | 0.8963003229 | −π/8 |
| `ε = −1` | **−1/4** | **0** | 0.8963003229 | +π/8 |

`ζ_{D²}(0) = 0`, `ζ'_{D²}(0) = 0.2189594807`. The `η` values are exact
(analytic continuation in alternating Hurwitz zetas), cross-checked against
the directly convergent mode sum at `s = 4, 6` to **2e-19**.

Two results matter:

  - **`h = 0` for both.** Every eigenvalue has `|λ| ≥ 3/2 > 0`: the genuine
    3D operator has **no zero mode at all**. This settles the earlier
    probe's kernel from the other side — its 1D "zero mode" was an artefact
    of the reduction, confirming that retraction independently.
  - **The moduli are identical.** For each `k` exactly one sign survives,
    with the same `|λ|` and multiplicity, so `ζ_{D²}`, `|det D|` and the
    vacuum energy **cannot distinguish the sectors**. The two differ *only*
    in the determinant phase, `π/4` apart.

**This retires the energetic route exactly.** The `π/(4C)` energetics of
the twist-selection arc were computing a quantity whose difference is
identically zero on the genuine 3D operator.

The instinct that the answer was index-theoretic was right; the object was
wrong. The discriminator is a spectral **asymmetry** (`η`), not a kernel
and not an index.

## 3. Probe C — bulk extension and inflow (5/5)

**The objects, concretely.** `Y³ = RP³ = L(2,1)`. `X⁴` = the disk bundle
over S² of Euler number `e = ±2`, whose boundary is `L(2,1)` and which is
spin precisely because `e` is even.

**APS.** `ind(D_X) = ∫_X Â − (h+η)/2` with `h = 0`, `η = ±1/4`, giving a
fractional boundary defect `∓1/8`, opposite for the two spin structures.

**Both extend.** `Ω₃^Spin = 0`, so every spin 3-manifold bounds; both spin
structures are null-bordant, the index is integral for each, and the
fractional defect is absorbed by the bulk `Â` integral. **Anomaly inflow
does not select** — and that is precisely the precondition triggering
Probe D.

### 3.1 Where a real selection rule does appear

The mouth is RP², carrying Pin⁻. `Ω₂^{Pin⁻} = Z₈` generated by RP², with
`ABK(RP²) = ±1 mod 8`. A throat has two mouths and the classes add:

| mouth A | mouth B | total mod 8 | bounds? |
|---:|---:|---:|---|
| +1 | −1 | 0 | **yes** |
| −1 | +1 | 0 | **yes** |
| +1 | +1 | 2 | no |
| −1 | −1 | 6 | no |

**The two mouths of a throat must carry opposite Pin⁻ structures** — which
the repo already asserts as a modelling convention (`ThroatDefect` /
`ConjugatePair` enforce `mouth_a.orientation == −mouth_b.orientation`).
Here it is *derived* from Z₈ bordism.

**Scope, plainly:** this is a **relative** constraint between the two
mouths. It says nothing about which absolute sector a network occupies, so
it does not answer the selection question either.

## 4. Probe D — semiclassical back-reaction (6/6)

Precondition met (Probe C). The result follows from Probe B without further
computation:

| quantity | differs? | why |
|---|---|---|
| radial semiclassical equations | no | identical `⟨T_AB⟩` source |
| regularity | no | same equations, same boundary data |
| on-shell action **modulus** | no | `\|det D\|` identical |
| on-shell action **phase** | **yes** | `arg det D = ∓π/8` |
| negative modes / stability | no | same real fluctuation operator |

  - `ΔE_vac = 0` **exactly** — the vacuum energy depends only on `|spec|`,
    and the two agree mode by mode.
  - `Δ⟨T_AB⟩_ren = 0` **pointwise** — RP³ is homogeneous, so the vacuum
    stress tensor is fixed up to a normalization set by the total energy;
    equal energy and equal volume force equality locally, not just on
    integration.

**A phase does not source the semiclassical Einstein equations** (`⟨T_AB⟩`
is built from the modulus), so the one surviving discriminator is invisible
to gravity at this order.

**What is not done.** The degeneracy argument leans on homogeneity of the
round RP³. BAM's actual exterior is the Tangherlini throat geometry, which
is not homogeneous, and there equal integrated energy does not force equal
pointwise `⟨T_AB⟩`. A genuine solve would need the mode functions on the
real radial background, the renormalized coincidence limit of the Green
function *difference*, and then the coupled radial system. None is
attempted. What is established is the **baseline**: any back-reaction
selection must come entirely from the inhomogeneity, as a correction on top
of an exact zero.

## 4a. Probe E — the #224 cycle does **not** generate π₁(RP³) (6/6)

`network_cycle_pi1_class_probe.py`. A loop in `RP³ = S³/{±1}` is classified
by lifting to `S³`: the lift either closes (class 0) or ends at the
**antipode** (class 1, the generator). For geodesic exterior arcs on the
unit sphere that is decided by total arc length mod 2π.

Read from the constructions themselves, not the prose:

| network | arcs | total exterior | `S³` lift | deck class | generator? |
|---|---:|---:|---|---:|---|
| **#224** `build_two_throat` | 2 | 1.999π | **closes** | **0** | no |
| **#223** single-bridge ring | 1 | 1.000π | antipodal | 1 | **YES** |

`#224` sets `Ltot = 2*_D_CH + 4*_SIG_M + 2*_LARC` with `_LARC = 3.14` —
**two** π-arcs, total 2π, so the lift closes. `#223` has **one** π-arc to
the antipodal source point.

**And the ambient space is not RP³ anyway.** A wormhole is a 1-handle, and
attaching `n` of them gives `π₁ = π₁(exterior) * Z * … * Z`, one free `Z`
per handle. The #224 cycle passes through both throats, so its class lives
in the **free (handle) part**, while its deck component — the only part B2
grades — is trivial. `W` on that cycle is a character on a free `Z`
generator; B2 is a character on the `Z₂` deck generator. Different
generators of different groups.

*Robust to the reading:* exterior as `S³` gives `π₁ = 0` so the deck
component is trivially 0; exterior as `RP³` gives `generator² = 0`. Both 0.

### What this does to Probes B–D

| probe | labelled by | bears on the #224 twist? |
|---|---|---|
| B (`η = ±1/4`, `h`, `det`) | the deck `Z₂` | **no** |
| C (inflow, ABK) | the deck `Z₂` (ABK: the mouth Pin⁻ structure) | **no** |
| D (back-reaction) | the deck `Z₂` | **no** |

**This is a scoping correction, not a refutation.** Nothing in B–D is
withdrawn — every number stands — but their target narrows: they answer
"how do the two RP³ **spin structures** differ?", and the #224 network
twist is not that label.

Two consequences follow, and both matter:

  1. The composition rule `θ = 0 iff W·η = +1` is justified on the **#223**
     ring and **not** on the **#224** ring — and #224 is exactly where the
     arc's headline exchange-freeze result lives.
  2. The "three selection mechanisms all fail" conclusion of §5 is a
     statement about the **spin-structure** sectors. **The network twist
     sector remains untested by any of them.**

> **FOLLOWED UP (PR #231, Probe F).** The machinery has since been pushed
> through the **#223** ring, the one network E found with deck class 1.
> There the `end` (source-point) parity **is** `W`, and the predicted
> moding shift is **measured**: the two `W` sectors interleave with
> fractional offset → **0.5**. So the composition rule is confirmed on an
> actual construction. And with the labels genuinely identified, B–D apply
> to `W` — and **still nothing selects**. The by-product is sharper than
> the result: #223 couples the labels but has **no mouth doublet**, while
> #224 has the doublet but decouples them — **the phenomenon and the
> coupling never co-occur**, so RP³ spin-structure data cannot select the
> exchange-freeze sector even in principle. See
> `docs/rp3_spin_structure_on_223_ring_research_plan.md`.

## 5. Where the selection question now stands

**Three independent selection mechanisms have been tested and none
distinguishes the SPIN-STRUCTURE sectors** (and per §4a, none of them bears
on the *network twist* sector at all):

| mechanism | verdict |
|---|---|
| energetics / vacuum energy | **exactly degenerate** (identical `\|spec\|`) |
| anomaly inflow / APS | **both extend** (`Ω₃^Spin = 0`) |
| semiclassical back-reaction | **exactly degenerate** (on the round exterior) |

The only surviving distinction in the entire arc is **`η = ±1/4`**, a
determinant phase. It is a genuine, exact, index-theoretic invariant — and
it selects nothing energetically.

Two candidate resolutions remain open, and they are different in kind:

  1. **The phase is physical but not selective** — it would show up in
     CP-type observables or interference between sectors, not in which
     sector is occupied. Then "why matter is twisted" is not an energetic
     question at all and needs a different frame (e.g. a nucleation-history
     argument, per the topological freeze).
  2. **The relevant Z₂ is not the one being computed** — Probe A's
     independence result, now **confirmed concretely by Probe E**: the
     network `W` and the RP³ spin structure ARE different data on the #224
     ring. This arc constrains B2 and leaves `W` untouched.
     *This is the branch the evidence actually supports.*

Settling (2) required determining whether the #224 network cycle generates
`π₁(RP³)`. **That computation has now been done — Probe E — and the answer
is no.** See §4a.

## 6. What this arc overturns in the earlier probes

| earlier claim | status |
|---|---|
| the twisted sector carries a Dirac zero mode | **artefact** — `h = 0` on the genuine operator |
| energetics favour untwisted by `π/(4C)` | **not a discriminator** — difference is exactly 0 in 3D |
| the moding swap from `W·η` | **conditional** — the composition rule is unjustified (Probe A) |
| paired mouths carry opposite orientation | **derived** — was an assertion, now forced by ABK/Z₈ |
| the spin-structure results apply to the #224 twist | **scoped out** — Probe E: the #224 cycle has trivial deck class, so B–D target a different label |

## 7. Reproduce

```bash
python -m experiments.closure_ledger.rp3_spin_lift_probe             # A, 6/6
python -m experiments.closure_ledger.rp3_dirac_eta_probe             # B, 8/8
python -m experiments.closure_ledger.bulk_aps_spin_structure_probe   # C, 5/5
python -m experiments.closure_ledger.twist_sector_einstein_dirac_probe  # D, 6/6
python -m experiments.closure_ledger.network_cycle_pi1_class_probe     # E, 6/6
```

## 8. Cited inputs (not derived here)

  - the S³ Dirac spectrum `±(3/2+k)` with multiplicity `(k+1)(k+2)`
  - `Ω₃^Spin = 0`; `Ω₂^{Pin⁻} = Z₈`; `ABK(RP²) = ±1 mod 8`
  - the disk bundle over S² is spin iff its Euler number is even

## 9. Cross-references

  - `docs/bam_spinor_spectrum_effective_action_research_plan.md` — the 1D
    construction this arc supersedes
  - `docs/twist_sector_selection_research_plan.md` — the energetics this
    arc retires
  - `docs/k1_zero_mode_index_mechanism_research_plan.md` — #195, genuine
    index-protected zero modes (S² + monopole), for contrast
  - `docs/pin_rp2_fermi_statistics_research_plan.md` — #170, Pin⁻ on the mouth
  - `geometrodynamics/embedding/topology.py` — `ThroatDefect`, whose
    opposite-orientation assertion Probe C derives

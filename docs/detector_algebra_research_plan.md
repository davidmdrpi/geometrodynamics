# Probe K — the detector algebra and the marginals (PR #238)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## 0. The question

#237 established that **correlator coverage is not the remaining
problem**: once convex mixing is allowed, the hull of BAM's tables is the
whole quantum body. What remains is the **detector marginals** and the
**measurement algebra**. This probe makes both precise — and they turn out
to be one question, with an uncomfortable answer.

## 1. The generated algebra, computed separately per pairing

| setting × detector | accessible (real dim) | generated algebra (complex dim) | abelian? | `= span{I, σ_z}`? |
|---|---:|---:|---|---|
| fiber U(1) × `σ_z` — **both derived** | **1** | **2** | **yes** | **yes** |
| fiber U(1) × `σ_x` *(assumed)* | 2 | 4 | no | no |
| `_rot` y-rotation × `σ_z` *(the code)* | 2 | 4 | no | no |

> **For the fully derived pairing the generated `*`-algebra is exactly
> `span_C{I, σ_z}` — complex dimension 2, and abelian.** The accessible
> observable set is a single point, because a fiber rotation cannot move a
> winding-diagonal detector.

That is the real explanation of the `CHSH = 2.000000` in §3, and it is a
**theorem rather than an observation**: a commutative algebra of
observables admits a joint distribution over all of them at once, hence a
local hidden-variable model, hence `CHSH ≤ 2` necessarily. The other two
pairings generate all of `M₂(C)`, which is what makes a violation possible
at all.

> **Correction to this probe's own first answer.** The first version
> computed the generated algebra only for the `fiber × σ_x` pairing, found
> dimension 4, and concluded *"the algebra is not the restriction — the
> dialable set is."* For the pairing that matters the algebra **is** the
> restriction, and it is abelian. Computing it per pairing, rather than
> once and generalising, is what shows this.

## 2. The docs and the code disagree about what a setting is

The documents describe the setting as *"the fiber-frame rotation before
the device"*. A rotation of the `χ` fiber multiplies the winding-`k`
channel by `e^{2πikδ/N}`, so on the `k = ±1` qubit it is

```
fiber:  diag(e^{iθ}, e^{−iθ})        — a z-rotation
```

But the committed `measurement_sector_probe._rot` is `[[c, −s], [s, c]]`,
which is

```
code:   exp(−iθ σ_y / 2)             — a y-rotation
```

verified identical to **0.0e+00** — with `_rot` **imported directly from
`measurement_sector_probe`**, not reimplemented, so the comparison is
against the code the repo actually runs. These are not the same operation,
and the difference is exactly the one that matters:

| | `‖[·, σ_z]‖` |
|---|---:|
| fiber rotation | **0.00e+00** (commutes exactly, every angle) |
| y-rotation | up to 1.95 |

**A fiber rotation cannot move a winding-diagonal detector.** Under the
setting the docs describe, the dial does nothing.

*Method check.* The committed code applies the setting by rotating the
**state** (`np.kron(_rot(a), _rot(b)) @ _SINGLET`), while this probe
conjugates the **observable**. The two are verified to give identical
correlators (max difference **3.3e-16**), so nothing in the comparison
turns on that choice.

## 3. The trilemma

Maximum CHSH on the committed #206 singlet, by (setting, detector):

| setting | detector | both derived? | dialable span | max CHSH |
|---|---|---|---:|---:|
| fiber U(1) *(derived)* | `σ_z` winding SG *(derived)* | **yes** | **1** | **2.000000** |
| fiber U(1) *(derived)* | `σ_x` transverse *(assumed by #236/#237)* | no | 2 | 2.828427 |
| y-rotation *(in #209 code)* | `σ_z` winding SG *(derived)* | no | 2 | 2.828427 |

> **The only pairing in which both the setting and the detector are
> physically derived by this repository gives exactly the local bound.**
> No Bell violation at all. The dialable span there is **1** — a single
> point — because the setting commutes with the detector.

The violation requires **exactly one** non-derived ingredient: a
transverse detector, or a winding-mixing setting. The repo supplies
neither. And §1 says why in one word: both of those are ways of making the
generated algebra **non-abelian**, which a violation requires and the
derived pairing does not have.

## 4. The marginals follow the same fork — and #237 is corrected

| `β` | `\|cos 2β\|` | fiber×`σ_z` | fiber×`σ_x` | y-rot×`σ_z` |
|---:|---:|---:|---:|---:|
| π/12 | 0.8660 | 0.8660 | 0.0e+00 | 0.8660 |
| π/8 | 0.7071 | 0.7071 | 0.0e+00 | 0.7071 |
| π/6 | 0.5000 | 0.5000 | 0.0e+00 | 0.5000 |
| π/4 | 0.0000 | 0.0000 | 0.0e+00 | 0.0000 |

#237 reported *"identically zero marginals, immune to convex mixing"* as
its one surviving falsifiable content. That measurement was made with the
`fiber × σ_x` pairing, and is correct **there**. Under the other two
pairings the marginals are `|cos 2β|` exactly — up to **0.866**.

> **So the zero-marginal restriction is an artifact of one of three
> modelling choices and is retracted as a general claim about BAM.** It
> holds only for the transverse-detector model. #237's falsifiable content
> was contingent on the same undetermined choice that §3 shows the Bell
> violation itself hangs on. **The marginals and the algebra are one
> question, not two** — which is what the review that prompted this probe
> said.

## 5. Both missing ingredients are the same requirement

`σ_x` in the winding basis is purely **off-diagonal** between `k = +1` and
`k = −1` — a coherence between distinct topological charges — and the
y-rotation setting mixes those same sectors. So the one thing BAM must
supply is:

> **an operation or observable at a mouth that coherently connects
> distinct winding sectors.**

And if winding is superselected, that is unavailable. Restricting both
parties to `k`-diagonal dichotomic observables gives, over 4000 random
draws,

```
max CHSH = 2.000000        exactly the local bound
```

The repo already carries a superselection structure in this sector —
#208's charged-GHZ zero-sum no-go — and nothing in it has been reconciled
with the Bell chain.

## 6. What this does and does not say

It does **not** show BAM's Bell violation is wrong. It shows that the
violation rests on one ingredient the program has not derived, that the
two candidate ingredients are the *same* physical requirement, and that
the requirement is in tension with a superselection structure the repo
committed elsewhere. Naming which single thing has to be supplied is the
deliverable.

| claim | status |
|---|---|
| #237: "identically zero marginals … the whole of the probe's falsifiable content" | **retracted as a general claim** — true only for `fiber × σ_x` |
| this probe's own first answer: "the algebra is full, so the algebra is not the restriction" | **corrected** — computed for one pairing only; the derived pairing generates the abelian `span{I, σ_z}` |
| BAM's Bell violation (#206/#209) | **not refuted, but located** — needs exactly one non-derived ingredient |
| the docs and the code agree on what a setting is | **false** — fiber rotation vs `exp(−iθσ_y/2)` |

## 7. Tests

| # | test | outcome |
|---|---|---|
| T1 | goal | posed |
| T2 | **generated algebra per pairing** | derived → `span{I,σ_z}`, dim 2, abelian; others dim 4 |
| T3 | `_rot` imported and tested directly | = `exp(−iθσ_y/2)` to 0.0e+00; state-rot ≡ obs-conj to 3.3e-16 |
| T4 | **the trilemma** | fully derived pairing → **2.000000** |
| T5 | marginals follow the same fork | `\|cos 2β\|` vs 0; #237 corrected |
| T6 | both ingredients = winding coherence | superselected CHSH = 2.000000 |
| T7 | consequences | 1 retracted, 1 rediagnosed, 1 located, 1 false |
| T8 | assessment | 8/8 |

## 8. Open

  - Whether **any** mouth operation or observable in BAM coherently
    connects distinct winding sectors. Supplying that settles the
    violation and the marginals at once.
  - How #208's charged-GHZ superselection structure is to be reconciled
    with a Bell chain requiring exactly the coherence it removes.
  - Whether relocating the Bell chain onto the **transported-frame / spin
    doublet** (#192/#197 Berger-squash, named in #209 as the analogous
    device for that carrier) avoids the problem. Untested here, and the
    natural next probe.
  - This is the two-input two-output **winding** sector on one lattice;
    the multipartite chain (#207/#208) is not analysed.

## 9. Reproduce

```bash
python -m experiments.closure_ledger.detector_algebra_probe
```

Expected verdict:
`THE_DERIVED_PAIRING_GENERATES_ONLY_THE_ABELIAN_ALGEBRA_SPAN_I_SIGMA_Z_SO_ITS_CHSH_IS_EXACTLY_TWO_BY_COMMUTATIVITY_AND_BOTH_VIOLATING_PAIRINGS_NEED_A_NONABELIAN_ALGEBRA_REQUIRING_WINDING_SECTOR_COHERENCE`, 8/8 PASS.

## 10. Cross-references

  - `docs/admissible_correlation_tables_research_plan.md` — #237, whose
    zero-marginal claim this corrects
  - `docs/tsirelson_ceiling_research_plan.md` — #236, whose transverse
    detector this identifies as an assumption
  - `docs/measurement_sector.md` — #209, the winding Stern–Gerlach and the
    CHSH run whose setting operation this examines
  - `docs/configuration_space_emergence.md` — #206, the pair state

# Probe J — which correlation tables does BAM's encoding admit? (PR #237)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## 0. The question

#236 located BAM's Tsirelson ceiling at one algebraic fact, `B² = I`, by
probing the CHSH **maximum**. A maximum is one number. The encoding
represents a whole correlation table

```
E(x, y) = ⟨ψ| A_x ⊗ B_y |ψ⟩ ,     A_x, B_y dichotomic
```

so the structural question is not how large CHSH can get but **which
tables the representation admits at all** — and, separately, which tables
the committed implementation can actually produce. Those turn out to be
very different sets, and the gap is the result.

## 1. Level A — the representation in the abstract

It admits **exactly the unit-vector Gram matrices** `E_xy = ⟨u_x, v_y⟩`,
equivalently the Tsirelson–Landau–Masanes body

```
max over the four sign variants of | Σ ± arcsin E_xy |  ≤  π
```

Checked in **both** directions rather than cited:

  - **reverse** — an independent unit-vector fit agrees with the TLM
    predicate on **60/60** random tables;
  - **forward** — 4000 tables built *from* random unit vectors never
    violate TLM (minimum slack **+2.6e-04**).

| body | fraction of the correlation cube |
|---|---:|
| local polytope | 66.83% |
| quantum (TLM) | 92.64% |
| no-signaling | 100% |

The PR box sits at TLM slack **−π**: not representable, as it must not be.

## 2. Level B — what BAM actually implements is strictly smaller

Read off the committed lattice rather than assumed: the pair state's
T-matrix has x–y block **exactly `c·I₂`**, with `c = −sin 2β`.

| `β` | `c` | `−sin 2β` | deviation from `c·I₂` |
|---:|---:|---:|---:|
| π/12 | −0.50000 | −0.50000 | 2.1e-15 |
| π/8 | −0.70711 | −0.70711 | 3.1e-15 |
| π/6 | −0.86603 | −0.86603 | 3.0e-15 |
| π/4 | −1.00000 | −1.00000 | 2.8e-15 |

Combined with #236's finding that the settings are a **U(1)** confined to
the x–y plane, BAM's achievable tables are **exactly**

```
E_xy = c · cos(θ_x − φ_y) ,     c = −sin 2β ∈ [−1, 1]
```

That parametrization, not a bound, is the answer to "which tables" for the
implementation.

## 3. Coverage — and a methodological trap worth recording

| | coverage of the quantum body |
|---|---:|
| maximal preparation `c = ±1` (`d = 2` Gram) | **0.0%** — a measure-zero 3-parameter surface |
| with the bridge preparation `β` free | **55.2% ± 2.0%** (2000-sample run: 54.6% ± 1.1%) |
| general `d = 3` | 100.0% |

**The bridge preparation is the knob that buys the missing dimension.**
Without it the implemented family is a 3-parameter surface in a
4-dimensional body — measure zero. With `β` free it becomes
full-dimensional and covers a little over half. `β` is not decorative.

Of the `c = 1` surface, **33.7%** lies exactly *on* the Tsirelson boundary
— which is why #236 found saturation there. Saturating the boundary and
filling the body are different things.

> **The trap.** Membership was first tested with generic optimizers, which
> got it wrong **in both directions**: L-BFGS-B reported 48.3% and
> multi-restart Nelder-Mead 63.3%. *A failed fit is not proof of
> non-representability.* The probe now uses an **exact** test: fixing the
> global rotation with `θ₀ = 0` leaves four unknowns for four equations,
> three of which solve in closed form, leaving `E₁₁` as a single
> consistency condition in `c` alone — a 1-D root find on each of 8 sign
> branches, bracketed and bisected. It is validated by recovering
> **100%** of tables *built* from the family (max residual 3.8e-14), which
> neither optimizer does, and by the check that nothing it admits lies
> outside the quantum body (min admitted TLM slack **+9.7e-03**).

**A table unreachable from a *single* preparation:**

```
[[ 0.7424, -0.9971],
 [-0.5465,  0.4867]]
```

TLM slack **+0.740** (comfortably quantum) · `d = 3` Gram residual
**4.3e-09** · BAM single-shot residual **4.6e-01**. §3a decomposes it.

## 3a. But convex mixing dissolves that gap entirely

Shared randomness over preparations is an ordinary operational resource,
and it makes the reachable set the **convex hull** of the `c = 1` cosine
tables. Tested by linear programming — `min ‖Gᵀλ − E‖₁` over the simplex.

| generators | coverage (ladder set) | max residual |
|---:|---:|---:|
| 2 000 | 99.2% | 1.2e-02 |
| 8 000 | **100.0%** | 0.0e+00 |
| 32 000 | **100.0%** | 0.0e+00 |

On **300** fresh random TLM tables: **98.7%**, with support size at most
**5** — exactly the Carathéodory bound for a 4-dimensional body. The
residual failures are boundary resolution, not a gap: refining to 128 000
generators *on those tables only* takes 4 failures → **1**, worst residual
**7.2e-03 → 4.2e-03**. The PR box stays outside (residual **1.18**).

**And the witness above decomposes into exactly 5 `c = 1` cosine tables**,
residual **0.0e+00**, reconstruction error **4.4e-16**, weights
`[0.2124, 0.3447, 0.0409, 0.0102, 0.3918]`.

> **So this PR's own first framing was wrong.** "~45% of quantum
> correlation tables are unreachable" is a **single-shot** statement
> dressed as an operational one. With mixing there is no
> correlation-table restriction to report at all. The hull is contained
> in the quantum body by convexity and covers it to sampling resolution,
> so the two coincide.

## 4. The one restriction that survives

A correlation table is only part of a behavior. Under BAM's own U(1)
settings the marginal is

```
max |⟨A_θ⟩| = 0.000e+00       across the sweep
```

— identically zero, so **no behavior with biased marginals is
representable at all**. This is not because the state lacks the structure:
the same states carry a z-marginal of `cos 2β` (0.866 at `β = π/12`,
matching to 1e-6) which the U(1) settings simply cannot access, because
they never leave the x–y plane. **The restriction is in the readout, not
the state — the same place #236 found the ceiling.**

And unlike the correlation gap, **this one survives convex mixing**: a
mixture of zero-marginal behaviours has marginals `Σᵢ λᵢ·0 = 0`,
identically. §3a dissolves one of the two restrictions and leaves this one
completely intact, which makes it the whole of the probe's falsifiable
content.

## 5. What this does to #236 — scope, not retraction

| claim | status |
|---|---|
| #236's T3: the U(1) is "exactly sufficient — neither a sub-quantum deficit nor any excess" | **scoped** — correct for the CHSH *maximum*, which is what it measured. For the *set of tables* the U(1) alone reaches a measure-zero surface, and even with `β` free reaches ~55% |
| the bridge preparation `β` is a state parameter | **promoted** — it is the one knob lifting the family from measure zero to full dimension |
| BAM can represent quantum correlations | **qualified, then largely restored** — ~55% single-shot, but **100%** under mixing; 0% of biased-marginal behaviours either way |
| this PR's own first framing: "~45% unreachable" | **corrected** — a single-shot statement, dissolved by §3a |

Extremal saturation and full representational power are different claims,
and only the first was established. The phrase "exactly sufficient"
invites over-reading and is scoped here.

## 6. Falsifiable content

One restriction survives and one does not. The correlation-table gap is
dissolved by mixing (§3a). The **zero-marginal** restriction is immune to
it, so BAM's implemented encoding cannot reproduce **any** quantum
behaviour with biased marginals. That is experimentally ordinary, so it is
a prediction that can fail rather than a safe one — and it is the whole of
what this probe claims.

Whether the restriction is fundamental or an artifact of the committed
readout is the open question. It would be lifted by any physical operation
at a mouth that rotates out of the x–y plane, and nothing here says
whether the throat geometry supplies one.

## 7. Tests

| # | test | outcome |
|---|---|---|
| T1 | goal | posed |
| T2 | Gram ⟺ TLM, both directions | 60/60; min slack +2.6e-04 |
| T3 | the three bodies | 66.8% / 92.6% / 100%; PR box at −π |
| T4 | the implemented family, from the lattice | x–y block `c·I₂` to 3.1e-15, `c = −sin 2β` |
| T5 | **single-shot** coverage | 0.0% at `c = ±1`, 55.2% with `β` free, 100% at `d = 3` |
| T6 | **convex hull by LP** | 100% coverage; witness = 5 generators; PR box excluded |
| T7 | marginals vanish — and survive mixing | 0.000e+00 |
| T8 | consequences; #236 and #237's own headline scoped | 1 scoped, 1 corrected, 1 promoted, 1 surviving |
| T9 | assessment | 9/9 |

## 8. Open

  - Whether the throat geometry supplies **any** physical mouth operation
    that rotates out of the x–y plane — that, and only that, would lift
    the surviving zero-marginal restriction.
  - Whether BAM actually has access to shared randomness over
    preparations. §3a assumes it does, which is standard operationally but
    is an assumption about the program rather than a result of it.
  - The multipartite and higher-input cases, where the quantum set has no
    single-inequality characterization and this analysis does not extend.
  - Whether the zero-marginal restriction is related to the #208
    charged-GHZ superselection no-go, which also removes marginal-carrying
    sectors. Untested here.

## 9. Reproduce

```bash
python -m experiments.closure_ledger.admissible_correlation_tables_probe
```

Expected verdict:
`THE_REPRESENTATION_ADMITS_THE_GRAM_TLM_BODY_AND_SO_DOES_BAM_ONCE_MIXING_IS_ALLOWED_BECAUSE_THE_CONVEX_HULL_OF_THE_C_EQUALS_ONE_COSINE_TABLES_IS_THE_WHOLE_BODY_SO_THE_ONLY_SURVIVING_RESTRICTION_IS_IDENTICALLY_ZERO_MARGINALS`, 9/9 PASS (~3 min).

## 10. Cross-references

  - `docs/tsirelson_ceiling_research_plan.md` — #236, whose ceiling this
    extends from the maximum to the whole body, and whose T3 this scopes
  - `docs/configuration_space_emergence.md` — #206, the lattice and the
    state
  - `docs/tensor_product_emergence.md` — #232, the encoding being
    characterized

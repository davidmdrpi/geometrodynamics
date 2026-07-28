# What selects the twist sector? Energetics select, topology freezes (PR #230)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

> ### ⚠️ PARTIALLY SUPERSEDED — same PR. Read this first.
>
> §1.1's rule "**the sign is the statistics of the link field**" and
> §1.3's cross-check are **wrong**, corrected by
> `docs/bam_spinor_spectrum_effective_action_research_plan.md`. That
> computation takes BAM's actual Pin⁻ spinor instead of asserting its
> statistics, and finds the network holonomy composes with the
> **intrinsic B2 spin structure** (`T² = −I`, so `η = −1`), which **swaps
> the spinor's moding** relative to the scalar. The statistics flip and
> the spin-structure flip then cancel: `ΔE = +π/(4C)` for the BAM spinor
> too — **untwisted favoured in both sectors**.
>
> The twisted sector does carry a zero eigenvalue (`dim ker D = 1`,
> `det D = 0`), but that is **not** an index: it is the constant mode,
> lifted by any mass or potential, with index 0. So why matter sits in the
> twisted sector is **open** — neither energetics nor the zero mode
> explains it.
>
> **What survives unchanged:** the magnitude `π/(4C)` (§1.1), and the
> entire topological-freeze half — the amplitude-zero gate (§1.2) and the
> deformation-invariance of `W`. Those never depended on statistics.

## 0. The last open item

PR #227 promoted the ER-link orientation Z₂ to a gauge field over the
network and took the link variables `ε_b` as **free labels**. #228 showed
the twist acts only through cycles. #229 made the Möbius identification
quantitative in its topological half. All three left the same thing open
(#227 §8.3): **nothing said which networks carry `W = −1`.**

This closes it, and with it the arc.

## 1. The answer, in two halves acting at different times

### 1.1 The twist costs energy — `π/(4C)` per degree of freedom

The #229 Casimir computation, re-read as a selection rule. The two
holonomy sectors differ in zero-point energy on a loop of circumference
`C` by exactly `π/(4C)` per degree of freedom (verified to 1e-12 at every
`C`). This plan originally argued the **sign** follows from the statistics
alone — `+½Σω` for a boson, `−½Σω` for a fermion:

~~| link field | `ΔE` | favoured |~~
~~| **bosonic** | `+π/(4C)` | untwisted |~~
~~| **fermionic** | `−π/(4C)` | twisted |~~

**The magnitude `π/(4C)` is solid and is what survives.** The *sign* rule
above is **withdrawn**: it flipped the zero-point sign while leaving the
**moding** alone, but for a spinor the holonomy composes with the
intrinsic spin structure. With the moding computed rather than assumed:

| field | `η` | `ΔE = E₀(W=−1) − E₀(W=+1)` | favoured |
|---|---:|---:|---|
| scalar | `+1` | `+π/(4C)` | **untwisted** |
| **BAM Pin⁻ spinor** | `−1` | `+π/(4C)` | **untwisted** |

Energetics favour the untwisted sector **universally**. See
`docs/bam_spinor_spectrum_effective_action_research_plan.md`.

### 1.2 But the sector cannot relax into the favoured one

`W` is a discrete invariant, and the only continuous path between the
sectors goes through a hole. Interpolating the closing-link amplitude
`t → t(1−2s)`:

| `s` | closing amplitude | `b₁` | holonomy |
|---:|---:|---:|---|
| 0.00 | +1.00 | 1 | `+1` |
| 0.25 | +0.50 | 1 | `+1` |
| **0.50** | **0.00** | **0** | **UNDEFINED** |
| 0.75 | −0.50 | 1 | `−1` |
| 1.00 | −1.00 | 1 | `−1` |

At `s = ½` the amplitude is exactly zero: the cycle is **cut**, `b₁` drops
1 → 0, the ring becomes an open chain, and there is no closed loop left to
measure a holonomy on. So `W` cannot change by continuous deformation
without a topology change — **the same amplitude-zero gate #175 found for
the winding number**, now for the orientation Z₂.

And `W` is exactly invariant under any deformation that keeps the links
nonzero. The diagnostic: the untwisted ring has the constant vector in its
Laplacian kernel for *any* positive weights; the frustrated ring has no
zero mode at all. Across 200 random re-weightings per sector, with weights
spanning a factor of ~7:

  - untwisted: lowest eigenvalue ≤ **1.6e-15** (machine zero), always
  - twisted: lowest eigenvalue ≥ **4.8e-02**, never closes

**The energetic preference therefore does not open a relaxation channel.**
It sets which sector is cheaper *at nucleation*; after that the sector is
frozen.

### 1.3 ~~The check that makes the sign non-trivial~~ — **superseded**

> The table below was the original cross-check: it read the fermionic sign
> as confirming that matter sits in the twisted sector. **The sign was
> wrong, so the check does not stand.** Both sectors still come out right,
> but by *two different mechanisms* — energy for the bosonic QCD sector,
> the **zero-mode index** for the fermionic matter sector. Retained as the
> superseded record.

| sector | link field | statistics | predicted | what the repo already does |
|---|---|---|---|---|
| **QCD Möbius** (#100/#103) | flux-tube transverse displacement (phonon) | bosonic | untwisted favoured ⟹ Möbius states are **excitations above** the ground | presents them as an extra tower above the orientable one: `+πσ` (glueball), `+2√σ = 849 MeV` (heavy baryon) |
| **BAM matter** (#170/#195/#202) | throat modes are Pin⁻ spinors | fermionic | twisted favoured ⟹ matter sits **in** the twisted sector | #202's Pin parity selects the electron `k=1` mode **odd** at the neck; #227's T6 finds constructive antipodal self-return in the twisted odd-`l` sector |

Neither reading was arranged for this: both predate the probe.

**Honest limit:** this is a consistency check on the *sign*, not an
independent derivation of either sector's content.

## 2. Tests

| # | test | outcome |
|---|---|---|
| T1 | goal: the last open item | stated |
| T2 | the twist costs energy | `π/(4C)` to 1e-12 (**sign rule withdrawn — see banner**) |
| T3 | the sector cannot relax | gate at `s=½`: `b₁` 1→0, holonomy undefined |
| T4 | `W` is deformation-invariant | 200 re-weightings/sector; gap never closes |
| T5 | ~~consistency via the statistics sign~~ | **superseded** by the spinor computation |
| T6 | assessment | 5/5 |

## 2a. The corrected picture: two mechanisms, not one sign

| sector | selected by | why |
|---|---|---|
| **QCD Möbius** (#100/#103) | **energy** | untwisted cheaper by `π/(4C)` ⟹ Möbius states are excitations *above* the ground, as #100/#103 present them |
| **BAM matter** (#170/#195/#202) | **unexplained** | energetics favour untwisted here too, so energy cannot explain it — and the zero eigenvalue in that sector is the constant mode (index 0, lifted by any perturbation), so it cannot either |

Two candidate explanations have been tried and both failed: the statistics
sign (wrong moding) and the zero-mode index (not an index). The question is
narrowed, not answered.

## 3. What selection does **not** explain

  - **Absolute populations.** Energetics give the *ordering*, not the
    abundance. The population of each sector is set by the nucleation
    ensemble, which this probe does not compute.
  - **Why a given network nucleated as it did.** Selection is a
    preference, not a mechanism.
  - **The coefficient.** `π/(4C)` is per degree of freedom on a free loop;
    counting conventions (polarizations, spinor components) change the
    coefficient.
  - **And the sign was not safe either** — see the banner. Reading a
    selection sign off "statistics" without computing the moding is
    exactly the error this plan originally made.

## 4. Cross-references

  - `docs/nonorientable_er_link_z2_gauge_research_plan.md` — #227, the
    gauge field and the §8.3 item this closes
  - `docs/mobius_identification_quantitative_research_plan.md` — #229, the
    Casimir computation re-read here as a selection rule
  - `docs/nonlinear_antipodal_focusing_pde_research_plan.md` — #175, the
    amplitude-zero gate this reproduces for the orientation Z₂
  - `docs/odd_k_generation_survival_research_plan.md` — #183, the
    deformation-invariance of the orientability class, reached here from
    the network side
  - `docs/glueball_closed_flux_loop_research_plan.md` — #100, the bosonic
    Möbius sector; `docs/pin_rp2_fermi_statistics_research_plan.md` — #170,
    the fermionic one

## 5. Reproduce

```bash
python -m experiments.closure_ledger.twist_sector_selection_probe
```

Expected verdict:
`THE_TWIST_SECTOR_IS_SELECTED_ENERGETICALLY_WITH_A_SIGN_SET_BY_THE_LINK_FIELD_STATISTICS_AND_THEN_TOPOLOGICALLY_FROZEN_AT_NUCLEATION`, 5/5 PASS.

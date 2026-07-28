# What selects the twist sector? Energetics select, topology freezes (PR #230)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## 0. The last open item

PR #227 promoted the ER-link orientation Z₂ to a gauge field over the
network and took the link variables `ε_b` as **free labels**. #228 showed
the twist acts only through cycles. #229 made the Möbius identification
quantitative in its topological half. All three left the same thing open
(#227 §8.3): **nothing said which networks carry `W = −1`.**

This closes it, and with it the arc.

## 1. The answer, in two halves acting at different times

### 1.1 The twist costs (or pays) energy — and the sign is the statistics

The #229 Casimir computation, re-read as a selection rule. The two
holonomy sectors differ in zero-point energy on a loop of circumference
`C` by exactly `π/(4C)` per degree of freedom (verified to 1e-12 at every
`C`). The **sign** depends on what lives on the link, because the zero
point is `+½Σω` for a boson and `−½Σω` for a fermion:

| link field | `ΔE = E₀(W=−1) − E₀(W=+1)` | favoured sector |
|---|---:|---|
| **bosonic** | `+π/(4C)` | **untwisted** |
| **fermionic** | `−π/(4C)` | **twisted** |

Same magnitude, opposite sign. So there is **no universal answer** to
"which sector is the ground state" — it is a property of the field on the
link. That is precisely what makes §1.3's cross-check meaningful rather
than circular: the sign is not tunable.

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

### 1.3 The check that makes the sign non-trivial

The statistics-dependent sign is not a free parameter, and the program
already contains **two Möbius sectors of opposite character**. One Casimir
sign gets both right:

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
| T2 | the twist costs energy; sign = statistics | `π/(4C)` to 1e-12; signs opposite |
| T3 | the sector cannot relax | gate at `s=½`: `b₁` 1→0, holonomy undefined |
| T4 | `W` is deformation-invariant | 200 re-weightings/sector; gap never closes |
| T5 | consistency with both Möbius sectors | both correct |
| T6 | assessment | 5/5 |

## 3. What selection does **not** explain

  - **Absolute populations.** Energetics give the *ordering*, not the
    abundance. The population of each sector is set by the nucleation
    ensemble, which this probe does not compute.
  - **Why a given network nucleated as it did.** Selection is a
    preference, not a mechanism.
  - **The coefficient.** `π/(4C)` is per degree of freedom on a free loop;
    counting conventions (polarizations, spinor components) change the
    coefficient. The **sign** — the part doing the work — does not.

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

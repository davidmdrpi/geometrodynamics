# Probe G — the freeze and deck class 1 are **independent** (PR #233)

_Run 2026-08-03T06:44:47.560325+00:00 · 7/7 PASS_

> **Correction.** The first version of this probe claimed a parity
> obstruction making the two mutually exclusive. That claim assumed
> every exterior arc has length π — which forces all throats to
> share one antipodal mouth pair — and is retracted here.

## The two laws

| law | status |
|---|---|
| full freeze ⟺ `N` even (`S^N = W`) | **stands** |
| deck class = (Σ arc lengths)/π mod 2 | **set by mouth placement, not by `N`** |

## The counterexample

| ring | arcs/π | Σ/π | deck | doublet | gap untwisted | gap twisted | frozen? |
|---|---:|---:|---:|---|---:|---:|---|
| `N` = 2 | 0.5 | 1.0 | 1 | yes | 1.908e-02 | 8.482e-13 | **YES** |
| `N` = 2 | 1.5 | 3.0 | 1 | yes | 8.266e-03 | 3.297e-12 | **YES** |
| `N` = 4 | 0.75 | 3.0 | 1 | yes | 4.705e-04 | 2.093e-12 | **YES** |

## The other evasion classes

| class | result |
|---|---|
| one throat, finite channel (`D > 0`) | nondegenerate interior ladder — no doublet, no freeze |
| inhomogeneous / internally structured | freeze **destroyed** even at even `N` |
| non-cyclic (theta network) | no freeze in any twist sector |

## Consequences

  - *the freeze law: full freeze iff N even (S^N = W)* — **STANDS**
  - *the earlier probe's deck-class law, deck class = N mod 2* — **CORRECTED**
  - *the earlier probe: freeze and deck class 1 are mutually exclusive by parity* — **RETRACTED**
  - *#231's recommended next construction (a geometry carrying both the doublet and deck class 1)* — **RESTORED, in corrected form**
  - *the freeze sector is cut off from RP^3 spin-structure data permanently, not contingently* — **RETRACTED**
  - *Probe F (#231): the phenomenon and the coupling never co-occur* — **DEMOTED back to an observation about the two rings the repo contains**

## Verdict

**THE_FREEZE_NEEDS_EVEN_N_BUT_THE_DECK_CLASS_IS_SET_BY_MOUTH_PLACEMENT_NOT_BY_N_SO_THEY_ARE_INDEPENDENT_AND_A_TWO_THROAT_RING_WITH_QUARTER_CIRCLE_ARCS_CARRIES_BOTH**

# The `k ≠ 0, n ≥ 3` quadrant: winding shell modes = leptoquark sector (PR #85)

PR #83 unified the lepton and quark mass operators into one
Bohr-Sommerfeld operator

```
m²(k, n)  =  (k·2π / L_throat)²  +  ((n+1)·π / L_cavity)²,   L_throat = √(2π)/k_5,
```

and flagged as open: *"prediction of new states — e.g. winding shell
modes with both `k ≠ 0` and `n ≥ 3`, if physical."* This probe maps
the full `(k, n)` lattice and shows the previously-empty
`(k ≠ 0, n ≥ 3)` quadrant is the **leptoquark sector**.

## The four quadrants

The two channels (winding `k`, radial overtone `n`) split the lattice
into four quadrants — one sector of each per generation `g ∈ {1,2,3}`:

| sector | (k, n) | winds? | resolves cavity? | character |
|---|---|---|---|---|
| neutrino candidate | (0, g−1) | no | partial | neither |
| quark | (0, g+2) | no | yes | cavity-resolving |
| charged lepton | (2g−1, 0) | yes | no | throat-winding |
| **leptoquark** | **(2g−1, g+2)** | **yes** | **yes** | **both** |

The raw-operator `m²` values (BAM units):

| g | neutrino (0,g−1) | quark (0,g+2) | charged lep (2g−1,0) | leptoquark (2g−1,g+2) |
|---:|---:|---:|---:|---:|
| 1 | 1.11 | 14.63 | 158.19 | 171.71 |
| 2 | 3.90 | 22.67 | 1414.83 | 1436.38 |
| 3 | 8.38 | 32.49 | 3928.10 | 3959.48 |

## The leptoquark quadrant

The `(k ≠ 0, n ≥ 3)` quadrant carries **both** throat-winding (lepton
character, `k ≠ 0`) and cavity-resolution (quark character, `n ≥ 3`).
In the unified operator both mass terms add, so the leptoquark is the
**heaviest state in each generation**.

This is exactly Pati-Salam SU(4)'s prediction: the `SU(4)/SU(3)×U(1)`
coset states are leptoquarks — quark↔lepton converters carrying both
quark color and lepton number. The BAM leptoquark quadrant is the
operator-level realization of the structural Pati-Salam bridge built
in PR #82.

## The candidate neutrino quadrant

The complementary `(k = 0, n < 3)` quadrant gives the lightest states
(non-winding, throat-region) — candidate neutrinos, partially closing
one of PR #82's three open extensions.

**Honest caveat:** the BAM-operator ν/charged-lepton mass ratio is
~0.07, far above the observed `< 10⁻⁶`. A genuine neutrino
identification needs an additional suppression (Majorana seesaw, or a
special `k=0, n=0` structure). Flagged as open.

## Mass ordering

Within each generation the robust, scaling-independent endpoints are:

  - **neutrino candidate** lightest (neither winds nor fully resolves
    the cavity);
  - **leptoquark** heaviest (does both).

The charged-lepton vs quark middle ordering is a raw-operator artifact
(the raw operator lacks the v3 `(k−3)²` uplift and proper scaling), but
the endpoints do not depend on scaling.

## What this establishes

  1. The unified `(k, n)` operator has a complete **four-quadrant
     interpretation** — neutrino / quark / charged lepton / leptoquark,
     one of each per generation, matching the Pati-Salam content.
  2. The `(k ≠ 0, n ≥ 3)` quadrant flagged by PR #83 is the
     **leptoquark sector** (both characters, heaviest per generation).
  3. Connects directly to the Pati-Salam SU(4) bridge (PR #82): the
     leptoquark quadrant is the `SU(4)/SU(3)` coset.
  4. Surfaces a **candidate neutrino sector** (the `k=0, n<3`
     quadrant), partially addressing PR #82's missing-neutrino
     extension (with the mass-scale caveat).

## Falsifiable prediction

BAM has a **4th matter sector** — leptoquarks at `(k ≠ 0, n ≥ 3)`,
heaviest in each generation because both mass terms add. Their
non-observation is consistent with them being heavy. The unified
operator extends from "leptons + quarks" to "neutrinos + quarks +
charged leptons + leptoquarks" — a complete generation multiplet
matching the Pati-Salam content.

## Honest scope

  - **Is:** the four-quadrant map of the unified `(k, n)` operator; the
    structural identification of `(k ≠ 0, n ≥ 3)` as the leptoquark
    sector (both characters, heaviest per generation); the Pati-Salam
    coset connection; the candidate neutrino quadrant; the
    within-generation mass ordering (robust endpoints).

  - **Is not:** a prediction of absolute leptoquark masses (the
    absolute scale needs the L_eff unification still open from PR #83,
    plus the B4 anchor; the raw operator without the v3 `(k−3)²` uplift
    does not give physical MeV values); a claim that these states are
    observed; a derivation of the neutrino mass scale (the `k=0, n<3`
    states are light but not neutrino-light); a spin/statistics
    assignment of the leptoquarks (odd-`k` = fermion by PR #67;
    Pati-Salam leptoquark gauge bosons would need the even-`k` bosonic
    sector, left open).

## B4 accounting

The quadrant map and orderings are scale-free. Absolute masses ride on
the L_eff unification (PR #83 open) plus the single B4 anchor (PR #53).

## Tests

  T1. Map the `(k, n)` lattice; four quadrants.
  T2. Recover charged leptons (winding-dominated) and quarks
      (cavity-only).
  T3. `(k ≠ 0, n ≥ 3)` leptoquark quadrant: both terms add, heaviest
      per generation.
  T4. Pati-Salam connection: leptoquark = `SU(4)/SU(3)` coset.
  T5. `(k = 0, n < 3)` candidate neutrino quadrant; mass-scale caveat.
  T6. Within-generation mass ordering: neutrino lightest, leptoquark
      heaviest.
  T7. Falsifiable predictions + honest scope.
  T8. Assessment.

## Verdict structure

  - **`WINDING_SHELL_QUADRANT_IS_LEPTOQUARK_SECTOR`** (expected): the
    `(k ≠ 0, n ≥ 3)` quadrant is the leptoquark sector; the unified
    operator has a complete four-quadrant interpretation.
  - **`QUADRANT_MAP_INCONCLUSIVE`**: a structural test fails.

## What this leaves open

  - **Absolute leptoquark masses** — L_eff unification (PR #83) + B4
    anchor.
  - **Leptoquark spin/statistics** — odd-`k` fermion vs even-`k` boson
    (PR #67).
  - **Neutrino mass scale** — extra suppression needed.
  - **Observability / stability** of the leptoquark sector.

## Cross-references

  - `docs/throat_shell_mass_operator_unification_research_plan.md` —
    PR #83, the unified operator and the flagged `(k≠0, n≥3)`
    prediction.
  - `docs/pati_salam_throat_shell_bridge_research_plan.md` — PR #82,
    the Pati-Salam SU(4) bridge and the missing-neutrino extension.
  - `docs/beta_lepton_derivation_research_plan.md` — PR #71, the
    charged-lepton winding sector.
  - `docs/qcd_shell_waveguide_scaffold_research_plan.md` — PR #77, the
    quark cavity sector.
  - `docs/even_k_absence_research_plan.md` — PR #67, odd-`k` = fermion
    (relevant to leptoquark spin/statistics).
  - `experiments/closure_ledger/winding_shell_quadrant_probe.py` —
    this probe.

# The absolute-coupling capstone: canonical Hopf–KK normalization, geometric α, and the global no-retuning holdout (PR #225, FINAL)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. The capstone question: *does
> canonical dimensional reduction and charge normalization uniquely
> determine the four-dimensional electromagnetic coupling, or does a
> continuous dimensionless modulus remain?* The companion probe
> machine-checks every claim (~10 s).

## 0. The answer, stated first

**A modulus remains.** The canonical chain determines the coupling
*exactly* as

```
α_k  =  4 k² (l_P / R_f)²  =  4 k² / ρ²,        ρ = R_f / l_P
```

— a function of exactly **one continuous dimensionless modulus**, the
Hopf-fiber radius in Planck units (the radion). Canonical kinematics is
invariant under rescaling ρ: **uniqueness fails at the canonical
level**, and α is not a kinematic pure number. No amount of
normalization discipline changes this — the honest general answer to
the capstone question. What the program adds is the *correctly located*
dynamical stabilizer: the EM cap selects ρ\* = 2/√α ≈ 23.4, and that
insertion passes the global no-retuning holdout over the entire arc.

## 1. The canonical chain, assembled

The 5D Einstein–Hilbert action on the Hopf-fibered ansatz

```
ds₅² = g_μν dx^μ dx^ν + R_f²(dθ + A_μ dx^μ)²
```

reduces to 4D gravity + Maxwell: G₄ = G₅/(2πR_f); the F² coefficient
R_f²/4·(2πR_f/16πG₅); canonical normalization A_c = A·R_f/√(16πG₄);
and the fiber-winding-k mode e^{ikθ} carries charge

```
q_k = k·√(16πG₄)/R_f      ⟹      α_k = q_k²/4π = 4k²·G₄/R_f².
```

Machine-checked structural content:

- **the fiber KK tower** m_k = k/R_f: FD ring spectra = the discrete
  closed form to 10⁻¹¹, continuum to 4×10⁻⁴;
- **the #193 weld**: the Berger closed form
  E(j,k;λ) = 4j(j+1) − k² + k²/λ² carries exactly the same fiber term,
  with the sector ground splitting as monopole zero-point (2k, the
  Wu–Yang charge q = k/2) plus the fiber tower — the identity is exact:
  **the #193 spectra are the KK tower of this chain**, and the
  half-charge is the Hopf half-radius fiber statement.

## 2. Charge = fiber winding, canonically — no freedom

Threading flux η through the fiber flows the spectrum as
((k + η)/R_f)²:

- the tracked-branch **spectral-flow slope** d(ω²)/dη = 2(k+η)/R_f²
  holds to 3×10⁻⁴ (discrete dispersion included) — the charge of the
  k-mode *is* k in canonical units;
- the spectrum is **exactly periodic in one flux quantum** (4×10⁻¹¹) —
  large gauge invariance: the normalization is fixed by topology, not
  by choice;
- the **dynamical force check**: a slowly threaded flux (an EMF around
  the fiber) chirps the k-mode at the predicted instantaneous
  frequency |k + η(t)|/R_f to 3×10⁻⁴ — the electric force acting on
  charge k, *evolved* rather than asserted.

## 3. The modulus

The twisted-tower-plus-flow data at two fiber radii are related by
**exact scaling** (machine zero): no canonical equation selects R_f.
This is the precise mirror image of #222's finding — there the
*coupled field equations* killed the rescale freedom (the weld); here
the *vacuum kinematics* keeps it (the radion flat direction). Both
statements are exact, and together they bracket what geometry alone
can and cannot fix.

## 4. The stabilizer and the guardrail

The EM cap (#55–#58; primordial per #222) is the program's dynamical
stabilization of the radion:

```
ρ*  =  2/√α  =  23.4125          (k = 1: the fiber at ~23 l_P)
```

The #165 guardrail scan finds **no closure-constant match** for ρ\*
(nearest: e^π at 1.2%, k₅² = 25 at 6.8% — both rejected): the
selection is dynamical, not numerological, and deriving ρ\* from the
cap's own equations is the program's open dynamical problem —
**correctly located outside canonical kinematics by this capstone**.

## 5. The global no-retuning holdout

Eight keystone constants of the arc, re-read from the **committed run
ledgers** and independently re-derived where exact:

| constant | ledger | independent | dev |
|---|---|---|---|
| z\* (z·J₁ = 3J₂, #223) | 2.299481 | fresh Bessel root | <10⁻⁶ |
| the π/2 interior-depth step (#223) | family values | π/2 | <10⁻³ |
| the quarter wave (#221) | π/2 closed form | π/2 | exact |
| r_sω supremum = √μ_crit (#222) | ledger pair | identity | 10⁻⁹ |
| the Rabi identity (#224) | ≤2% stored | 1 | ✓ |
| the Weyl commutator (#160) | — | e^{2πi/5}·I live | 10⁻¹⁶ |
| φ_h = π/k₅ | exact | exact | 0 |
| β_lepton = 50π | exact | exact | 0 |

All are **ratios, roots, or topological phases — none is a function of
ρ**: fixing ρ\* = 2/√α by the EM cap retunes *nothing*. The holdout
passes globally: the arc's books are balanced.

## 6. Honest scope

- The answer is kinematic and tree-level; quantum corrections and
  running are not addressed — α here is the boundary coupling that
  #184 protects.
- The modulus verdict is the honest one; the EM cap is the named
  dynamical stabilizer, its derivation inheriting #55–#58's scope.
- The Hopf reduction is at the round point; the Berger squashing is
  part of the modulus space (the #192–#197 dynamics on it is the
  already-mapped story).
- l_P enters through G₄ (the B4 single-anchor discipline unchanged);
  ℏ is still not derived.
- Factor conventions are fixed against the #193 machine-validated
  spectra, not chosen.
- The holdout is an integration test of mutual consistency under the
  α(ρ\*) insertion, not a re-derivation of every probe.

## 7. What would falsify this

- A canonical equation fixing R_f — the modulus claim would be wrong.
  (Checked: exact rescale invariance of the kinematics.)
- A charge normalization with residual freedom — flux quantization or
  the flow slope would drift. (Checked: 10⁻¹¹ periodicity; 3×10⁻⁴
  slope.)
- A ledger constant depending on ρ — the holdout would fail on
  insertion. (Checked: eight constants, all ρ-independent, all
  matching their independent re-derivations.)
- A closure-constant match for ρ\* — the stabilizer would be
  kinematic after all. (Checked: none within 1%.)

## 8. Companion probe

`experiments/closure_ledger/absolute_coupling_capstone_probe.py`
(T1–T9, ~10 s): the canonical chain and the #193 weld; the exact
spectral flow, flux quantization, and the adiabatic-ramp force check;
the α law and the modulus; the stabilizer and the guardrail scan; the
global holdout over the committed ledgers.

**Verdict:**
`CANONICAL_HOPF_KK_NORMALIZATION_GIVES_ALPHA_EQUALS_4K2_OVER_RHO2_EXACTLY_A_CONTINUOUS_RADION_MODULUS_REMAINS_THE_EM_CAP_IS_ITS_DYNAMICAL_STABILIZER_AND_THE_GLOBAL_NO_RETUNING_HOLDOUT_PASSES`

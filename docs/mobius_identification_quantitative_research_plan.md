# Is the Möbius identification quantitative, or only structural? (PR #229)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## 0. The question

PR #227 established that the ER-network **link twist** and the QCD
**Möbius label** are the same Z₂ *structurally*: same holonomy, same
integer → half-integer comb, and — the sharpest part — two
implementations the repo had been carrying independently without
connecting them (`make_glueball_ring`, `topology_kind='periodic'`,
`orientation=+1`; `make_mobius_tube`, `'mobius'`, `−1`).

#227 §8.2 left the follow-up open: is the **`+πσ`** glueball tower shift
of #100 *quantitatively* the same object, or only structurally?

## 1. The answer: partially — and the missing part is a gap in #100

The shift factorizes, and the two factors have very different standing:

```
Δ M²  =  (1/2)              ×   (2πσ)
         moding shift            closed-string M² spacing
         DERIVED, shared         INPUT, not supplied by the network
```

### 1.1 The `1/2` is derived, shared, and topological — **CONFIRMED**

#100's own T4 *asserts* the integer → half-integer shift: it returns
`pass: True` with no computation, recording only the string
`'integer n → half-integer n+½ (antiperiodic)'`. This probe computes it,
three independent ways, all agreeing on exactly ½ (T2):

| implementation | result |
|---|---|
| the #227 twisted ring | comb `0.5, 0.5, 1.5, 1.5` vs periodic `0, 1, 1, 2` |
| the repo's own `MobiusSpectrum` | `ω_mob/π = 0.5, 1.5, 2.5, 3.5` vs `ω_ord/π = 1, 2, 3, 4`; its numerical solve confirms |
| a direct antiperiodic closed-loop solve | `1/2, 1/2, 3/2, 3/2`, independent of circumference |

And it is **topological**, in the strict sense that it depends on neither
the geometry nor the discretization (T3):

| grid | antiperiodic ground `n` | dev from ½ |
|---:|---:|---:|
| 120 | 0.499985277 | 1.5e-05 |
| 240 | 0.499996430 | 3.6e-06 |
| 480 | 0.499999108 | 8.9e-07 |
| 960 | 0.499999777 | 2.2e-07 |

Successive ratios **4.00, 4.00, 4.00** — exactly `O(1/M²)`, the
finite-difference dispersion of the discrete Laplacian and nothing else.
Across a factor-of-5 change in circumference at fixed grid the deviation
moves by **1.4e-11**, four orders below the discretization error itself.
So the ½ is exact in the continuum and geometry-independent — which is
what makes the identification quantitative rather than merely suggestive.

*Measurement note:* the deviation is taken on the antiperiodic ground mode
directly, not as a difference against the periodic zero mode. That mode's
eigenvalue is ~1e-13 and the `sqrt` amplifies eigensolver noise into an
erratic ~1e-6 unrelated to the shift — the only ill-conditioned entry in
either comb.

### 1.2 The `2πσ` is **not** supplied by the network

The `M²` spacing comes from the closed-string Regge slope `1/(4πσ)`, a
string-theoretic relation #100 imports, and `σ` itself is the
lattice-calibrated QCD anchor. Nothing in the ER network supplies either.

So the honest claim is neither "the network derives `πσ`" nor "the
identification is only structural", but precisely: **the topological half
transfers exactly; the dimensionful half does not come from the network.**

### 1.3 #100's tower carries an unjustified constant — **a real gap**

#100 writes both towers with the **same** intercept:

```
orientable:  M²ₙ = M₀² + 2πσ·n
Möbius:      M²ₙ = M₀² + 2πσ·(n + ½)
```

But antiperiodic moding shifts the zero-point (Casimir) energy as well as
the level, so the Möbius intercept cannot be the orientable one. Computed
two independent ways — ζ-regularization using `ζ(−1) = −1/12` and
`ζ(−1,½) = +1/24`, and an exponential-cutoff extraction with the `1/ε²`
divergence removed analytically — agreeing to **4.1e-07** (T5):

| `C` | `E₀` periodic | `E₀` antiperiodic | `ΔE₀` | `π/(4C)` |
|---:|---:|---:|---:|---:|
| 1.0 | −0.523599 | +0.261799 | **+0.785398** | +0.785398 |
| 2.0 | −0.261799 | +0.130900 | **+0.392699** | +0.392699 |
| 3.5 | −0.149600 | +0.074800 | **+0.224399** | +0.224399 |

Per transverse polarization on a loop of circumference `C`:

```
E₀ periodic     = −π/(6C)
E₀ antiperiodic = +π/(12C)
ΔE₀             = +π/(4C)      (exact, both methods)
```

The Möbius zero point is **higher**, and it scales as `1/C` — verified by
halving `C` and watching it double — so it cannot be absorbed into a
common `M₀²`. **#100's Möbius ground state is therefore not
`M₀² + πσ`.**

Quantifying the resulting `M²` intercept needs the full closed-loop
dynamics that #100 itself defers, so this probe **reports the gap and
does not publish a corrected tower**. The level shift `+πσ` stands; the
common intercept does not.

## 2. Scope of the correction — it does **not** hit the search table

This touches the **glueball** tower of #100 only, which is built from the
closed-string mass formula. It does **not** touch the heavy Möbius
**baryon** predictions of #103/#109/#114 — Λ_c 3135, Λ_b 6469, the 849 MeV
dipion endpoint — which are built from the flux-tube quantum
`Δ = 2√σ` rather than the closed-loop intercept. **The search table stands
as published.**

Glueballs are unobserved (they mix with `qq̄` mesons), so the corrected
item is precisely the one with no experimental exposure. The net effect is
to **sharpen #100, not overturn it**.

## 3. Tests

| # | test | outcome |
|---|---|---|
| T1 | the factorization under test | stated |
| T2 | the `1/2` from three independent implementations | all exactly ½ |
| T3 | the `1/2` is topological | `O(1/M²)` exactly; circumference spread 1.4e-11 |
| T4 | the `πσ` ledger: derived factor vs inherited scale | `0.5 × 1.130973 = 0.565487 GeV²` |
| T5 | the intercept gap in #100 | `ΔE₀ = +π/(4C)`, two methods to 4.1e-07 |
| T6 | assessment and scope | 6/6 |

## 4. What this leaves open

  - **The corrected Möbius `M²` intercept.** Needs the full closed-loop
    dynamics #100 defers. Not published here — reporting a gap is not the
    same as fixing it.
  - **Whether the ER network can supply the `2πσ` spacing at all.** If it
    could, the identification would become fully quantitative rather than
    quantitative-in-one-factor.
  - **#227 §8.3** — what selects the twist sector dynamically — untouched.

## 5. Cross-references

  - `docs/nonorientable_er_link_z2_gauge_research_plan.md` — #227, the
    structural identification this makes quantitative (§8.2 of its roadmap)
  - `docs/glueball_closed_flux_loop_research_plan.md` — #100, the `+πσ`
    tower whose `1/2` is confirmed and whose intercept is corrected
  - `docs/heavy_mobius_baryon_search_table.md` — #114, the search table
    this correction explicitly does **not** touch
  - `geometrodynamics/qcd/spectrum.py` — `MobiusSpectrum`, the repo's own
    half-integer machinery, used as an independent check
  - `geometrodynamics/qcd/topology.py` — `make_glueball_ring` /
    `make_mobius_tube`, the two holonomy sectors

## 6. Reproduce

```bash
python -m experiments.closure_ledger.mobius_identification_quantitative_probe
```

Expected verdict:
`THE_MOEBIUS_IDENTIFICATION_IS_QUANTITATIVE_IN_ITS_TOPOLOGICAL_HALF_THE_SCALE_IS_INHERITED_AND_THE_MOEBIUS_INTERCEPT_IS_AN_OPEN_CORRECTION`, 6/6 PASS.

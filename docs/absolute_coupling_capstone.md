# The absolute-coupling capstone: canonical Hopf–KK normalization, geometric α, the Einstein-frame radion, and the α-dependent holdout (PR #225, FINAL; revised)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. The capstone question: *does
> canonical dimensional reduction and charge normalization uniquely
> determine the four-dimensional electromagnetic coupling, or does a
> continuous dimensionless modulus remain?* The companion probe
> machine-checks every claim (~5 s), and now runs in CI.

## 0. The answer, stated first

**A modulus remains.** The canonical chain determines the coupling
*exactly* as

```
α_k  =  4 k² (l_P / R_f)²  =  4 k² / ρ²,        ρ = R_f / l_P
```

— a function of exactly **one continuous dimensionless modulus**, the
Hopf-fiber radius in Planck units (the radion). Canonical kinematics is
invariant under rescaling ρ: **uniqueness fails at the canonical
level**, and α is not a kinematic pure number. The revision sharpens
every leg of this answer: the modulus is now a *canonical 4D field*
(the Einstein-frame radion, flat at tree level, with α(φ) = α(0)e^{√3φ}
its dilaton coupling); the fiber is *explicitly* the committed
geometry (R_f = λR_u); the rank audit proves the EM cap fixes *exactly
the one α-coupled direction*; the implied KK scale is lepton-compatible
by decoupling; and the holdout is now genuinely α-dependent — twelve
committed constants, one common α.

## 1. The canonical chain, on the committed geometry

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
  Wu–Yang charge q = k/2) plus the fiber tower — **the #193 spectra
  are the KK tower of this chain**;
- **the explicit geometry map** (revision): on the committed Berger
  metric g = (R_u²/4)[σ₁² + σ₂² + λ²σ₃²] the Hopf fiber has proper
  length 2πλR_u by quadrature — **R_f = λ·R_u**, the round point
  R_f = R_u — while the KK-monopole base is the *half-radius* sphere
  (area πR_u², the geometric home of the Wu–Yang half-charge). With
  the committed Tangherlini conventions R_u = r_h = 1 (#216–#224),
  the EM cap converts the model unit:
  **1 model unit = ρ\*·l_P = 23.41 l_P** — the committed throat is a
  ~23-Planck-length object.

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

## 3. The full Einstein-frame reduction with the radion (revision)

Promote the fiber radius to a field with the two-parameter Weyl ansatz

```
ds₅² = e^{2aφ} g_μν dx^μ dx^ν + e^{2bφ} R₀²(dθ + A_μ dx^μ)²
```

and reduce **symbolically** (sympy, 5×5 Ricci scalars on test metrics
that isolate each term of the reduced action). Machine-checked:

- the 4D curvature terms carry the prefactor **e^{(2a+b)φ}**: the
  Einstein frame is **b = −2a** and nothing else;
- the radion kinetic coefficient is **−6a²** (the φ″ remainder an
  exact total derivative) — canonical −½ at **a = 1/(2√3)**;
- the gauge-kinetic coefficient is **−(R₀³/2)e^{3bφ}** — the standard
  KK dilaton coupling **e^{−√3φ}F²**;
- the tree-level radion potential is **exactly zero**: the modulus of
  §0 is the flat direction of a genuine canonical 4D field;
- the closing identity: the dilaton exponent √3 **equals** the
  geometric-law exponent 6a of 1/R_fE² (the fiber proper radius
  referred to Einstein-frame lengths, R_fE = e^{−3aφ}R₀):

```
α(φ) = α(0)·e^{√3 φ}  =  4k²(l_P/R_fE)²
```

**The α law and the dilaton coupling are one statement.** The EM cap
is thereby a genuinely *dynamical* claim — where ⟨φ⟩ sits on its flat
tree-level potential — not a kinematic assertion.

## 4. The modulus, and the full rank audit (revision)

The twisted-tower-plus-flow data at two fiber radii are related by
**exact scaling** (machine zero): no canonical equation selects R_f.
The **modulus-rank audit** makes the counting exact. The committed
arc's continuous dimensionless knobs, in log space:

| knob | meaning |
|---|---|
| ln(R_u/l_P) | the universe scale in Planck units |
| ln λ | the Berger squashing (R_f = λR_u) |
| ln(r_h/R_u) | the throat-to-universe ratio |
| ln(R_RMS/r_s) | the soliton-to-throat ratio |
| ln A | the dressing amplitude (#223) |

with α = 4k²/ρ², ρ = λ·R_u/l_P, so ∇ln α = (−2, −2, 0, 0, 0).
Committed constraint rows: the **#222 weld** (fixes ln(R_RMS/r_s));
the **EM cap** (fixes the combination ln λ + ln(R_u/l_P) = ln ρ).
Numpy rank:

- **before the cap**: rank 1, four flat directions, and ∇α is *not*
  annihilated by the flat space (projection 2.83) — α undetermined;
- **after the cap**: rank 2, three flat directions, and ∇α is
  annihilated to machine zero (2×10⁻¹⁶) — **the cap fixes exactly
  the one α-coupled direction**, and every remaining flat —
  squashing-at-fixed-fiber (λ observationally pinned to (0.986, ≥3]
  by #192, not dynamically fixed), r_h/R_u (a background datum), and
  A (a matter dof) — is **α-decoupled**. This is the rank-proven
  form of no-retuning.

## 5. The stabilizer, the guardrail, and the lepton sector

The EM cap (#55–#58; primordial per #222) is the program's dynamical
stabilization of the radion:

```
ρ*  =  2/√α  =  23.4125          (k = 1: the fiber at ~23 l_P)
```

The #165 guardrail scan finds **no closure-constant match** for ρ\*
(nearest: e^π at 1.2%, k₅² = 25 at 6.8% — both rejected): the
selection is dynamical, not numerological, and deriving ρ\* from the
cap's own equations is the program's open dynamical problem.

**Lepton-sector compatibility** (revision): the KK mass the cap
implies is

```
m₁ = 1/R_f = (√α/2)·m_P = 5.2×10¹⁷ GeV
```

— **5×10¹⁸ above m_μ** (10²¹ above m_e): leptons are *not* fiber KK
modes, and do not need to be. The committed lepton mechanism is the
throat-cavity sector (#221/#223), which carries α in *ratios*
(m_e/m_μ = α/X); and m_lepton ≪ m_KK is exactly the zero-mode
decoupling the 4D Maxwell chain assumes — **compatible by
decoupling**. Structurally: the #197 winding-ladder slope 1/λ *is*
the fiber-tower slope 1/R_f (same fiber radius, machine-checked to
4×10⁻⁵ across λ), and β_lepton/L_fiber = 50π/2π = 25 = k₅² exactly —
a topology statement, ρ-independent.

## 6. The α-dependent holdout (revision)

The earlier holdout checked ρ-independent constants (roots, ratios,
phases); per the revision it is **replaced** by a holdout over twelve
committed ledger constants that are *genuine functions of α* —
linear, quartic, and inverse-quartic — each re-derived from its
independent ledger inputs plus α = 7.2973525693×10⁻³:

| constant | α-dependence | ledger | dev |
|---|---|---|---|
| #221 X_required (conv B) | α·(m_μ/m_e) | 1.508861 | 0 |
| #221 m_e/m_μ hard wall | 2α/π | 0.0046456 | 0 |
| #221 landing (hard wall) | (2α/π)(m_μ/m_e) | 0.96057 | 10⁻¹⁶ |
| #221 landing (reference) | α/X_ref·(m_μ/m_e) | 0.95600 | 10⁻¹⁶ |
| #222 anchor r_sω | α itself | 0.0072974 | 0 |
| #222 anchor exclusion | band_lo/α | 209.89 | 0 |
| #222 X exclusion | X₂₀₆/(α·m_μ/m_e) | 380.14 | 0 |
| #223 X_required (conv B) | α·(m_μ/m_e) | 1.508861 | 0 |
| #223 m_e/m_μ class max | 2α/π | 0.0046456 | 0 |
| #223 structural residual | 1 − 2α(m_μ/m_e)/π | 0.039429 | 0 |
| #223 transit protection | ∝ α⁴ | 3.51×10⁻¹² | 0 |
| #224 frozen period | ∝ α⁻⁴ | 1.19×10¹¹ | 0 |

All twelve re-derive at machine zero; **inverted, the twelve inferred
α's agree to 4×10⁻¹⁶ relative** — one common α runs the whole arc,
and fixing ρ\* = 2/√α retunes nothing that carries α. (The Bessel
universal of #223 remains as committed: z\* = 2.29991…, the root of
z·J₁(z) = 3·J₂(z) — the earlier draft of this document misprinted
its digits.)

## 7. Honest scope

- The answer is tree-level: canonical reduction and canonical charge
  normalization, now including the full Einstein-frame radion with
  V_tree = 0. Quantum corrections, a radion potential from loops or
  matter, and running are not addressed — α here is the boundary
  coupling that #184 protects.
- The Einstein-frame reduction is machine-checked on test metrics
  (a warped 4D probe; a constant-radion gauge sector) that isolate
  each term; it is not a symbolic reduction of the fully generic 5D
  metric. The checked coefficients are exactly the standard
  KK-dilaton results.
- The rank audit enumerates the knobs the committed arc actually
  carries; λ is pinned observationally (#192), r_h/R_u is a
  background datum — both α-decoupled, neither derived.
- Leptons are compatible by decoupling, not identification: no
  committed number identifies a lepton with a fiber KK level.
- l_P enters through G₄ (the B4 single-anchor discipline unchanged);
  ℏ is still not derived.
- Factor conventions (fiber period, q = k/2, R_f = λR_u with the
  half-radius *base*) are fixed against the #193 machine-validated
  spectra and the explicit Berger quadrature, not chosen.
- The holdout re-derives the α-dependent ledger entries as committed;
  it is an integration test of mutual consistency in α, not a
  re-derivation of every probe.

## 8. What would falsify this

- A canonical equation fixing R_f — the modulus claim would be wrong.
  (Checked: exact rescale invariance; V_tree = 0 symbolically.)
- A charge normalization with residual freedom — flux quantization or
  the flow slope would drift. (Checked: 10⁻¹¹ periodicity; 3×10⁻⁴
  slope; the evolved force check.)
- An Einstein-frame reduction with a different dilaton exponent — the
  α(φ) law would decouple from the geometric law. (Checked: 3b = −6a
  = −√3 symbolically, equal and opposite to the geometric exponent.)
- A committed constraint hitting one of the "flat" directions with an
  α-coupled component — the rank audit would show a nonzero residual
  projection. (Checked: 2×10⁻¹⁶.)
- A lepton identified with a fiber KK level — the 5×10¹⁸ hierarchy
  would be a contradiction rather than a decoupling.
- Two committed α-dependent constants implying *different* α's — the
  holdout inversion would split. (Checked: 4×10⁻¹⁶ relative spread
  across twelve constants.)

## 9. Companion probe

`experiments/closure_ledger/absolute_coupling_capstone_probe.py`
(T1–T9, ~5 s; run in CI by
`tests/test_absolute_coupling_capstone.py`): the canonical chain, the
#193 weld, and the geometry map; spectral flow, flux quantization,
and the adiabatic-ramp force check; the symbolic Einstein-frame
reduction; the α law, the modulus, and the rank audit; the stabilizer,
the guardrail, and the lepton-sector compatibility; the α-dependent
holdout over the committed ledgers.

**Verdict:**
`CANONICAL_ALPHA_EQUALS_4K2_OVER_RHO2_THE_EINSTEIN_FRAME_RADION_IS_ITS_E_SQRT3_PHI_DILATON_THE_EM_CAP_FIXES_THE_ONLY_ALPHA_COUPLED_DIRECTION_AND_THE_ALPHA_DEPENDENT_HOLDOUT_PASSES_ON_ONE_COMMON_ALPHA`

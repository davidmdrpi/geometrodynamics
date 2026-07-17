# The soliton as perturbative dressing, and the two-mouth network port (PR #223)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. #222 proved the frozen primordial
> background is the only consistent reading of the anchor; this PR
> returns to it with the two requested steps: the scalar as a **light
> perturbative dressing** of the fixed-μ bridge, and the **port of the
> coupled solver to the full two-mouth wormhole topology**, mapping the
> scale-welded soliton under global, non-local network transits. The
> companion probe machine-checks every claim (~3 min).

## 0. The geometry

The **ultrastatic MTY bridge** (the frozen-transit reading): spatial
geometry = the #202 t = 0 bridge, ρ(σ) = √(r_s² + σ²), lapse = 1
(traversable, no horizon). The exact reduction u = ρ^{3/2}φ gives

```
u'' + [ω² − V_b] u = 0,
V_b = ℓ(ℓ+2)/ρ² + (3/2)ρ''/ρ + (3/4)(ρ'/ρ)²
```

with the neck height **(ℓ(ℓ+2) + 3/2)/r_s²** in closed form and the far
tail **(ℓ(ℓ+2) + 3/4)/σ² — exactly the #215 form** (machine-checked to
10⁻³): the 5D ρ³ measure is carried exactly by the reduction, answering
the audit's measure item on the true geometry. The network ring: neck →
bridge → mouth → S³ exterior arc (π R_u to the antipodal source point),
with the S³ zonal reduction u = R sin χ·ψ exact on the arcs; #202's Pin
parity selects the electron (k = 1) mode **odd at the neck** (the
near-neck φ ∝ σ law reproduced, drift 2%).

## 1. The port: the spectrum and the network comb

The two-mouth ring solved by shooting from the neck with both parities
(neck and source-point): the odd-neck (electron-parity) modes of the
two source parities interleave, and the same-parity spacing is the
**network closure comb π/L_half** to 2% — the #217 resonance-comb
structure realized on the true bridge geometry.

## 2. The new universal: X = 2.2995, in closed form

On the true bridge (no interior channel), the #202 matching radius —
the first antinode of the *field* φ = u/ρ^{3/2} beyond the neck — gives

```
X = σ_match · ω = 2.2995
```

**invariant** across modes (three modes agree to 4 digits), r_s (the
anchor limit r_s·ω → 0), exterior curvature (flat vs S³ arcs
identical), and the mouth seam. And it has a **closed form**: near the
neck the odd solution is u = √σ·J₂(ωσ), so φ = J₂(z)/σ and the first
antinode solves

```
z·J₁(z) = 3·J₂(z)      →      z* = 2.29948…
```

— the 5D two-mouth analog of the quarter wave (#221's π/2 was the
1D-cavity case; z* is its Bessel-index-2 counterpart on the true
bridge).

## 3. The interior-depth family: #221 and #223 unified — and a structural bound

Inserting a flat interior channel of depth D at the neck (the
#215 horizon-tortoise reading, which #216–#221 modeled with regulated
depth) gives a sharp step:

| D | X_match |
|---|---|
| 0 (the pure bridge) | 2.2978 |
| 2, 4, 8, 12 (any channel past the quarter wave) | **π/2 exactly** (to 10⁻⁴) |

The #221 geometry and the true bridge are **one family**, discriminated
by the interior depth — and across the *whole* family

```
X ≥ π/2      ⟹      m_e/m_μ = α/X ≤ 2α/π = 0.004646
```

vs observed 0.004836: the required conv-B value 1.5089 sits **3.9%
below the class infimum**. The residual is **structural**, not a
modeling choice: either the EM-cap anchor r_s = α·λ̄_C carries a ~4%
correction, or the #201 S₁ = m_μ convention does. (The #165 guardrail
forbids matching the 3.9% to any constant without a derivation —
deriving it is the named successor, alongside the physical interior
depth that adjudicates between 2.2995 and π/2.)

## 4. The network non-locality map: transit protection

The even/odd frequency splitting through the neck — the through-bulk
identity exchange of the two-mouth network acting on the mode — obeys
the Bessel-index tunneling law:

```
Δω  ∝  (ω·r_s)^p,       p = 4.15 / 4.04   (predicted 4 = 2ν)
```

Extrapolated to the primordial anchor (r_s·ω = α): the electron mode's
non-local exchange is **O(α⁴) ~ 10⁻¹²** of its frequency. **The dressed
soliton is transit-protected**: it sits still on the network while the
even/ℓ = 0 carrier channels — the transactional arc's waves
(#213–#221) — do the transiting. The scale-welded soliton's behavior
under global network transits is: protected in identity, coupled only
through the carrier comb.

## 5. The perturbative dressing, quantified

The first-order back-reaction — the linearized 5D mass function
δμ(σ) = (2κ/3)∫ρ³ρ_E dσ from the mode's stress on the fixed-μ bridge:

- **exactly quadratic in the amplitude** (δμ(A)/δμ(2A) = 1/4 to
  machine precision) — the dressing is genuinely perturbative;
- **the throat-local share is 0.2%** of the cloud total
  (δμ(σ < 3r_s)/μ = 0.002 at A = 0.05): the neck geometry is
  undisturbed — **μ remains the primordial datum**, exactly the #222
  division of roles;
- the cloud energy is the particle's contribution to the exterior mass
  (the exterior Tangherlini mass beyond the cloud reads μ + E_cloud);
- the perturbative window: A ≈ 0.014 for a 10% μ-shift — the "light
  dressing" regime quantified.

## 6. Honest scope

- The ultrastatic lapse is the MTY frozen-transit reading; the #215
  horizon reading carries the tortoise interior instead. The two are
  now one family bracketing X (2.2995 vs π/2); deriving the physical
  interior depth — or adjudicating the readings — is the named
  successor.
- The 3.9% one-sided residual is structural; matching it to a constant
  without derivation is forbidden (#165).
- The zonal 1D reduction is exact in the measures, but the brane/bulk
  mouth seam is modeled as transparent (u, u′) continuity
  (seam-position robustness checked).
- The back-reaction is first-order; second-order geometry shifts are
  (δμ/μ)² and unresolved.
- Classical, frozen background; α enters only through the anchor.

## 7. What would falsify this

- The bridge potential missing the #215 far form or the #202 near-neck
  law — the reduction would be wrong. (Checked: 3.75/σ² to 10⁻³; φ ∝ σ
  to 2%.)
- X depending on the mode, r_s, curvature, or seam — no universal.
  (Checked: 4-digit invariance; the closed form z·J₁ = 3J₂ hit to
  10⁻³.)
- The interior channel failing to reproduce π/2 — #221 and #223 would
  be inconsistent models. (Checked: π/2 to 10⁻⁴ for every D past the
  quarter wave.)
- A splitting exponent far from 2ν = 4 — the transit coupling would
  not be the Bessel tunneling. (Checked: 4.15/4.04.)
- δμ non-quadratic in A, or throat-dominant — the dressing would not
  be perturbative, contradicting #222's forced reading. (Checked:
  exact quadratic; 0.2% local share.)

## 8. Companion probe

`experiments/closure_ledger/bridge_dressing_network_probe.py` (T1–T9,
~3 min): the bridge closed forms and welds; the two-mouth spectrum and
comb; the universal X with its Bessel closed form; the interior-depth
family and the structural bound; the transit-protection law; the
quantified perturbative dressing.

**Verdict:**
`THE_TWO_MOUTH_PORT_YIELDS_A_CLOSED_FORM_UNIVERSAL_ZJ1_EQUALS_3J2_THE_INTERIOR_DEPTH_STEPS_IT_TO_PI_OVER_2_THE_ELECTRON_MODE_IS_TRANSIT_PROTECTED_AS_ALPHA_TO_THE_4_AND_THE_DRESSING_IS_EXACTLY_PERTURBATIVE`

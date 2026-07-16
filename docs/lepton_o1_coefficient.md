# The derived O(1) lepton mass coefficient (PR #221)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. #210 closed the mass-ladder thread
> onto one number: the O(1) coefficient σ_mode/λ̄_C, required
> 88.6α = 0.6467 (convention A) or 206.8α = 1.5089 (convention B) —
> "constrained, not derived", a ×2.3 band. This PR **derives it** —
> no fit anywhere — on the spatially converged, spectrally stable PDE
> eigenhistory background of #220. The companion probe machine-checks
> every claim (~2 min).

## 0. Why the coefficient is derivable at all

The eigenhistory's **mass is its frequency** (m = ℏω*), so
λ̄_C = 1/ω*. And σ_mode is a length on the *same* object. The
coefficient

```
X  =  σ_mode / λ̄_C  =  σ_mode · ω*
```

is therefore a ratio of two properties of **one geometric mode** —
there is nothing left to fit. Everything reduces to: *which mode, and
which length?* Both answers are repo canon, fixed before this PR:

1. **Which mode.** The transactional arc (#217/#218) puts completed
   eigenhistories on the throat's **interior resonance comb**; the
   ring interior of the #220 background is the bridge coordinate of
   #202's 5D throat solve, with the neck (interior midpoint) as the
   cross-cap. #202's machine-verified **Pin-twisted boundary
   condition** — odd-k modes have a node at the cross-cap — forces the
   electron (k = 1) transit mode to be the **odd interior
   fundamental**. (The even mode, which touches the neck, is the k = 0
   uncharged channel — #202's dichotomy, realized here: the measured
   odd mode has its node at the neck to 10⁻¹⁵, and its near-neck
   profile is exactly #202's regular solution φ ∝ σ.)
2. **Which length.** #202's suppression law ε₁ = r_s/σ_mode defines
   σ_mode as the **matching radius** — where the interior linear
   growth φ ∝ σ turns over into the mode cloud: on the mode, the
   **antinode distance from the neck**.

## 1. The hard-wall theorem: X = π/2 exactly

For the odd fundamental of a cavity the antinode sits at the quarter
wave: σ_match = L_eff/4 with ω = 2π/L_eff, so

```
X = ω · σ_match = π/2      (exactly, for ANY cavity length)
```

— the analytic reason the coefficient is O(1) and universal: the
cavity length cancels. The definitional family is closed-form on hard
walls (all machine-verified by sampled quadrature to 10⁻⁵):

| definition | even (k = 0 channel) | odd (k = 1, the electron) |
|---|---|---|
| amplitude RMS (the #203 convention) | π·√(1/12 − 1/2π²) = 0.5679 | √(π²/3 − ½) = 1.6703 |
| **matching radius (the #202 convention)** | — (antinode at neck) | **π/2 = 1.5708** |
| energy RMS | π/√12 = 0.9069 | 2π/√12 = 1.8138 |

## 2. The measurement on the real Tangherlini barriers

On the #220 background (the glued finite-width greybody barriers of
#215; S³ exterior arc 2π; reference interior depth D = 12, where the
odd mode is well-trapped, interior fraction 0.95):

- **X_match = 1.579** (reference), converged in dx to < 10⁻³;
- **regulator scan** (the interior depth is the #215
  horizon-continuum regulator): depth 8 → 48 gives
  X_match ∈ [1.579, 1.638] — a 4% band, all within 6% of π/2;
- **exterior independence**: S³ arc π → 4π moves X_match by < 0.04 —
  and the matching radius stays clean exactly where the RMS definition
  is contaminated by interior–exterior hybridization (it is a local
  feature of the profile, not a second moment);
- the parity dichotomy: the even mode's antinode is *at* the neck
  (u_neck/u_max = 1), the odd mode's node is at the neck to 10⁻¹⁵.

## 3. The eigenhistory carries it

The full #220 Gauss–Newton machinery (periodic orbit of the complete
source–field state, energy closure H = E₀, phase condition) seeded on
the odd interior fundamental:

- converged to residual **3×10⁻¹³**;
- **the source decouples exactly** (q* = p* ≈ 10⁻¹⁷): the odd/charged
  transit mode has u(0) = 0 by parity — it is *source-transparent at
  the antipodal crossing* (the phase condition becomes ⟨u, π⟩ = 0);
- the complete monodromy — source in the tangent space — is on the
  **unit circle to 10⁻¹⁴**: the coefficient lives on a spectrally
  stable object;
- the orbit-averaged profile carries **the same X_match** as the
  linear mode (to 10⁻³), and the decoupled sector is exactly linear:
  X is **amplitude-independent** (the half-amplitude state is
  periodic at the same T to 10⁻¹⁰).

## 4. The confrontation — zero fitted numbers

Inputs: the throat geometry (r_h = 1) and α (the #184-protected
boundary invariant), through #210's primordial anchor r_s = α·λ̄_C and
#202's exact k = 1 law ε₁ = r_s/σ_mode. Then

```
m_e/m_mu = alpha / X
```

| | required X | derived X | verdict |
|---|---|---|---|
| convention A (S₁ = (3/7)·m_μ) | 0.6467 | 1.579–1.638 | **excluded ×2.4** |
| convention B (S₁ = m_μ) | 1.5089 | 1.579–1.638 | **selected** |

```
m_e/m_mu = alpha/X = (2 alpha/pi)·(1 + throat shift)
         = 0.004455 – 0.004622      (hard-wall anchor 2alpha/pi = 0.004646)
observed = 0.004836
```

**The derivation lands at 92–96% of the observed ratio.** The #201
convention ambiguity is resolved *by the measurement* (B selected, A
excluded), the #210 O(1) band ×2.3 collapses to a derived value with a
4% regulator band, and the neck aspect becomes
c = ln(X/α) = 5.38 vs the convention-B-observed ln(m_μ/m_e) = 5.33 —
#210's "c = ln(1/α) + O(1)" with the O(1) now **computed** (ln X ≈
ln π/2 = 0.45).

The residual 4–8% is real and one-sided (the throat correction pushes
X *up* from π/2 while the observation sits 4% *below* it) — owned by
the named successor below.

## 5. Honest scope

- The ring-transit ↔ #202-bridge identification (ring interior = the
  bridge σ coordinate, neck = cross-cap) is structural and
  machine-consistent here (exact node at the neck; the φ ∝ σ near-neck
  law), but the 1D transit measure is not the 5D ρ³ bridge measure —
  **redoing X on the true 5D radial operator is the named successor**
  and the leading candidate for the 4–8% residual.
- The winding k is not represented in the zonal scalar; the k = 1
  parity is imported from #202's theorem, not re-derived.
- m_μ enters only as the anchor scale of the #201 law; **m_e is used
  nowhere in the construction** — only in the final comparison.
- The interior depth is the #215 horizon-continuum regulator; X_match
  is stable to 4% over depth 8–48 (band carried); the physically
  capped depth (the α-capped tortoise depth) is a successor
  refinement.
- The even-parity/conv-A alternate lands in-band numerically but is
  **parity-excluded** by #202; it is recorded, not used.
- Classical, zonal scalar, frozen geometry; ℏ enters only through
  λ̄_C = ℏ/mc (the B4 anchor), as always.

## 6. What would falsify this

- An even-parity (node-free) mode at the neck for k = 1 — the #202
  boundary condition would fail on the transit background. (Checked:
  node to 10⁻¹⁵.)
- X depending on the cavity length, the regulator depth, or the
  exterior arc — the coefficient would not be an O(1) of the throat.
  (Checked: hard-wall exact cancellation; 4% regulator band; < 0.04
  exterior spread.)
- X landing outside both conventions' requirements — the primordial
  anchor r_s = α·λ̄_C would be wrong. (It lands 5–8% above conv B and
  ×2.4 above conv A: the anchor survives, and the convention is
  resolved.)
- The eigenhistory orbit shifting X or going spectrally unstable.
  (Checked: same X to 10⁻³; monodromy unit-circle to 10⁻¹⁴.)

## 7. Companion probe

`experiments/closure_ledger/lepton_o1_coefficient_probe.py` (T1–T8,
~2 min): the hard-wall closed forms; the parity dichotomy and the
measurement; the regulator/exterior/grid robustness; the
source-decoupled eigenhistory orbit with complete monodromy; the
confrontation.

**Verdict:**
`THE_REMAINING_O1_IS_DERIVED_THE_ODD_THROAT_FUNDAMENTAL_GIVES_X_EQUALS_PI_OVER_2_PLUS_THROAT_SHIFT_AND_ME_OVER_MMU_EQUALS_ALPHA_OVER_X_AT_92_TO_96_PERCENT_NO_FIT`

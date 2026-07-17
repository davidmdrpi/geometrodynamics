# The derived O(1) lepton mass coefficient (PR #221)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. #210 closed the mass-ladder thread
> onto one number: the O(1) coefficient σ_mode/λ̄_C, required
> 88.6α = 0.6467 (convention A) or 206.8α = 1.5089 (convention B) —
> "constrained, not derived", a ×2.3 band. This PR derives the
> **quarter-wave cavity invariant** — no fit anywhere — on the
> spatially converged, spectrally stable PDE eigenhistory background
> of #220, demonstrates it invariant under energy budget, cavity
> depth, source coupling, mode amplitude, and branch choice, and
> lands m_e/m_μ at 92–96% **conditional on the soliton↔cavity weld**
> (§5b: the independent identifiability audit that narrows the claim).
> The companion probe machine-checks every claim (~2 min).

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

## 4. The invariance suite

The coefficient is demonstrated invariant under every knob of the
construction:

| knob | range | result |
|---|---|---|
| **energy budget E₀** | ×16 (0.178 → 2.844), *independent* Gauss–Newton solves | X and T identical to 10⁻⁷/10⁻⁹ |
| **source coupling g** | 0 → 4× the #220 value | the *same* orbit periodic to 10⁻¹² at every g (u(0) = 0 by parity); complete monodromy at 4g unit-circle to 10⁻¹⁴ |
| **mode amplitude** | ×¼ → ×2 | one knob with the budget for the exactly linear decoupled sector; the half-amplitude state periodic at the same T to 10⁻¹³ |
| **cavity depth** | 8 → 48 (the #215 regulator) | X_match band 4% (§2) |
| **branch choice** | the first three odd interior branches | X_match = 1.620 / 1.575 / 1.571 → **π/2**, spread < 0.05 |

The branch row is the sharpest: each branch's first antinode is a
quarter of *its own* wavelength, so X_match is branch-invariant
(hard-wall exactly π/2 for every branch, and the throat values
converge *toward* π/2 up the ladder as the barriers look harder to
shorter wavelengths) — while the amplitude-RMS definition **grows
×3.3** from branch 1 to branch 3. Branch invariance singles out the
#202 matching radius as the physical definition of σ_mode; the RMS
convention is a fundamental-only surrogate.

**The alternative branch is run, not just measured**: the second odd
branch's complete Gauss–Newton orbit (deep cavity D = 24; residual
6×10⁻¹³; the source still exactly decoupled, q* ≈ 10⁻¹⁵; complete
monodromy unit-circle to 2×10⁻¹³) carries

```
X_match(branch 2 orbit) = 1.57070   —   π/2 to 1×10⁻⁴
```

— the alternative rung of the eigenhistory ladder is a genuine,
spectrally stable, energy-closed orbit carrying the *same* derived
coefficient, closer to the hard-wall theorem than the fundamental.

## 5. The confrontation — zero fitted numbers

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

**The alternate branch, run.** The even (k = 0 channel) branch behind
the conv-A alternate is put through the *same* full machinery — the
source-**coupled** Gauss–Newton orbit (residual 4×10⁻¹³, q* =
−6.2×10⁻³ genuinely nonzero, the coupling contrast with the odd
branch; complete monodromy unit-circle to 3×10⁻¹²) carries
X_u = 0.6395, landing the parity-excluded alternate at

```
m_e/m_μ (even + conv A) = 0.004891 = 101.1% of observed
```

— **numerically closer than the primary**. Stated plainly: the
alternate's exclusion is purely structural (#202's parity theorem
forces the k = 1 mode odd), not numerical. The parity identification
is therefore the sharpest falsification target for the 5D
bridge-measure successor, which adjudicates between the two readings
— either outcome is decisive.

## 5b. The independent identifiability audit — the claim, narrowed

`lepton_o1_identifiability_audit_probe` (independent companion, 8/8
PASS) audits the *logical status* of the coefficient, and its verdict
is adopted here:

- **Unconditional**: the quarter-wave invariant — X_cavity = π/2 +
  throat shift, regulator/exterior/amplitude/coupling/branch-robust —
  is a validated property of the throat cavity.
- **Conditional**: the m_e/m_μ landing. Within the *existing*
  equations, #202/#203's σ_mode is the **soliton** IR localization
  scale (the vacuum throat problem has no bound state; it supplies
  only the suppression law), and no current equation welds the soliton
  length unit to the cavity/bridge unit: the audit shows the soliton
  length family (r50/r90/r95/RMS/R*) spans **×4** under the cavity
  frequency, and a rescaling of the soliton radial unit moves X
  linearly while leaving all normalized soliton shape data unchanged.
  Identifying σ_mode with the cavity antinode — the
  eigenhistory-particle identification this PR is built on — is
  therefore a *posit* the current equations neither derive nor
  exclude.
- **The successor contract** (the audit's executable list): one
  coupled Pin-Dirac/soliton state on the Tangherlini bridge; the ρ³
  5D radial measure in the normalization; the 3D-soliton-radius →
  5D-bridge-coordinate map derived from the action; the
  neck-to-asymptotic overlap functional computed, not defined by an
  antinode; convergence under grid/boundary/matching refinement;
  Floquet stability used only to reject unstable branches; the
  coefficient locked *before* comparison with the observed masses.

The residual 4–8% is real and one-sided (the throat correction pushes
X *up* from π/2 while the observation sits 4% *below* it) — owned by
the named successor below.

## 6. Honest scope

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
- The even-parity/conv-A alternate is run through the full machinery
  (§5) and lands at 101.1% of observed — *closer than the primary* —
  but is **parity-excluded** by #202; the exclusion is structural, not
  numerical, and the 5D bridge-measure successor adjudicates.
- Classical, zonal scalar, frozen geometry; ℏ enters only through
  λ̄_C = ℏ/mc (the B4 anchor), as always.

## 7. What would falsify this

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
- X depending on the energy budget, the source coupling, the mode
  amplitude, or the branch — the coefficient would be a property of
  the *state*, not of the throat geometry. (Checked: invariant under
  all four — §4; the budget ×16 through independent solves moves X by
  < 10⁻⁶, and the branch ladder returns π/2.)

## 8. Companion probe

`experiments/closure_ledger/lepton_o1_coefficient_probe.py` (T1–T9,
~2 min): the hard-wall closed forms; the parity dichotomy and the
measurement; the regulator/exterior/grid robustness; the
source-decoupled eigenhistory orbit with complete monodromy; the
five-knob invariance suite with the alternative branches run; the
confrontation with the conditionality stated.

`experiments/closure_ledger/lepton_o1_identifiability_audit_probe.py`
(independent, 8/8): the definition ledger; the quarter-wave benchmark
classification; the soliton localization family; the cross-sector
non-uniqueness; the scale-weld falsification; the stability role; the
successor contract.

**Verdict:**
`THE_QUARTER_WAVE_INVARIANT_IS_DERIVED_PI_OVER_2_PLUS_THROAT_SHIFT_AND_THE_ME_OVER_MMU_LANDING_AT_92_TO_96_PERCENT_IS_CONDITIONAL_ON_THE_SOLITON_CAVITY_WELD`

**Independent audit:**
`lepton_o1_identifiability_audit_probe` (8/8) —
`QUARTER_WAVE_INVARIANT_VALIDATED_BUT_LEPTON_COEFFICIENT_NOT_IDENTIFIABLE_UNTIL_THE_SOLITON_TO_5D_SCALE_WELD_IS_DERIVED`

# The 5D throat-core solve: the exact suppression law on the Tangherlini bridge (PR #202)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. This PR executes the next #200
> register item (via #201's closing sentence): solve the k-winding mode
> problem on the actual 5D throat geometry, so that #201's fitted mouth
> coupling becomes an exact geometric law. The companion probe verifies
> every claim to stated precision.

## 0. Results, up front

On the t = 0 slice of the J-quotiented Tangherlini throat (#167–#169),
the k-winding mode problem reduces to a 1D problem on the bridge
ρ(σ) = √(r_s² + σ²) (σ = proper distance from the neck; the coordinate
identity dσ = dρ/√f verified to machine precision), with three exact
results:

1. **The Pin-twisted boundary condition.** The deck map acts as
   (σ → −σ) × (fiber phase (−1)^k), so the quotient projects odd-k
   modes onto φ(0) = 0 — **fermionic (Pin-twisted) modes have a node at
   the cross-cap** — while even-k modes obey the Neumann condition and
   remain O(1) there. Machine-checked: the half-bridge Dirichlet
   (Neumann) spectrum equals the full-bridge odd-parity (even-parity)
   spectrum to 10⁻⁶. This is the **geometric realization of #195's
   forbidden one-mouth mass lift**: the self-pairing channel dies at
   the identification locus for exactly the odd k that carry the
   fermions.
2. **The exact radial law.** At E ≈ 0 (the physical regime — bound
   energies are tiny on throat scales) the regular solution behaves as
   φ ∝ σ^k far from the neck; **for k = 1 the regular mode is exactly
   φ(σ) = σ** on the whole bridge (closed form:
   (ρ³σ′)′ = 3ρσ identically). Far-field exponents verified to 10⁻⁴;
   a genuinely shallow bound state's interior tail matches the regular
   solution pointwise to ~2% (the tail theorem).
3. **The suppression law is a power law in the scale ratio.** The
   amplitude suppression between a mode of extent σ_mode and the neck is

   ```
   ε_k  =  (r_s / σ_mode)^k · e^{c₀(k)} ,
   c₀ = { 0 (exact),  −0.405,  −0.783 }  for k = 1, 3, 5 ,
   ```

   verified: d ln ε_k / d ln(σ_mode) = −k to four decimals. This is the
   5D derivation of #201's ε_k = e^{−kc}: **c = ln(σ_mode/r_s)** — the
   "neck aspect" is the logarithm of the throat-to-mode scale ratio.

## 1. What this does to the naturalness ledger

#201 reported the rebuilt electron's worst Barbieri–Giudice sensitivity
as Δ = c ≈ 4.48. The 5D law shows that number was an artifact of
parametrizing by the *exponent*: in the physical parameter — the scale
ratio x = σ_mode/r_s — the electron mass is **linear**,

```
m_e ∝ (r_s/σ_mode)¹    ⟹    |d ln m_e / d ln x| = k = 1  exactly.
```

The electron level's sensitivity to the actual geometric quantity is
**1 — maximal naturalness**. The full rebuilt-ladder table (re-derived
in the probe): every physical-parameter sensitivity ≤ 1 for the
electron chain, heavy sector unchanged (≤ 0.89, #201).

## 2. The inversion, now exact

With c₀(1) = 0 the inversion has no O(1) fudge: ε₁ = 0.01129
(convention A, #201) gives

```
σ_mode / r_s  =  1/ε₁  =  88.6      (convention B: 206.8) ,
```

i.e. **the throat's GR core radius is ~0.5–1% of the particle's wave
extent** — the electron/muon mass ratio is, structurally, a direct
measurement of the throat-to-mode scale hierarchy, entering linearly.
This is physically the right shape (a point-like core inside a
Compton-scale cloud); its *microphysical* value is regime-dependent and
its independent determination — the coupled 5D + #180-soliton solve —
is the one number that still separates this from an outright
m_e/m_μ prediction. The unknown has been reduced to **one dimensionless
ratio governed by an exact law**.

## 3. The impossibility bound, with exact exponents

ε₃/ε₁ = (r_s/σ_mode)² e^{−0.405} ≈ 8.5×10⁻⁵ — versus the ε₃/ε₁ ≈ 88.6
that a pairing-generated μ/e would require (#201). The 5D exponents
*steepen* the suppression with k (power k, plus the neck constants), so
the conclusion hardens: **the inter-generation hierarchy cannot come
from mouth pairing at any scale ratio**; it stays with the dynamical
uplift (fitted and natural, #194).

## 4. The zero-winding contrast

The parity condition splits by k (odd: Dirichlet node; even: Neumann),
and *every* winding sector k ≥ 1 is additionally barrier-suppressed as
σ^k. The clean contrast is the **zero-winding channel**: for k = 0 the
E ≈ 0 Neumann solution is exactly flat ((ρ³φ′)′ = 0 admits φ = const),
so the windingless (vacuum/boson) mode keeps an O(1) amplitude at the
neck while the k = 1 (electron) mode is suppressed to ~(r_s/σ_mode) —
measured contrast ×40 at σ_mode = 40 with equally shallow wells. Only
the uncharged channel touches the identification locus: the
multiplicative mass protection is specific to the charged/winding
sectors, with the odd (Pin-twisted) ones carrying the node on top — the
#183/#193 grading with its dynamical face on.

## 5. Honest scope

- Rigid vacuum Tangherlini background (no back-reaction); the E → 0
  regime (verified appropriate: the shallow-well test at
  |E| ≪ k(k+2)/σ_mode²).
- σ_mode (the IR localization) is supplied by the soliton sector
  (#180/#185), not by the vacuum geometry — as it must be (the vacuum
  problem has no bound state; the winding barrier is repulsive).
- The scalar reduction of the winding sector (the spinor case adds the
  #195/#197 structure to the angular part; the radial barrier and the
  parity boundary condition are as computed here).
- The pairing constant between two geons at brane separation R
  (#185/#201's kernel) is the *exterior* overlap; this PR's law governs
  the interior/neck side and the k-scaling. The two match at the mode
  boundary — the matching constants c₀(k) are exactly the O(1) numbers
  computed.

## Reproduce

```bash
python -m experiments.closure_ledger.throat_core_5d_solve_probe
# Verdict: FIVE_D_THROAT_CORE_SOLVED_ODD_K_NODE_AT_THE_CROSS_CAP
#          _SUPPRESSION_IS_A_POWER_LAW_IN_THE_SCALE_RATIO_ELECTRON_SENSITIVITY_ONE
```

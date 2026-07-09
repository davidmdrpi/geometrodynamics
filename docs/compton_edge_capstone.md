# The Compton-edge capstone: Release II (PR #211)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. This is the closing PR of the
> #204–#210 arc: the whole chain re-verified green in one run (as #200
> did for the previous arc), the register updated, and one new
> parameter-free law that sharpens the arc's final unknown. The
> companion probe checks every claim (~1 min).

## 0. What the arc did, and what closes it

The #204–#210 arc walked the program's two named frontiers to their
edges and closed both:

- **the nonlinear measurement theory** (#198's condition 2): audited for
  no-signaling (#204), the committed Φ[ρ] confronted with the laboratory
  (#205), configuration space derived from the bridge (#206), swapping
  from surgery (#207), GHZ from valence (#208), measurement made
  operational (#209);
- **the strong-field core contraction** (#203's pre-registered target):
  measured and adjudicated (#210) — the collapse reading refuted
  constructively, the mass anchor relocated to α.

This capstone consolidates all of it and carries **one new result**
aimed at the single number the arc left open.

## 1. The Compton-edge law (new)

#210 left the mass ladder with one O(1): σ_mode/λ̄_C. Rather than *dial*
that number (which would undo what twenty PRs removed), this PR
**bounds** it with a parameter-free computation. Add the mass term to
#202's bridge equation:

```
(ρ³φ′)′ = [ k(k+2)ρ + m²ρ³ ] φ ,     ρ(σ) = √(r_s² + σ²) ,   k = 1 .
```

Measured consequences:

1. **The massless limit re-derives #202's exact law**: φ(σ) = σ to
   machine zero (|φ − σ|/σ = 0), i.e. c₀(1) = 0 exactly — the
   suppression law's anchor, recovered.
2. **Universality**: the deformation depends on the single scaling
   variable **x = σ_mode/λ̄_C** alone. The m·r_s = 0.02 and 0.05 curves
   coincide to 2.6×10⁻⁴ across the whole grid — one law, not a family.
3. **The law**: ε₁ = (r_s/σ_mode)·D(x), with the deformation D(x)
   computed (D = 0.966 at x = 0.647, 0.831 at x = 1.509, falling
   through 0.50 at x = 3).
4. **The sensitivity window**:
   S(x) = |d ln ε₁ / d ln σ_mode| is **exactly 1 below the Compton
   scale** — the #202 naturalness result, now with its validity domain
   located — and grows beyond it:

| x = σ_mode/λ̄_C | 0.25 | 0.647 | 1.0 | 1.509 | 2.0 | 3.0 | 5.0 |
|---|---:|---:|---:|---:|---:|---:|---:|
| S(x) | 1.01 | 1.07 | 1.16 | 1.36 | 1.62 | 2.28 | 3.95 |
| D(x) | 0.995 | 0.966 | 0.921 | 0.831 | 0.726 | 0.501 | 0.179 |

   **Naturalness caps the O(1)**: S ≤ 2 confines x ≤ **2.58**, and the
   ladder's own worst-tolerated sensitivity (4.48, from #201/#202)
   confines x ≤ **5.58**. Both convention anchors sit *inside* the
   natural window (S = 1.07 at conv A, 1.36 at conv B).

5. **The tightened band**: solving the deformation-corrected
   self-consistency x = (needed·α)·D(x) gives
   **σ_mode/λ̄_C ∈ [0.626, 1.307]** — tighter than #210's [0.647, 1.509],
   and **still bracketing 1**. The successor derivation now has a
   compact, derived target window at the Compton edge of the natural
   domain, rather than an unbounded O(1).

The reading is physically transparent: the winding mode is a
Compton-scale cloud (σ_mode ≈ λ̄_C) around a throat whose GR core is the
classical electron radius (r_s ≈ α·λ̄_C, #55–#58), and the suppression
exponent is c = ln(1/α) + O(1) — the electron/muon ratio as a
fine-structure phenomenon, with the correction to the exact law
controlled and small precisely where the anchor sits.

## 2. The ledger, re-verified green in one run

**Ledger 1 — the commitment chain (#204/#205).** A kicked, *retarded*
evolution with real potentials conserves the norm to 6×10⁻¹⁴ and closes
the continuity equation at 2×10⁻³ (integrator error) — no-signaling's
causal completion costs the Born rule nothing. The classical mean-field
channel keeps a product state at entanglement entropy 0 while the
quantized pairwise operator entangles it (S = 0.026) — the BMV-null
stands. The lab card: SN phase at the interferometry record 4.8×10⁻¹⁷
rad, the record mass ~10⁵ below the SN inhibition scale, the BMV witness
phase 0.79 rad where BAM predicts strictly zero.

**Ledger 2 — the topology chain (#206–#208).** ψ_eff = (I⊗T)|Φ⁺⟩ is the
singlet (fidelity 1.0) with E(a,b) = −cos(a−b) to 2×10⁻¹⁶; the swapping
law leaves (1,4) in Ψ_{a+b+c} (fidelity 1.0, probability ¼) with a
separable unconditioned holonomy mixture (negativity 2×10⁻¹⁷); the
enumerated local bounds CHSH = 2 and Mermin = 2 against the bridge
sector's 2√2 and the Y-junction's Mermin = 4, whose pairwise marginals
are exactly unentangled (negativity 0). Bridge valence is the
entanglement class.

**Ledger 3 — measurement and mass (#209/#210).** The k-odd dispersion
identity holds to 2×10⁻¹⁶ (the pointer coupling is the KK connection
gradient); a live Stern–Gerlach Born check lands within 0.0024 of
cos²β; a live Kaup point gives M = 0.632 (the strong-field solver that
refuted the collapse reading); and the α chain re-verifies —
σ_mode/λ̄_C ∈ [0.647, 1.509] brackets 1, and m_e/m_μ = (3/7…1)·α ∈
[0.00313, 0.00730] brackets the observed 0.00484.

## 3. The register after the arc

**Open items (all named and bounded):**
1. Derive σ_mode/λ̄_C within the §1 Compton-edge window (≤ 2.58 natural
   / ≤ 5.58 tolerated) — a *bounded* target, no longer an open number.
2. Connect the #55–#58 equilibrium radius R* to the 5D bulk mass
   μ = r_s².
3. The 5D pants (trousers) nucleation — its Sorkin causal-continuity
   class and rate (the GHZ carrier, #208).
4. W-class reachability via bridge networks + surgery (#208).
5. Registration/irreversibility beyond branch separation (radiative
   decoherence; the #204 dissipative machinery).

**Standing negative kept on the books:** the cosmological constant /
one-R failure (#165).

**The falsification card:** two near-term nulls — SN-scale signatures
(a detection refutes the committed sourcing) and a BMV entanglement
witness (a detection refutes classical Φ outright) — plus the
m_e/m_μ = (3/7…1)·α prediction at the ×1.5 level with its O(1) now
confined, and the unchanged neutrino cards (m_ββ ≈ 1.5–6 meV,
Σm_ν ≈ 59 meV).

## 4. What the arc did *not* do (honest scope)

- Full-GR no-signaling remains a completion argument (the wave-equation
  Φ is the minimal causal completion; the constraint analysis is not
  run).
- The emergence of the configuration-space pilot wave from the 3-space
  field is derived only in its gravitational sourcing (#205) and its
  entangled structure (#206–#208); the full field-level account remains
  a program.
- The pinch and the pants are represented by their topological content;
  the 5D dynamics (Sorkin classes, rates) is open.
- The α anchor is **constrained** (the Compton-edge window), not
  derived; conv A vs B is a stated convention band; the #55–#58 R*
  itself carries the program's one dimensionful anchor (per B4/#184).
- Equilibrium is a hypothesis with a mechanism (the #198 relaxation)
  throughout.

## 5. Assessment — Release II

When the arc opened, the program's deepest interpretive machinery rested
on two unexecuted items: a measurement theory that might not survive
Gisin, and a mass anchor that might not survive GR. Both edges were
walked to, and both fired exactly where the program had pre-registered
them — the Newtonian model signals (and the causal completion fixes it
for free); the cores collapse (and collapse cannot make them light).
What emerged is stronger than what was risked: Bell violation as the
pointer statistics of classical beables with every ingredient derived —
the tensor product from the bridge, the singlet from the Pin⁻ transport,
the outcome from a holonomy, the statistics from equivariance — and an
electron/muon ratio that is a fine-structure phenomenon with its one
remaining O(1) confined, by a parameter-free law, to the Compton edge of
the natural domain. Two lab programs can now falsify BAM within the
decade; one bounded number separates its mass ladder from an outright
prediction. The arc closes green, its conditions stated, its falsifiers
armed.

## References

- The BAM arc: PR #204 (no-signaling audit), #205 (SN phenomenology),
  #206 (configuration space from the bridge), #207 (swapping as
  surgery), #208 (GHZ as valence), #209 (the measurement sector), #210
  (the strong-field core solve).
- The suppression law: PR #202 (the exact 5D throat-core law);
  #201 (the Dirac-tower ladder); #55–#58 (the EM-capped radius);
  #184 (α as a protected boundary invariant).

## Reproduce

```bash
python -m experiments.closure_ledger.compton_edge_capstone_probe
# Verdict: RELEASE_II_GREEN_THE_ENTANGLED_SECTOR_CLOSED_OPERATIONALLY
#          _THE_MASS_LADDER_ON_THE_ALPHA_ANCHOR_THE_O1_CONFINED_TO
#          _THE_COMPTON_EDGE_OF_THE_NATURAL_DOMAIN
```

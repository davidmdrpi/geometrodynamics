# ℏ-origin research plan

The closure-ledger experiment (`experiments/closure_ledger/`) and the
five-probe quark β-derivation (`docs/quark_beta_status.md`) have
brought the BAM framework to the point where the **action quantum
itself** becomes the next sharp target. This document lays out the
problem structure, identifies the sub-targets, and points at the
concrete numerical tests each sub-target supports.

## Statement of the problem

BAM produces dimensionless ratios (mass ratios m_μ/m_e, m_τ/m_e at
sub-percent agreement; closure-quantum integers N_l = 100, N_q = 466
identified at the locked baseline). The absolute MeV scale is set by
**anchoring m_e** — there is currently no first-principles derivation
of m_e c² in geometric units. Equivalently, there is no
first-principles derivation of ℏ.

The closure-ledger machinery has, however, **localized the natural
geometric object** that an ℏ-derivation should produce. Layer 1 of
the ledger is universal mod 2π (odd-k closure lemma), and the locked
constants — `action_base = 2π`, the Hopf+throat sum = 2π,
`4β_lepton = 100·(2π)` — are all integer multiples of 2π. The
"2π" unit IS the candidate geometric ℏ.

Stating this carefully: the closure-cycle phase is

```
Φ_cycle  =  ∮(L / ℏ) dt   [along the closed worldline]
         =  2π · N         [integer N from quantization]
```

In geometric units (R_MID = c = 1), the integral on the LHS becomes
a pure number, and the RHS is `2π · N`. The conversion to SI gives

```
ℏ_SI  =  (∮L dt)_geom  ·  (m_e · R_MID · c) / (2π · N)
```

So the question "where does ℏ come from" decomposes into:

  (a) Is the closure cycle phase exactly `2π · N` for some integer
      N? (Layer 1 says yes mod 2π. Layer 2 is open.)

  (b) What sets the conversion factor `m_e · R_MID · c`? In CGS,
      `ℏ / (m_e c) = 3.86 × 10⁻¹¹ cm` is the reduced Compton
      wavelength of the electron. Identifying `R_MID = ℏ/(m_e c)`
      makes the conversion factor self-consistent — but currently
      this is an identification, not a derivation.

  (c) Is m_e itself derivable? (THESIS.md "Self-consistent throat
      radius" target. Currently impossible in this framework.)

## Sub-targets, in order of tractability

### (1) Layer 2 closure of the closure ledger (most concrete)

**Status:** Open. Best Layer-2 candidates: C1 at 0.326 rad circular
spread, D1 at 0.577 rad. Neither closes mod 2π.

**Question:** Does the radial bulk channel admit a Layer-2 closure
formula that makes the **full** closure cycle (Layer 1 + Layer 2)
universally `2π · N` for some N(species)?

**Sub-probes already attempted (and what they found):**
- 11 S(k) → {(l, n)} bridge candidates A/B/C/D/maslov: none close.
- Dynamic phases (moving throat, Hopf loop, antipodal back-and-forth,
  worldline crossings): cannot redistribute the residual ~0.3 rad.
- Geometric Hamiltonian + radial matrix elements: insufficient
  dynamic range to reproduce the lepton mass ladder, let alone
  close the ledger.
- Composed Hamiltonian (closure quantum + radial matrix elements):
  τ row is geometrically explained, μ row has factor-100 gap.
- Pinhole-origin chain: μ-row pinhole γ ≈ 22.5 IS the QCD-style
  barrier sum `Σ_{l=0..5} V_max(l)` extended to include the 5D-
  specific l=0 channel.

**Next sub-probe candidate:** the closure cycle treated as a single
worldline integral, NOT as Layer 1 + Layer 2 added separately. The
WKB action `∮p·dq` over the full closed orbit (angular + radial)
might close mod 2π even when the per-channel decomposition does not.

### (2) Aharonov-Bohm form for the action quantum

**Status:** Conceptual; numerical comparison is straightforward.

The Hopf connection `A = (1/2) cos(χ) dφ` has holonomy `π · cos(χ)`
around a fibre at hyper-latitude χ. At χ = 0 (canonical fibre), the
holonomy is π per loop. A spinor going around twice (4π, the
double-cover closure) accumulates `2π` — exactly the action quantum
candidate.

**Sub-probe:** Compute the closure-cycle action for an Aharonov-Bohm
loop around the Hopf fibre and verify it equals `2πℏ` per spinor
double-cover closure. Identify what part of the closure ledger
corresponds to the AB contribution.

### (3) Tangherlini eigenvalue discreteness

**Status:** Conceptual.

The radial Schrödinger problem on Tangherlini has discrete spectrum
{ω_n}. The lowest eigenvalue ω_0 ≈ 1.05 in geometric units sets a
natural frequency scale. Identifying `m_e c² = ℏ ω_0 · (some natural
factor)` could give ℏ as a derivative of m_e and ω_0.

**Sub-probe:** Test whether a natural relationship like `ℏ ω_0 = m_e c²`
or `ℏ (ω_1 - ω_0) = m_e c²` is consistent across species (electron at
k=1, muon at k=3, tau at k=5).

### (4) R_MID self-consistency (deepest, least tractable)

**Status:** Open per `docs/THESIS.md` "Open problems".

R_MID is currently imposed as a geometric constant. The deeper
version determines R_MID as the equilibrium throat radius for the
locked mass spectrum. With R_MID derived, the conversion `ℏ = m_e ·
R_MID · c` becomes a prediction in physical units.

This is unlikely to be solved in the closure-ledger framework alone;
it requires throat dynamics that the present codebase doesn't
implement.

## What success looks like

Three falsifiable predictions a successful ℏ-derivation should yield:

  **P1.** The closure cycle phase Φ_cycle = `2π · N(species)` with
  N a small integer derivable from species-level inputs (k, color,
  Z₂ partition).

  **P2.** A natural identification `R_MID = ℏ/(m_e c)` (the reduced
  Compton wavelength of the electron) emerges from the closure-cycle
  quantization, NOT as a convention.

  **P3.** The Tangherlini lowest eigenvalue ω_0 ≈ 1.05 has a
  geometric meaning in units of m_e c², predicting an absolute mass
  scale rather than relying on the m_e anchor.

**P1** is the most-tractable test. It refines what the closure-
ledger experiment was already doing, by treating the cycle as a
single worldline integral rather than Layer 1 + Layer 2.

## What failure looks like

If no Layer-2 form closes the ledger AND no natural ω_0 ↔ m_e
relationship emerges, the honest verdict is: **BAM cannot derive ℏ
from its current geometric input set**. The framework is dimensional-
ratio-complete (it predicts m_μ/m_e, m_τ/m_e at sub-percent) but
dimensional-scale-incomplete (it cannot predict the absolute MeV
scale or ℏ from geometry alone).

Recording that as a clean negative result is also progress, per the
THESIS.md framing.

## Cross-references

- `docs/odd_k_closure_lemma.md` — the Layer-1 closure mod 2π
  invariance that motivates the action-quantum framing.
- `docs/quark_beta_status.md` — closure of the quark β-derivation
  thread, which produced the `n_part = N_q/2` structural reading.
- `docs/THESIS.md` "Open problems" — original framing of the ℏ-origin
  question, refined here.
- `experiments/closure_ledger/` — all eleven probes building up to
  this research target.

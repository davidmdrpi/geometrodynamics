# Transport / resistance derivation — research plan

Opens the follow-up thread to the ℏ-origin sequence (probes 1–9 under
`experiments/closure_ledger/`). The closure-cycle integer-quantization
and the R_OUTER self-consistency picture are complete; the closure
note (`docs/hbar_origin_note.md` §3, "What this leaves open") flagged
two residual phenomenological dependencies in the locked lepton block:

  - `transport_strength` ≈ 25.1   (off-diagonal coupling magnitude)
  - `resistance_scale`   ≈ 0.2179 (diagonal action coefficient)

The R_OUTER fixed point is sensitive to these at the ~1–7 % level
(probe 8: 7 % R-shift per 1 % transport change, 3 % per 5 % resistance
change). Removing this residual dependency is what "lifting R_OUTER
from phenomenological to fully geometric" means in concrete terms.

## Opening probe — transport / resistance origin probe

`experiments/closure_ledger/transport_resistance_origin_probe.py` is the
first probe in this thread. It mirrors the structure of
`pinhole_origin_probe.py` (which closed γ ≈ 22.5 to Σ V_max). Six
candidate categories are scanned for each target:

  (A) Closure-quantum integers (N·π, N·2π, β/m).
  (B) Tangherlini barrier sums Σ V_max(l) on the canonical grid.
  (C) Eigenfrequency invariants ω(l, n)² and reciprocals.
  (D) WKB tunneling integrals κ_{l1, l2} / dk between shells.
  (E) Cross-shell overlap integrals ⟨u_l | w | u_{l'}⟩ (transport).
  (F) Inverse and log forms (resistance).

A candidate "explains" a parameter if it is within ≤ 1 % of the locked
value AND keeps m_μ, m_τ errors within ≤ 5 % when substituted into the
locked block.

## Findings (opening probe, 2026-05-13)

**Transport.** Best candidate is `8·π = 4·(2π) = 25.1327`, agreeing
with the locked 25.1 to +0.130 %. Three closure-quantum readings
collapse onto the same value (`4·2π`, `8·π`, `16/2·π`). The barrier
sum `Σ V_max[1..5] + π ≈ 25.15` (+0.199 %) lands almost on top.

The reading `transport = 8π` is structurally distinguished: it is the
**4th closure quantum**, parallel to the antipodal k·2π, the Hopf-throat
1·2π, and the τ-uplift 100·2π that organise the Layer-1 ledger
(`hbar_origin_note.md` §2). Transport is naturally a coupling between
4 separate winding sectors (one per pair (1,3), (3,5), (1,5) plus the
self-channel) — the closure-quantum reading aligns with this counting.

**Resistance.** No single candidate hits within 0.5 % of the locked
0.2179. Two readings hit within 1 %:

  - `7π / 100 = 0.2199`  (+0.937 %) — a small closure-quantum fraction.
  - `4·(ω(1,0) − 1) = 0.2189` (+0.477 %) — the 1.054-factor gap, scaled.

The closure-fraction reading and the eigenfrequency reading cannot be
distinguished at this probe's resolution. Both fit the same scale; the
mass-ratio sensitivity tests cannot resolve which is the structural
origin.

**Joint test.** Substituting BOTH transport → 8π AND resistance → 7π/100
into the locked block reproduces the lepton ladder at err μ = 0.59 %,
err τ = 0.65 % — well within the 5 % envelope. Single-parameter
substitutions break the ladder by ~8 % (transport-only) or ~3 %
(resistance-only); the joint reading partially cancels. This is the
non-trivial structural test the closure-quantum reading passes.

## Sub-targets, in order

### (1) Closed-form sharpening of `resistance`

**Status.** Two within-1 % candidates; neither survives mass-sensitivity
at sub-percent precision.

**Question.** Is resistance a closure-quantum fraction (`7π/100` or
similar small form), or is it the radial-Tangherlini eigenfrequency
quantity `4·(ω(1,0) − 1)` evaluated at the locked R_OUTER ≈ 1.262?

**Concrete probes.**
- Re-run with R_OUTER bisected as in probe 8, with each candidate
  resistance reading held fixed. If one of them yields the same
  cross-species fixed-point R* ≈ 1.262 at the original 0.008 %
  tolerance and the other shifts, that selects the geometric origin.
- Extend the candidate set: barrier-integrals with explicit
  ω-dependence (∫ V·sgn(V − ω²) dr*), normalized WKB κ integrals on
  alternative dk conventions, and 5D-specific l = 0 corrections to
  Σ V_max-style operators.

### (2) Transport overlap interpretation

**Status.** Closure-quantum reading 8π works to ~1 % via cancellation,
but no Tangherlini overlap integral hit within 1 % in the opening
probe.

**Question.** Is the transport closure-quantum 4·(2π) the right
structural object, or does it secretly decompose as a sum of
cross-shell overlaps ⟨u_l | V_{l,l'} | u_{l'}⟩ that just happens to
sum near 8π?

**Concrete probes.**
- Identify the OPERATOR whose 1, 3, 5-overlap sum is closest to 8π.
- If multiple operators land near 8π only after re-scaling, that is
  evidence the 8π reading is structural rather than emergent.
- The `experiment_transport_overlap.py` script (already in `scripts/`)
  is for the quark sector; mirror it for the lepton sector and
  cross-reference both readings.

### (3) Close the R_OUTER self-consistency loop on principled inputs

**Status.** The R_OUTER probe (probe 8) bisects with the locked
phenomenological transport / resistance. With principled inputs from
(1) and (2), re-run the bisection.

**Question.** Does the cross-species fixed point R* still land at
1.262 with 0.008 % cross-species tolerance, OR does it shift?

**Concrete probes.**
- Re-run probe 8 with transport = 8π, resistance = 7π/100.
- Re-run probe 8 with transport = 8π, resistance = 4·(ω(1,0) − 1).
- Compare R* values; the one closer to 1.262 selects the principled
  resistance reading.

### (4) Structural unification with γ and 1.054

**Status.** The γ ≈ 22.5 origin is `Σ V_max[1..5] ≈ 22.0` to ~2 %
(pinhole-origin probe); the 1.054 factor is ω(1, 0) at R* ≈ 1.262
(factor-1054 probe, negative closed-form result); the resistance
candidate `4·(ω(1,0) − 1) = 0.219` is structurally the 1.054 factor
in disguise (`4·0.054 = 0.218`).

**Question.** Are γ, the 1.054 factor, and resistance three
projections of one Tangherlini matrix-element family on the same
geometry?

**Concrete probes.**
- Build the explicit linear combination Σ V_max[1..5] + α·(ω(1,0) − 1)
  + β·(something) and ask whether the locked γ, transport, and
  resistance ALL emerge for one consistent (α, β).
- This is the closure-ledger version of the "principled coupling"
  question that opened the ℏ-origin thread for the lepton spectrum.

## Stopping condition

The thread closes when:

  (a) Transport and resistance are both expressible in closed form
      from (Tangherlini grid quantities, closure-quantum integers,
      π, R*) with mass-fit precision matching the locked baseline
      to ~0.5 % or better, OR

  (b) A clean negative result is established: neither parameter
      admits a closed form at the precision the lepton ladder
      requires. In that case the thread documents the
      dimensional-scale-incompleteness as a structural property of
      the framework, not a contingent calibration.

## Cross-references

- `docs/hbar_origin_note.md` — closure-cycle picture and the
  "What this leaves open" agenda this thread inherits.
- `docs/hbar_origin_status.md` — eight-probe summary table.
- `experiments/closure_ledger/transport_resistance_origin_probe.py` —
  opening probe (this thread).
- `experiments/closure_ledger/pinhole_origin_probe.py` — structural
  template (γ → Σ V_max).
- `experiments/closure_ledger/R_outer_self_consistency_probe.py` —
  source of the sensitivity numbers (probe 8).
- `experiments/closure_ledger/factor_1054_search_probe.py` —
  parallel attempt for the 1.054 factor; negative closed-form result.
- `scripts/experiment_resistance_wkb.py`,
  `scripts/experiment_transport_overlap.py` — quark-sector
  precedents for the same operators.

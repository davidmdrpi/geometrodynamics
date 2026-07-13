# The wormhole-network confirmation (PR #216)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. #213 derived the Feynman propagator
> from retarded + advanced completions, but the advanced half sat where
> Wheeler left it: waves propagate forward in time, and the engine's
> `advanced_confirm_amplitude` assigns the confirmation's conjugated
> phase by fiat. This PR supplies the mechanism: an explicit
> wormhole-network traversal, everywhere locally future-directed, whose
> exterior projection **is** the advanced kernel. New module
> `geometrodynamics/transaction/network.py`; the companion probe
> machine-checks every claim (~30 s).

## 0. The mechanism

A wormhole network can carry a retarded wave into the global past
without any locally-backward propagation:

1. **emit** one retarded C-wave (t = 0, ψ = 0);
2. **propagate forward** to the future antipode (t = π — the #166/#213
   refocusing);
3. **transmit** into the throat with the #215 greybody amplitude
   t_AB(ω);
4. **traverse the throat forward** (τ_th > 0 in mouth-clock time);
5. **emerge** through a clock-offset mouth at
   t = π + τ_th + Δ_BA **< 0** — the global past, reached through the
   frozen desynchronization of the mouths (differential aging: one
   mouth parked in a deep gravitational well — the
   Morris–Thorne–Yurtsever mechanism);
6. **propagate forward again** on S³;
7. **intersect the original particle crossing** (t_return = t_emit).

The traversal ledger (probe T2, ω = 3, τ_th = 0.8, Δ_BA = −(2π+0.8)):

| leg | global start | global end | local duration |
|---|---:|---:|---:|
| offer: source → mouth A | 0.000 | 3.142 | **+3.142** |
| throat | 3.142 | **−3.142** | **+0.800** |
| return: mouth B → source | −3.142 | 0.000 | **+3.142** |

Every leg is future-directed in its own clock; the only pastward jump
is the mouth embedding offset — a property of the frozen geometry, not
of the propagation.

## 1. The data structures

`NetworkMouth`: `psi`, `link_id`, `clock_rate`, `clock_offset`,
`orientation`, `transfer_phase`. `NetworkThroat`: `mouth_A`, `mouth_B`,
`tau_th`, `t_AB(ω)` (the #215 complex greybody), with `delta_BA` and
`U_BA(ω)` **derived**, not free:

**The central identity.** A wave entering mouth A at global t exits
mouth B at t′ = t + τ_th + Δ_BA carrying local phase e^{−iωτ_th};
re-expressed against the global form e^{−iωt′}, the transfer factor is
forced by clock continuity to be

```
U_BA(ω) = e^{iωΔ_BA} × (orientation & mouth phases)
```

— **the clock offset *is* the transfer phase** (machine-checked to
10⁻¹² across frequencies, |U| = 1). Nothing about the traversal's phase
is assigned; it is all bookkeeping of one wave read in two clock
systems.

## 2. The projection is advanced

The absorption → return segment, in exterior labels, spans
D′ = τ_th + Δ_BA + d_B **< 0** and carries the kernel

```
K = Λ(ω) · e^{−iωD′} ,        Λ(ω) = t_AB(ω) e^{iωΔ_BA} × decorations
```

— **the retarded phase rule analytically continued to a negative
exterior interval: the advanced kernel**, with greybody weight
|Λ| = √T (machine identity to 10⁻¹²). The response precedes its cause
in global time while every local segment is retarded. This is the
theorem the PR exists for: *a clock-offset wormhole traversal projects
onto the exterior as advanced propagation.*

**Phase closure is #213's coherence condition.** The confirmation
enters coherently (the #213 relative phase φ = 0) exactly when

```
arg t_AB(ω) + ωΔ_BA = 0   (mod 2π)
```

— satisfied on a **discrete frequency comb** with spacing
2π/|Δ_BA + τ_W| (τ_W = the throat's Wigner delay, measured 2.30 → 0.08
from below to above the barrier); measured comb spacing within 4% of
the prediction. The transaction *selects its modes* — transactional
selectivity from geometry. And the weld to the existing engine: the
repo's `retro_phase_match` accepts the network loop amplitude at a comb
frequency (weight 1.000) and rejects it at quarter detuning (< 0.5).

## 3. The live packet

A Gaussian packet (ω₀ = 3, above the barrier) run through the closed
network, component by component with the full complex greybody:

- **emerges in the global past**: envelope peak at t = −3.06 (the
  geometric t_e = −π plus the Wigner delay);
- **intersects the original crossing**: return envelope peaks at
  t = 0.088 vs the Wigner delay 0.085 — the packet lands on its own
  emission event;
- **energy closure**: returned fraction 0.9998 = mean T over the packet
  (flux conservation |r|² + |t|² = 1 at 3×10⁻⁴; unitary transfer
  |U| = 1);
- a **below-barrier** packet (ω₀ = 0.5) returns only T ≈ 0.33 —
  **partial confirmation**: the #215 IR transparency reappears as weak
  transactional confirmation of soft modes.

## 4. The pole structure

The two Compton completion families are now **two globally causal path
classes**: the direct exterior path (the offer ordering I₊) and the
network path (the confirmation ordering I₋), with weights 1 and Λ(ω):

- **at closure** (a comb mode above the barrier): |Λ| = 0.999996,
  arg Λ = 2×10⁻⁶, and K = I₊ + ΛI₋ reproduces the covariant pole to
  3×10⁻⁶ — #213's G_F, now assembled from two future-directed paths;
- **detuning the network** (Δ_BA → Δ_BA + δ) produces relative phase
  φ = ωδ *exactly* (10⁻¹⁵), and the pole-form deviation matches #213's
  T5 deform diagnostic to 10⁻⁵ relative: **the #213 refutation edge now
  has a geometric knob** — the deform test measures mouth-clock
  detuning;
- **below the barrier**, K = I₊ + √T·I₋ splits exactly into a coherent
  pole share (1+√T)/2 and an ordering-asymmetric deficit (1−√T)/2 —
  partial confirmation with the deficit controlled by the greybody.

## 5. Honest scope

- The network **breaks global hyperbolicity by construction** (a
  time-machine geometry). The closure condition is exactly the Novikov
  self-consistency fixed point — transactions are its solutions. No
  chronology-protection dynamics is addressed; the network is frozen
  background, per the program framing.
- The network geometry is **posited, not solved**: the mouth pair from
  the #58/#200 nucleation channel, Δ_BA from differential aging (MTY) —
  cited as mechanism, not derived from the 5D field equations.
- Single-mode scalar, zonal reduction; greybody = the ℓ = 0 channel.
- Elastic case (equal mouth clock rates) analyzed; unequal rates
  red/blueshift the traversal (`emergent_frequency`) — unit-tested
  only.
- `advanced_confirm_amplitude` is retained; `network_confirmation` is
  the mechanism-level replacement, and the T4 weld shows the engine's
  phase closure accepts it. Full engine migration is staged.
- Which tower modes land on the closure comb depends on network
  parameters — selectivity is demonstrated, not matched to a spectrum.

## 6. What would falsify this

- A locally-backward segment anywhere in the ledger — the "globally
  causal" claim would be empty. (Checked: all local durations > 0.)
- A projected kernel that is not Λ·e^{−iωD′}, or |Λ| ≠ √T — the
  advanced identification would fail. (Checked: identities to 10⁻¹².)
- A packet missing the crossing by more than the Wigner delay — time
  closure would fail. (Checked: 0.088 vs 0.085.)
- Network detuning ≠ the #213 deform phase — the two completion
  families would not be these path classes. (Checked: φ = ωδ exact.)

## 7. Deliverables

`geometrodynamics/transaction/network.py` (NetworkMouth, NetworkThroat,
traversal + projection functions; 9 unit tests in
`tests/test_transaction_network.py`);
`experiments/closure_ledger/wormhole_network_confirm_probe.py` (T1–T8,
~30 s) with the full #215 greybody.

**Verdict:**
`THE_ADVANCED_HALF_IS_A_FUTURE_DIRECTED_NETWORK_TRAVERSAL_THE_CLOCK_OFFSET_MOUTH_IS_THE_ANALYTIC_CONTINUATION_AND_PHASE_CLOSURE_IS_THE_COHERENCE_CONDITION`

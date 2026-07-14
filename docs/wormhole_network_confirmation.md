# The wormhole-network confirmation (PR #216)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. #213 derived the Feynman propagator
> from retarded + advanced completions, but the advanced half sat where
> Wheeler left it: waves propagate forward in time, and the engine's
> `advanced_confirm_amplitude` assigns the confirmation's conjugated
> phase by fiat. This PR supplies the mechanism: an explicit
> wormhole-network traversal — through a **two-port throat with
> repeated-loop interior physics** — everywhere locally
> future-directed, whose exterior projection **is** the advanced
> kernel. New module `geometrodynamics/transaction/network.py`; the
> companion probe machine-checks every claim (~90 s).

## 0. The mechanism

A wormhole network can carry a retarded wave into the global past
without any locally-backward propagation:

1. **emit** one retarded C-wave (t = 0, ψ = 0);
2. **propagate forward** to the future antipode (t = π — the #166/#213
   refocusing);
3. **transmit** through mouth A's interface with the #215 greybody
   amplitude t_A(ω);
4. **traverse the throat forward** (τ_th > 0 in mouth-clock time),
   **looping k times** between the two interior barrier faces — each
   loop future-directed, adding 2τ_th of local time;
5. **exit through mouth B's own interface** (t_B(ω)) and **emerge**
   through the clock-offset mouth at
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
`orientation`, `transfer_phase`. `MouthPort`: the full one-sided
S-matrix of a mouth's barrier — `t(ω)`, `r_out(ω)`, `r_in(ω)` (the
face the loops bounce off), from the #215 solver run from **both
sides**: reciprocity t_in = t (2×10⁻⁴), interior flux closure (10⁻⁶),
and the unitarity relation r_in = −r̄·t/t̄ (6×10⁻⁵) all measured.
`NetworkThroat`: `mouth_A`, `mouth_B`, `tau_th`, `port_A`, `port_B`,
with `delta_BA`, `U_BA(ω)` **and the composite transmission t_AB(ω)
derived**, not free:

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

## 1b. The two-port throat: repeated loops, validated directly

The composite transmission sums the interior Fabry–Pérot loops:

```
t_net(ω) = t_A t_B / (1 − r_inA r_inB e^{2iωτ_th})
```

- **transparent-port reduction**: open interfaces give t_net = 1
  exactly (the v1 single-interface model is the r_in → 0 limit);
- **loop-series convergence**: the echo expansion
  t_A t_B (r_inA r_inB e^{2iωτ})^k partial-sums to the closed form
  (10⁻¹⁰);
- **direct validation**: a from-scratch solve of the *glued
  two-barrier potential* V(y) = V_x(x_pk−(y+c)) + V_x(x_pk+(y−c))
  matches the composition to **3×10⁻⁴ in |T|²** across a sweep through
  the resonances — the loop physics is not a model ansatz, it is what
  the differential equation does;
- **resonant confirmation** (new physics): identical mouths transmit
  **perfectly on the interior resonance comb even deep below the
  barrier** — at single-port T = 0.33, the resonant T_net = 0.9996
  versus 0.040 off resonance (the measured off-resonance floor equals
  T²/(1+R)² exactly). The #215 IR transparency is punctured by the
  throat's own interior modes;
- **the Airy identity**: the frequency-averaged T_net equals the
  incoherent echo-train sum T²/(1−R²) (1.5×10⁻⁴) — coherent loops and
  the echo picture are the same physics.

**The echo train** (below-barrier, ω = 0.5): emergences at
t_absorb + (2k+1)τ_th + Δ_BA, geometrically damped by |r_in|² = 0.667,
every echo future-directed in the throat clock:

| k | local duration | global emergence | amplitude |
|---:|---:|---:|---:|
| 0 | 0.80 | −3.142 | 0.333 |
| 1 | 2.40 | −1.542 | 0.222 |
| 2 | 4.00 | +0.058 | 0.148 |
| 3 | 5.60 | +1.658 | 0.099 |

(Time closure is tuned to the primary; the k-th echo returns 2kτ_th
after the crossing — a network can instead close on any single echo.)

## 2. The projection is advanced

The absorption → return segment, in exterior labels, spans
D′ = τ_th + Δ_BA + d_B **< 0** and carries the kernel

```
K = Λ(ω) · e^{−iωD′} ,        Λ(ω) = t_net(ω) e^{iωΔ_BA} × decorations
```

— **the retarded phase rule analytically continued to a negative
exterior interval: the advanced kernel**, with the composite two-port
weight |Λ| = |t_net| (machine identity to 10⁻¹²). The response precedes its cause
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

- **emerges in the global past**: envelope peak at t = −2.96 (the
  geometric t_e = −π plus the composite Wigner delay);
- **intersects the original crossing**: return envelope peaks at
  t = 0.180 vs the composite Wigner delay 0.170 (both interfaces
  contribute) — the packet lands on its own emission event;
- **energy closure**: returned fraction 1.0000 = mean T over the
  packet (flux conservation at both interfaces; unitary transfer
  |U| = 1);
- **resonant confirmation, live**: a narrow below-barrier packet
  (ω₀ = 0.5) centered **on** an interior resonance confirms at 0.88 —
  stored in the throat cavity and released over the storage time (the
  composite Wigner delay ≈ 58; return peak at +48) — while **off**
  resonance it confirms at only 0.040 (**22× contrast**). Soft modes
  confirm through the throat's interior comb or barely at all.

## 4. The pole structure

The two Compton completion families are now **two globally causal path
classes**: the direct exterior path (the offer ordering I₊) and the
network path (the confirmation ordering I₋), with weights 1 and Λ(ω):

- **at closure** (a comb mode above the barrier): |Λ| = 0.9997 and
  K = I₊ + ΛI₋ reproduces the covariant pole at the |1−Λ| level —
  #213's G_F, now assembled from two future-directed paths;
- **detuning the network** (Δ_BA → Δ_BA + δ) produces relative phase
  φ = ωδ *exactly* (10⁻¹⁵), and the pole-form deviation matches #213's
  T5 deform diagnostic to 10⁻⁵ relative: **the #213 refutation edge now
  has a geometric knob** — the deform test measures mouth-clock
  detuning;
- **below the barrier**, K = I₊ + |t_net|·I₋ splits exactly into a
  coherent pole share (1+|t_net|)/2 and an ordering-asymmetric deficit
  (1−|t_net|)/2 — partial confirmation with the deficit controlled by
  the composite greybody, i.e. by where the mode sits relative to the
  interior comb.

## 5. Honest scope

- The network **breaks global hyperbolicity by construction** (a
  time-machine geometry). The closure condition is exactly the Novikov
  self-consistency fixed point — transactions are its solutions. No
  chronology-protection dynamics is addressed; the network is frozen
  background, per the program framing.
- The network geometry is **posited, not solved**: the mouth pair from
  the #58/#200 nucleation channel, Δ_BA from differential aging (MTY) —
  cited as mechanism, not derived from the 5D field equations.
- The interior transit τ_th and the mapping of the two mouths'
  near-horizon regions onto a single 1D channel are model-level (the
  glued-potential direct solve validates the composition on exactly
  that model); the true 5D interior is #168/#200 territory.
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

`geometrodynamics/transaction/network.py` (NetworkMouth, MouthPort,
NetworkThroat with the derived two-port composition, traversal +
projection + echo-train functions; 14 unit tests in
`tests/test_transaction_network.py`);
`experiments/closure_ledger/wormhole_network_confirm_probe.py` (T1–T8,
~90 s) with the full two-sided #215 greybody and the direct
glued-barrier validation.

**Verdict:**
`THE_ADVANCED_HALF_IS_A_FUTURE_DIRECTED_NETWORK_TRAVERSAL_THE_CLOCK_OFFSET_MOUTH_IS_THE_ANALYTIC_CONTINUATION_AND_PHASE_CLOSURE_IS_THE_COHERENCE_CONDITION_NOW_THROUGH_A_TWO_PORT_THROAT`

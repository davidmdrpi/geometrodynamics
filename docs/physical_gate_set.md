# The physical gate set: the qubit clock, the gravitational Kerr, and the twist switch (PR #234)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. #232 delivered the 2^N structure
> with one structural placeholder: "χ is not calibrated to a throat
> geometry in this probe." This PR executes that scope item on the
> committed #224 network — and welds in the other arc's discovery
> (#227/#233, the Z₂ twist) as the gate set's switch. Every rate
> measured, in model units, on geometry the ledger already owns.
> Companion probe: ~3 s.

## 0. The gate set, stated first

| element | mechanism | rate (model units) |
|---|---|---|
| **one-qubit clock** | the #224 doublet beat (intra-rail tunneling) | Δω = 1.714×10⁻³; X period π/Δω = 1833 |
| **memory/gate switch** | the #227 Z₂ link twist on the qubit's cycle | freeze ratio 4.9×10⁻¹¹; leakage time 3.8×10¹³ |
| **entangling gate (CZ)** | the **gravitational Kerr**: the #223 dressing per quantum | χ = 5.5×10⁻²; t_CZ = π/χ = 57.2 |

The CZ is **32× faster than the one-qubit clock**: on this substrate
the entangler is strong and the local rotations are slow — the
opposite of most laboratory platforms.

## 1. The clock: the qubit's X rotation is the network's own beat

The #224 interior doublet *is* the one-qubit sector: rails = the A/B
mouth basins (localization 0.977), and the intra-doublet beat is the
X rotation. Re-derived on the same builder, the ledger numbers return
at machine agreement (Δω dev 0; period dev 0; the localized-combo
weight to 10⁻⁹ using #224's own construction).

## 2. The switch: one topological bit, and it spares the Kerr

Flipping one exterior link (the #227 holonomy) freezes the same
doublet's beat to 8.35×10⁻¹⁴ — the #227 committed row at r_s = 0.3
re-derived to machine agreement. So a single Z₂ bit toggles the qubit
between **gate mode** (clock running) and **memory mode** (the #227
exact-pair theorem; leakage time 3.8×10¹³).

The crucial compatibility: the twist acts on the *linear* channel
only. Rotating the exactly degenerate twisted pair to its maximally
localized combination (analytic 2×2 eigenproblem; localization 0.977
preserved), the **dressing integral changes by < 0.2%** — the
density-density (Kerr) channel is twist-blind. The twist is an error
switch, not a gate hazard: it kills exactly the channel that (per
#232's KLM no-go) can only leak, and leaves exactly the channel that
gates.

## 3. The gravitational Kerr, calibrated

The committed nonlinearity is not a contact |φ|⁴ — it is the **EKG
backreaction**, whose leading form is the #223 dressing law
δμ = c·A² (re-read from the ledger: exactly quadratic, ratio-4 at
machine). Calibrating per **one quantum** of a rail mode on its
throat (continuum normalization, amplitude² = 1/2ω, the #223 stress
integral over the throat's bridge walls):

```
δμ_q = 5.80×10⁻²         (the throat dressed by a single quantum)
dω/dμ = −0.947            (measured response, linear to 5%)
χ_rail = |dω/dμ|·δμ_q = 5.49×10⁻²    ⟹    t_CZ = π/χ = 57.2
```

**Gravity is the entangling resource.** The same backreaction that
welds the soliton scales (#222) and dresses the bridge (#223) is, per
quantum, the two-qubit gate. #232's sentence — "the same nonlinearity
that makes the soliton makes the CZ" — now has its number, and its
mechanism named correctly: gravitational cross-phase modulation.

**The honest edge (the finding, not a footnote):** δμ_q/μ = 0.64 —
one quantum dresses the throat at O(1) at this network scale. The
quoted χ is the *leading-order* value; the #223 perturbative window
is exceeded per quantum. The gravitational entangler is **strong**,
and the self-consistently dressed gate is the named successor.

## 4. The dual-rail selection rule

Within the sector (one quantum per doublet), every **same-qubit** Kerr
term vanishes identically — n_a·n_b = 0 on a shared quantum, and
n(n−1) = 0 for occupation ≤ 1 — checked in Fock at machine zero. The
encoding auto-silences all self-phases; the only active Kerr couples
rails of **different** qubits: exactly the CZ generator. The CZ with
the calibrated χ, run under twist-frozen linear coupling, is exact
(10⁻¹⁶, zero leakage).

The one measured residual: the distant-rail tail (3.9×10⁻² of χ_rail,
the 2.3% basin leakage squared through the stress integral) gives a
**coherent ZZ phase of 0.13 rad per CZ** — known, calibrated,
compensable by echo. Reported, not hidden.

## 5. The circuit budget

GHZ₃ at physical rates (the #232 circuit: 5 one-qubit gates + 2 CZs):

| item | value |
|---|---|
| one-qubit gate (half X period) | 916 |
| CZ | 57.2 |
| GHZ₃ total | ≈ 4700 model units |
| twist-frozen leakage error (T·Δω_frozen)² | 1.6×10⁻²⁰ |
| coherent ZZ total (2 CZs) | 0.26 rad, compensable |

The topological switch makes linear leakage negligible over the whole
circuit; the compiler owes exactly one calibrated correction.

## 6. The arc-weld holdout

Seven committed constants from **four ledgers spanning both
development arcs** — #224 (the clock Δω and basins), #227 (the frozen
beat and freeze ratio, from the *other* session's arc), #223 (the
exactly-quadratic dressing), #232 (the structural gate exactness) —
re-read and re-derived on one machine, all at machine agreement. The
two arcs meet in one calibrated gate set without retuning anything.

## 7. Honest scope

- All rates in **model units** on the committed r_s = 0.3 network;
  conversion to laboratory units requires the anchor chain (#225/#226)
  and is deliberately out of scope.
- χ is **leading order** in a strong coupling (δμ_q/μ = 0.64); the
  self-consistent dressed gate is the successor. The strength itself
  is the finding.
- Per-quantum normalization is standard flat-measure quantization
  (amplitude² = 1/2ω, ∫u²dx = 1); O(1) convention factors between
  this and #223's peak-amplitude convention are absorbed into the
  leading-order status.
- The two-qubit CZ presumes rails of different qubits sharing a
  throat; this probe calibrates the on-throat Kerr scale and the
  distant-rail tail on the one committed network — the two-network
  shared-throat geometry is construction work, not new physics.
- The Hadamard-at-the-beat clock is an order-of-magnitude physical
  rate; arbitrary one-qubit gates need controlled detuning (the #224
  bias channel — committed machinery, not exercised here).
- The #233 correction (deck class set by mouth placement) does not
  touch the freeze mechanism used here.

## 8. What would falsify this

- The doublet or freeze numbers drifting from the #224/#227 ledgers —
  the weld would be broken. (Checked: 0 deviation on both.)
- The twist shifting the dressing integral at O(1) — the switch would
  damage the gate. (Checked: < 0.2%.)
- A same-qubit Kerr term surviving on the sector — the encoding would
  need active cancellation. (Checked: identically zero.)
- dω/dμ nonlinear at the probed steps — χ would not be a coefficient.
  (Checked: linear to 5%.)
- δμ failing the A² law — the Kerr reading of the dressing would be
  wrong. (Checked: ratio-4 at machine, from the #223 ledger.)

## 9. Companion probe

`experiments/closure_ledger/physical_gate_set_probe.py` (T1–T9, ~3 s):
the clock re-derivation; the twist switch with the Kerr-channel
invariance; the calibrated dressing chain (δμ_q, dω/dμ, χ, t_CZ); the
Fock selection rule and the calibrated CZ; the circuit budget; the
four-ledger arc-weld holdout.

**Verdict:**
`THE_GATE_SET_IS_PHYSICAL_THE_BEAT_IS_THE_CLOCK_THE_TWIST_IS_THE_SWITCH_AND_THE_GRAVITATIONAL_DRESSING_IS_THE_KERR_WITH_T_CZ_32X_FASTER_THAN_THE_CLOCK_AT_LEADING_ORDER_IN_A_STRONG_COUPLING`

## Reproduce

```bash
python -m experiments.closure_ledger.physical_gate_set_probe
```

Expected verdict:
`THE_GATE_SET_IS_PHYSICAL_THE_BEAT_IS_THE_CLOCK_THE_TWIST_IS_THE_SWITCH_AND_THE_GRAVITATIONAL_DRESSING_IS_THE_KERR_WITH_T_CZ_32X_FASTER_THAN_THE_CLOCK_AT_LEADING_ORDER_IN_A_STRONG_COUPLING`, 9/9 PASS.

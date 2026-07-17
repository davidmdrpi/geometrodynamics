# Mouth-exchange dynamics: P_other-mouth(t) and its laws (PR #224)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. The #223 successor: the five
> requested exchange observables — P_other-mouth(t), the maximum
> transferred probability, the exchange period π/Δω, the survival
> probability under realistic asymmetry, and the complex decay width
> when open channels exist — computed on the genuine two-level system
> of the network. The companion probe machine-checks every claim
> (~2 min).

## 0. The system — and an honest topological finding first

On a **single** two-mouth bridge ring the exterior is *one* connected
cavity (the two brane arcs join at the antipode), so no mouth-localized
doublet exists on it: #223's even/odd splitting is the neck-BC
sensitivity — the correct *transit-coupling measure* — but not a
physical two-state beat. The genuine mouth-to-mouth two-level system
lives on the **two-throat network**: two #221 interior-channel throats
on the shared S³ exterior. The dressed interior eigenhistories of
throat A and throat B are the localized basins; the grid is kept
segment-commensurate so the discrete ring carries the *exact* A↔B swap
symmetry (without it, the tiny coupling loses to grid detuning — a real
lesson of the construction, machine-visible).

The working doublet (r_s = 0.3, interior mode n = 3 at ω ≈ 2.19):
interior fraction 0.98, eigenmodes exactly half-half, basin combos
(e ± o)/√2 localized to **0.977** in one throat.

## 1. P_other-mouth(t)

For the doublet the regional norm obeys the **exact** identity

```
P_B(t) = P_B(0) + [P_max − P_B(0)]·sin²(Δω·t/2)
```

(cosine algebra of the two-mode cross term — machine zero), and the
**direct leapfrog evolution** of the wave equation from the throat-A
state reproduces it: throat A depletes to 10⁻⁴, the transfer maximum
0.977 arrives at t/T = 0.998, and the half-period point sits at
P_max/2 to 2%.

## 2. The maximum transferred probability and the exchange period

```
P_max  =  4a²b²  =  the localization weight of the prepared state
       =  0.977      (limited only by the 2.3% exterior dressing tail)
T_exchange  =  π/Δω  =  1833      (evolution: 1829 — 0.2%)
```

In the symmetric network the transfer is complete up to the dressing
tail; nothing else caps it.

## 3. The exchange-coupling law — frozen at the anchor

The doublet splitting falls with the neck radius with exponent
3.3 → 3.7 across r_s = 0.25–0.5, rising toward the #223 amplitude law
(ω·r_s)⁴ as the neck deepens. Extrapolated along the law to the
primordial anchor (r_s·ω = α): the exchange period scales as α⁻⁴ —
**mouth-to-mouth transfer is dynamically frozen for the physical
electron** (period ~10¹⁰⁺ in throat units even before asymmetry).

## 4. Survival under realistic asymmetry

Realistic networks have no identical throats. A clock-rate bias ε on
throat B (the MTY differential-aging proxy; continuous, unlike
grid-quantized length asymmetry) detunes the basins by
δ = ε·f_B/(2ω), and the textbook two-level laws hold to 1%:

| ε | Δω(ε) | P_max | Lorentzian pred. | survival floor |
|---|---|---|---|---|
| 0 | 1.71×10⁻³ | 1.000 | 1.000 | 0 |
| 0.004 | 1.94×10⁻³ | 0.784 | 0.785 | 0.216 |
| 0.016 | 3.98×10⁻³ | 0.186 | 0.186 | 0.814 |
| 0.064 | 1.43×10⁻² | 0.014 | 0.014 | **0.986** |

- the **Rabi identity** P_max·Δω² = Δω₀² holds to 1% across the scan;
- the **Lorentzian** P_max = Δω₀²/(Δω₀² + δ²) matches pointwise to 1%;
- a **1.5% frequency asymmetry already pins survival above 98.6%**:
  *localization by asymmetry* — the #223 transit protection completed
  dynamically. The exact-degeneracy loophole (complete transfer in a
  perfectly symmetric network) closes: no two physical throats are
  identical, so the dressed eigenhistory stays home.

## 5. The complex decay width

- **Compact network: Γ = 0 exactly.** The closed-ring operator is real
  symmetric — every frequency is exactly real. No open channels exist
  on the compact network, and the eigenhistory persists (#218).
- **Open channel** (a semi-infinite lead beyond throat B — the
  wider-network continuation): the doublet becomes a QNM pair with

```
Γ = 4.0×10⁻⁴ / 4.2×10⁻⁴     (outgoing-BC complex Newton;
                              lead-length-independent to < 1%)
```

  in the **strong-coupling regime** J > Γ/4 (J/(Γ/4) ≈ 8): the
  exchange beat survives the opening under the decay envelope
  e^{−Γt}. The width is shared by hybridization:
  Γ_pair ≈ Γ_direct/2, with Γ_direct the single-throat-on-lead width
  (machine-checked ratio ≈ 2).

## 6. Honest scope

- The topological finding is stated as such: the single-bridge ring
  carries no mouth doublet; the two-throat network does. #223's
  splitting remains the correct transit-coupling measure.
- The interior-channel (#215/#221) reading supplies the bound basins;
  the pure ultrastatic bridge (#223) has no interior state to
  exchange.
- "Probability" is the classical normalized field-norm fraction, as
  everywhere in the program; ℏ enters nowhere.
- The asymmetry is a clock-rate bias (MTY-aging proxy); length
  asymmetry is first-order equivalent but grid-quantized.
- The open channel is a model lead; the compact statement Γ = 0 is
  exact and structural.
- Classical, frozen background; the working point is chosen so the
  exchange period fits direct evolution; the anchor statement is the
  measured-law extrapolation.

## 7. What would falsify this

- A beat deviating from sin² beyond the dressing tail — the two-level
  reduction would fail. (Checked: exact identity + evolution to
  1%/3%.)
- A period ≠ π/Δω. (Checked: 0.2%.)
- The Rabi identity or the Lorentzian failing under bias — the basins
  would not be a two-level system. (Checked: 1% across the scan.)
- A nonzero width on the compact network — a hidden open channel.
  (Structural: real symmetric operator.)
- A lead-dependent QNM width — the outgoing condition mis-imposed.
  (Checked: < 1% between lead lengths.)

## 8. Companion probe

`experiments/closure_ledger/mouth_exchange_dynamics_probe.py` (T1–T9,
~2 min): the network and its doublet ladder; the exact beat and the
direct evolution; the two headline numbers; the coupling law; the
asymmetric survival laws; the closed/open width.

**Verdict:**
`THE_MOUTH_EXCHANGE_IS_A_TEXTBOOK_TWO_LEVEL_BEAT_WITH_PERIOD_PI_OVER_DW_ASYMMETRY_LOCALIZES_BY_THE_RABI_LAW_AND_THE_WIDTH_IS_ZERO_UNTIL_A_CHANNEL_OPENS`

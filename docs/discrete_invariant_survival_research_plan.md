# Discrete invariant survival on the ψ–Φ–q throat-soliton (PR #181)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The continuous geometry carries the discrete charge

The arc has built two things: a **continuous** object — the self-consistent
two-way ψ–Φ–q throat-soliton (#179/#180) — and a **discrete** one — the
winding charge `k` of the throat-order field (the #178 vortex charge, the
#174 odd-k ladder). This probe shows the continuous soliton **carries** the
discrete invariant: the winding charge `Q = (1/2π)∮∇φ·dl` rides the continuous
ψ–Φ–q evolution untouched, because `Q` is a **homotopy invariant** of the map
`q: (loop) → ℂ∖{0}` and the self-consistent soliton keeps `|q| > 0` on the
loop.

The survival criterion is exact and sharp: **`Q` can change only where
`|q| = 0`**. The exceptional event — `|q| → 0` and `Q` jumps — is the phase
slip of **PR #182** (next).

## The setup (on the actual #180 soliton)

Take the self-consistent soliton (the #180 solver, run here, at `M = 3.5`),
and an equatorial loop of radius `R` in its ordered core, where
`ρ = |ψ|² > ρ_c` (so the order field is on, `|q| > 0`). Dress it with winding
`k`: `q = |q| e^{ikφ}`.

A winding-`k` vortex is **sustained** on the loop when the well depth beats
the centrifugal cost:

```
A² = (gρ − a₀) − (κ/R²) k² > 0 .
```

With the #180 constants and the best loop (`R = 0.75`, `ρ = 0.36`):

| k | winding Q | A² = (gρ−a₀) − (κ/R²)k² | sustained? |
|---:|---:|---:|---|
| 1 | 1.0 | +0.153 | yes |
| 3 | 3.0 | +0.082 | yes |
| 5 | 5.0 | −0.061 | **no** (centrifugal beats the well) |

## Survival under continuous evolution

### Wave evolution (norm-conserving / unitary) — all k

Starting from a smoothly perturbed winding-`k` state and evolving by a
split-step wave equation `i∂_t q = −(κ/R²)∂²_φ q + (a₀ − gρ)q + λ|q|²q`, the
charge is conserved to **machine precision** for every `k ∈ {1, 3, 5}`:

| k | max winding deviation ΔW | min \|q\| |
|---:|---:|---:|
| 1 | 3×10⁻¹⁶ | 0.27 |
| 3 | 4×10⁻¹⁶ | 0.33 |
| 5 | 9×10⁻¹⁶ | 0.23 |

The norm-conserving dynamics keeps `|q|` off zero, so the discrete charge
rides the continuous evolution untouched — even the unsustained `k = 5`
survives unitary evolution (the wave does not relax it to a zero).

### Dissipative relaxation (the order field's own dynamics) — sustained k

Under the gradient flow `∂_t q = (κ/R²)∂²_φ q − (a₀ − gρ)q − λ|q|²q` (the same
relaxational dynamics that built the #180 soliton), the windings the soliton
**sustains** (`k = 1, 3`) survive: a perturbed winding-`k` state relaxes back
to the clean vortex, `Q` conserved to ~10⁻¹⁵, `min|q| > 0`, final `Q = k`.

## The criterion — Q changes only through |q| = 0 (→ #182)

The unsustained `k = 5` (`A² < 0`) has **no** ordered vortex on the loop, so
under dissipation the amplitude is driven to `|q| → 0` (`min|q| ≈ 10⁻⁴`); at
that zero the phase unwinds and the charge **slips**, `Q: 5 → 2`. Contrast the
sustained `k = 3`, which keeps `|q| > 0` and holds `Q = 3`.

So **survival ⟺ `|q| > 0`**, exactly: the soliton maintains `|q| > 0` for the
windings it sustains, so those survive; an unsustained winding reaches a zero
and slips. That slip — *exactly how the invariant changes when `q` hits zero*
— is the topology-change event of **PR #182**.

## Rigidity under homotopy

Under 40 random smooth homotopies per sector (independent amplitude and phase
perturbations of the winding-`k` vortex, all keeping `|q| > 0`), the charge is
**unchanged in every case** (40/40 for `k ∈ {1, 3, 5}`). No continuous
deformation that avoids `|q| = 0` can touch the winding — the dynamical
realization of the #173/#174 rigidity (the discrete charge lives outside the
rank-controlled continuous moduli), now on the self-consistent soliton.

## Honest scope

- **Homotopy-invariance is exact** (topological): the winding survives any
  continuous ψ–Φ–q evolution while `|q| > 0`, and the #180 soliton maintains
  `|q| > 0` for the windings it sustains.
- **Reduced geometry**: the amplitude envelope is the #180 **radial** soliton;
  the winding is azimuthal on an equatorial loop. The full 2D/3D
  self-consistent **vortex-line** soliton — `|q| = 0` on the axis, the winding
  back-reacting on `ψ` and `Φ` — is a follow-up.
- **Which rungs survive** is set by the soliton's capacity (`κ`, `R`, well
  depth): here `k = 1, 3` survive the dissipative dynamics, `k = 5` exceeds the
  loop and phase-slips (#182); a wave (norm-conserving) evolution preserves all
  three.
- The realized **physical** ladder is odd-`k` {1, 3, 5} by the #174
  orientability grading (even-`k` is topologically conserved too but excluded
  by orientability); its survival under a **deformed bulk geometry** is
  **PR #183**.
- Weak-field, semi-dynamical, effective constants.

## Reproduce

```bash
python -m experiments.closure_ledger.discrete_invariant_survival_probe
# Verdict: DISCRETE_WINDING_INVARIANT_SURVIVES_CONTINUOUS_EVOLUTION_ON_THE_SOLITON_WHILE_Q_NONZERO
```

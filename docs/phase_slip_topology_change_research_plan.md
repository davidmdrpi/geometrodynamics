# The phase-slip / topology-change event (PR #182)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The one event that changes the invariant

PR #181 showed the discrete winding charge `Q = (1/2π)∮∇φ` survives the
continuous ψ–Φ–q evolution while `|q| > 0`, and that it can change **only**
where `|q| = 0`. This probe dissects that exceptional event — the **phase
slip** — and shows *exactly* how the invariant changes when `q` hits zero.

## The obstruction — changing Q forces a zero

`Q` is a homotopy invariant of `q: S¹ → ℂ∖{0}`; it is locally constant on
nowhere-zero fields. To change `Q` the field **must** leave `ℂ∖{0}` — pass
through a configuration with a zero. This is a *topological obstruction*, not
a dynamical accident.

The straight homotopy `(1−s)·[winding 1] + s·[winding 0]` between the two
sectors is **forced** through an exact zero:

| quantity | value |
|---|---|
| global `min|q|` over the path | `2.5×10⁻¹⁷` (machine zero) |
| at parameter | `s* = 0.5` |
| located at | `φ* = π` |
| winding | `Q: 1 → 0` |

There is **no** nowhere-zero path between the winding sectors. (This is the
dynamical/topological content of the #175 gate: the discrete sector is
reachable only through an amplitude zero.)

## The quantum — ΔQ = ±1 (a 2π phase kink)

Just below the slip (`s = 0.48`) the field has winding `Q = 1`; just above
(`s = 0.52`), `Q = 0`: `ΔQ = −1`. Equivalently the integrated winding density
`∮∇φ` changes by exactly **`−2π`** — one full turn of phase removed as the
field passes through the zero (a 2π phase kink through the zero point,
concentrating at `φ* = π` as `s → 0.5`).

A generic simple zero carries **unit** topological charge in
(space × parameter), so the elementary phase slip is `ΔQ = ±1`. The winding is
an integer before and after; the only non-integer instant is exactly at the
zero.

## The dynamics — a single slip and a quantized staircase

In a genuine ψ–Φ–q evolution on the soliton (the order field's own
dissipative gradient flow), `Q(t)` is piecewise-constant and changes by `±1`
**exactly** at the instants `min|q|(t) → 0`.

- **Single event.** The unsustained winding `k = 5` holds flat at 5, then the
  first slip steps `ΔQ = −1` to 4 exactly when `min|q| → 0`
  (`min|q| = 5×10⁻⁴` at the step).
- **Staircase.** A strongly over-wound state `k = 8` cascades down the
  quantized staircase `[8, 7, 5, 4, 3, 2]`, **every** step coinciding with a
  `min|q| → 0` event — the over-wound throat sheds winding one quantum at a
  time through amplitude-zero nodes. (A recorded `−2` step is two elementary
  slips unresolved in sampling time, not a single `ΔQ = 2` event; the
  elementary quantum is `±1`, established above.)

## Localization — Q is an exact integer except at the zeros

Along the cascade, at every recorded time with `min|q| > 0.1` (away from a
slip), the unrounded winding `(1/2π)∮∇φ` equals an integer to **`10⁻¹⁵`** —
`Q` is exactly quantized while `|q| > 0`. The only times the winding is
ambiguous / changes are the slip instants, where `min|q| → 0`. The
non-integer-ness of the invariant is confined to the **measure-zero** set of
amplitude zeros — between them `Q` is a rigid integer (#181), at them it jumps
by `±1` (#182).

## Physical meaning

The phase slip is the throat changing its winding / generation sector
(`k → k∓1`) **through** the amplitude-zero node — the #175 antipodal node, the
#178 defect core. The #175 "gate" (the discrete sector is reachable only
through an amplitude zero) is here the sector-**changing** event itself,
resolved: each transition removes one unit of winding at a node. Together with
#181 (the sector survives between events), the throat's winding is a conserved
topological charge that transitions **only** at nodes — and the realized
ladder is odd-`k` by the #174 orientability grading (its survival under a
**deformed bulk geometry** is **PR #183**).

## Honest scope

- The obstruction and the `±1` quantum are **exact** (topological).
- The dynamical staircase is on the reduced **vortex-on-soliton** loop
  (amplitude from the #180 radial soliton, winding azimuthal); the elementary
  slip is `ΔQ = ±1`, and a recorded `−2` step is two unresolved-in-time
  elementary slips, not a `ΔQ = 2` event.
- The full 2D/3D **vortex-line reconnection** (the zero as a moving vortex line
  threading the loop) is a follow-up.
- Weak-field, effective constants.

## Reproduce

```bash
python -m experiments.closure_ledger.phase_slip_topology_change_probe
# Verdict: TOPOLOGY_CHANGE_OCCURS_ONLY_AT_AMPLITUDE_ZEROS_EACH_ELEMENTARY_SLIP_CHANGES_Q_BY_PLUS_MINUS_ONE
```

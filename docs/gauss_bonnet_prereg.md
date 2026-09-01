# Pre-registration: does Gauss–Bonnet reopen the throat?

**Status: frozen before the module exists.** Committed alone, as with
`docs/finite_mouth_prereg.md` and `docs/source_audit_prereg.md`.

---

## Why this branch and not the other four

`docs/source_audit_prereg.md` closed with five open branches, each a *premise*
of the flare-out theorem rather than a consequence. Gauss–Bonnet is the one
worth computing first:

- in `D = 5` the Gauss–Bonnet invariant is **dynamical**, not topological — it
  is the natural higher-curvature term exactly here;
- it is the only branch that keeps both a **classical** geometry and a
  traversable throat, so it needs no quantum stress and no ghost;
- it is a calculation, not a philosophical choice.

The field equations are `G_ab + α_GB H_ab = 8πG₅T_ab` with the Lanczos tensor

```
H_ab = 2[R R_ab − 2R_ac R^c_b − 2R^{cd}R_acbd + R_a^{cde}R_bcde] − ½ g_ab L_GB
```

Under null contraction the `g_ab L_GB` term drops, and `G_kk = R_kk`, so

```
R_kk + α_GB H_kk = 8πG₅ T_kk
```

**Raychaudhuri is untouched.** Flare-out is pure geometry, so `R_kk = −3f″/f`
at a neck regardless of the gravitational action. What changes is only *which
tensor* must supply it. The hope for this branch is that `α_GB H_kk` supplies
the negative part geometrically, leaving `T_kk ≥ 0`.

---

## Verified before freezing

The Lanczos formula was validated first, because everything below depends on
it. In `D = 4` the Gauss–Bonnet term is topological, so `H_ab` must vanish
identically:

| check | result |
|--|--|
| `H_kk`, 4D Schwarzschild | `0` |
| `H_kk`, 4D **general** `A(r)` (non-vacuum, so not trivially zero) | `0` |

Then, for `−N(s)²dt² + ds² + f(s)²dΩ₃²` in `D = 5` with `k = (1/N, 1, 0,0,0)`
verified null:

```
R_kk = −3 (N f″ − N′f′) / (N f)
H_kk = 12 (f′² − 1)(N f″ − N′f′) / (N f³)
```

---

## The six frozen predictions

**G1 — the lapse drops out at the neck, again.** `N′` multiplies `f′`, and
`f′(0) = 0` is what makes `s = 0` a neck. So for **any** `N, N′, N″`:

```
R_kk = −3f″/f₀ ,      H_kk = −12f″/f₀³
```

**G2 — Gauss–Bonnet has the same sign as Einstein at every neck.**

```
H_kk / R_kk = 4(1 − f′²)/f²   →   4/f₀²  > 0   at a neck
```

It **reinforces** the Einstein term rather than cancelling it. This is the
prediction that decides the branch, and it is the opposite of what the branch
was invoked to do.

**G3 — the ratio is the Misner–Sharp parameter.** Since `μ = f²(1−f′²)` is the
same quantity P2 showed is continuous across the seam,

```
H_kk = (4μ / f⁴) · R_kk
```

so `α_GB H_kk` reinforces `R_kk` wherever `μ > 0`. For this geometry `μ = b²`
everywhere, so the reinforcement holds along the whole throat, not only at the
neck.

**G4 — the matter NEC requires a negative coupling.**

```
8πG₅ T_kk|₀ = −3f″(4α_GB + f₀²)/f₀³
```

With flare-out `f″ > 0`, `T_kk ≥ 0` demands

```
α_GB ≤ − f₀²/4
```

**G5 — the numerical threshold for the finite mouth.** With `f₀ = b = R sin²a`:

```
α_GB ≤ −R² sin⁴a / 4
```

At `R = 1, a = 0.3` that is `α_GB ≤ −0.001906728`.

**G6 — at threshold the expansion has broken down.** The GB contribution
relative to Einstein is `α_GB H_kk / R_kk = 4α_GB/f₀²`, which at the threshold
is exactly `−1`: the "correction" equals the leading term. Truncating Lovelock
at Gauss–Bonnet is then unjustified — the whole tower contributes.

---

## The frozen verdict

**Gauss–Bonnet does not reopen the throat**, and it fails three ways at once:

1. **Wrong sign.** Heterotic string theory gives `α_GB = α′/8 > 0`; the
   string-motivated sign makes the NEC violation *worse*, by the factor
   `(1 + 4α_GB/f₀²) > 1`.
2. **Wrong magnitude even if the sign is granted.** `|α_GB| ≥ f₀²/4` ties the
   Gauss–Bonnet length to the throat radius.
3. **Outside its own regime.** At that magnitude the correction is 100% of the
   leading term, so the derivative expansion has failed.

---

## Falsifiers

1. Any `H_kk` computation disagreeing with the closed forms above falsifies G1
   or G2 — an independent numerical differentiation of the curvature invariants
   must reproduce them.
2. A non-vanishing `H_ab` in `D = 4` falsifies the Lanczos implementation, and
   with it everything else.
3. `H_kk/R_kk < 0` at any neck of any admissible profile falsifies G2, and
   would **reopen** the branch — the single most consequential possible
   outcome.
4. A positive `α_GB` satisfying `8πG₅T_kk ≥ 0` at the neck falsifies G4.
5. If `μ < 0` were attainable on an admissible profile, G3 would let the sign
   flip; the module must check `μ` rather than assume it.

---

## Scope, and what remains after this

This treats **Einstein–Gauss–Bonnet with a constant coupling**. It does not
address: a *dilatonic* EGB coupling `α(φ)L_GB`, where the scalar's own stress
enters and known 5D solutions exist; the full Lovelock tower; or `f(R)`-type
modifications. Those are separate premises and are not refuted here.

If G1–G6 hold, four branches from the source audit remain open — accept the
horizon, ghost, quantum stress, reinterpretation — and the honest position is
that the classical-geometry escape has now been checked and closed for the
natural `D = 5` higher-curvature term.

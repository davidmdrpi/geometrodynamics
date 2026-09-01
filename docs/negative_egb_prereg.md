# Pre-registration: does negative-coupling EGB actually work?

**Status: frozen before the module exists.** Committed alone, as with the three
previous pre-registrations.

---

## The branch this opens

PR #277 left negative-coupling Einstein–Gauss–Bonnet as the one branch
**narrowed but not closed**: `α_GB ≤ −R²/4` does satisfy the matter NEC along
the throat, and the ledger recorded *"global existence and stability are open"*.
This round asks the first half of that.

The key point the previous round did not use: **`α_GB` is a coupling constant in
the action, so the same value acts in the exterior.** The throat was analysed in
isolation. It should not have been.

---

## Verified before freezing

For the ultrastatic exterior `−dt² + R²(dχ² + sin²χ dΩ₃²)` — that is,
`ℝ × S⁴_R` — with `k = (1, 1/R, 0,0,0)` verified null, computed symbolically and
confirmed against a numerically differentiated Riemann tensor:

| quantity | value |
|--|--|
| `R_scalar` | `12/R²` |
| `L_GB` | `24/R⁴` |
| `R_kk` | `3/R²` |
| `H_kk` | `12/R⁴` |
| `H^t_t` | `−12/R⁴` |
| `H^i_j` | `0` |

*A diagnostic error worth recording:* a first pass mis-raised the leading index
when forming the scalar `R_abcd R^abcd`, giving `L_GB = 8.177` against the exact
`24/R⁴ = 8.403`. The tell was that **the error did not shrink under refinement**
— stuck at `2.26e-01` for `h = 4e-3 … 5e-4` — so it was structural, not
discretisation. The published `gauss_bonnet` module is unaffected: `g_ab L_GB`
drops under null contraction, and its `term4` contracts both leading indices
with `kᵃkᵇ`, which correctly keeps them down.

---

## The six frozen predictions

**E1 — the exterior obeys the same ratio law.** `H_kk/R_kk = 4μ/f⁴` holds
outside as well as inside, and in the exterior `μ/f⁴ = R²sin⁴χ/(R⁴sin⁴χ) =
1/R²` — **independent of `χ`**. So `H_kk/R_kk = 4/R²` throughout the exterior.

**E2 — the exterior's NEC constrains `α` in the *opposite* direction.**

```
8πG₅ T_kk^ext = R_kk + α_GB H_kk = 3(R² + 4α_GB)/R⁴  ≥ 0   ⟺   α_GB ≥ −R²/4
```

**E3 — the two constraints meet at a single point.** The throat needs
`α_GB ≤ −R²/4`; the exterior needs `α_GB ≥ −R²/4`. The only value satisfying
both is

```
α_GB = −R²/4        exactly
```

a **measure-zero** solution. There is no open set of couplings for which this
spacetime satisfies the matter NEC everywhere.

**E4 — the mechanism, and why it is not a coincidence.** Both regions share the
identical bracket, since `T_kk = R_kk(1 + 4α μ/f⁴)`:

- in the throat `R_kk < 0` (flare-out), so the NEC needs the bracket `≤ 0`;
- in the exterior `R_kk > 0` (ordinary matter), so it needs the bracket `≥ 0`.

**Same bracket, opposite required signs.** And `μ/f⁴` is continuous across the
seam at `1/R²` — the *same* Misner–Sharp continuity as P2 — so the bracket is
continuous there too and must therefore be exactly zero. The throat's binding
constraint is at the mouth, which is the seam; the exterior's is everywhere.
They meet because they are evaluated on the same continuous function.

**E5 — the exterior pressure is untouched by Gauss–Bonnet.** `H^i_j = 0` for a
maximally symmetric spatial slice, so

```
8πG₅ ρ_ext = 6/R² + 12α_GB/R⁴          8πG₅ p_ext = −3/R²   (α-independent)
```

**E6 — at the critical coupling the observed universe must be empty.** At
`α_GB = −R²/4`:

```
8πG₅ ρ_ext = 3/R² ,    8πG₅ p_ext = −3/R² ,    w = p/ρ = −1   exactly
```

The exterior is forced to be **pure vacuum energy with no ordinary matter at
all**, at half the Einstein-gravity density.

---

## ⚠ Erratum (post-review): E6's interpretation is withdrawn

Three claims below overreached, and the branch closes for a different reason.

**"The exotic matter is merely relocated" — false at the critical point.** The
throat stress there satisfies NEC, WEC **and** DEC. With `A = b²/f⁴` and
`q = R²b²/f⁴ ≥ 1`: `ρ = 3Aq`, `p_s = −3A`, `p_Ω = A`, so `ρ+p_s = 3A(q−1) ≥ 0`
(saturated only at the mouth), `ρ+p_Ω = A(3q+1) > 0`, and `ρ ≥ |p_s|, |p_Ω|`.
Relocation happens only for `α` strictly below critical.

**"The observed universe must be empty" — overreached.** `w = −1` describes the
*total 5D bulk stress* supporting `ℝ × S⁴_R`, not the observed `S³`, which is
its equator and carries a different stress tensor; and a homogeneous `−Λg_ab`
can be moved to the gravitational side. The defensible statement is that the
branch forces a **vacuum-energy-like 5D exterior**.

**"Measure zero, therefore refuted" — not by itself.** The relation is equally
read as the field equations selecting `R² = −4α_GB`, and gravity routinely ties
a vacuum curvature radius to a coupling. Calling it tuning needs an
independently fixed `R`.

**What closes the branch instead: the graviton.** Linearising the full
`G_ab + α H_ab` on this background for a TT mode gives

```
C_kin = −½(1 + 4α_GB/R²)        c² = 1/(1 + 4α_GB/R²)
```

*derived* for this product spacetime, since the familiar
`1 + 2α(D−3)(D−4)K` is for a maximally symmetric one. It vanishes **exactly** at
`α = −R²/4` — the same value the NEC pins — while the spatial `κ²` coefficient
is entirely coupling-independent. So the tensor cone is **superluminal for every
`−R²/4 < α < 0`**, diverges at criticality, and there the `ω²` term disappears
while `κ²` survives: the equation degenerates from evolution to constraint. The
mass term stays finite, so it is the kinetic operator that dies.

**"Global existence" is also too strong.** This determines the stress the metric
requires; it does not exhibit fields obeying their own equations that supply the
throat's anisotropic stress.

**E5's reason was wrong, and the true statement is stronger.** E5 above, and
three docstrings, attributed `H^i_j = 0` to the `S⁴_R` slice being *maximally
symmetric*. That is the right value for the wrong reason. In `D = 5` the spatial
block of `H^a_b` for **any** ultrastatic product `−dt² + h₄` is the
four-dimensional Gauss–Bonnet (Euler) tensor of `h₄`, which vanishes identically
because Gauss–Bonnet is topological in `D = 4`. No symmetry is used. Checked in
`gauss_bonnet.measure_the_spatial_block_vanishes` on the throat slice (which is
*not* maximally symmetric) and on a generic lumpy 4-slice with no symmetry at
all — both converging to zero at ratio 16 under a 4× step reduction — against a
nonconstant-lapse control that stays at `2.3e-2` and does **not** shrink.

The consequence is not cosmetic: Gauss–Bonnet cannot touch the pressures
*anywhere* on this geometry, throat as well as exterior, so the entire
`α`-dependence lands in the density everywhere. That is what makes the exact
throat stress `8πG p_s = −3b²/f⁴`, `8πG p_Ω = +b²/f⁴` the Einstein ones at every
coupling, and `H_tt = −12b⁴/f⁸` — both confirmed numerically to `1e-4`.

---

## The frozen verdict

**The negative-coupling branch closes on physical grounds, not on a regime
complaint.** It survives only at one exact value of a continuous coupling, and
at that value BAM's exterior — which is supposed to *be* the observed closed
universe — contains nothing but a cosmological constant. Push `α_GB` any more
negative and the **exterior** violates the NEC, so the exotic matter is merely
relocated from the throat to the universe around it.

---

## Falsifiers

1. An exterior `H_kk` disagreeing with `12/R⁴`, or `H^i_j ≠ 0`, falsifies E1/E5.
2. A `χ`-dependent `μ/f⁴` in the exterior falsifies E1 and would decouple the
   two constraints — the most consequential possible outcome, since it would
   reopen an interval of `α`.
3. Any `α_GB` making **both** `T_kk^throat ≥ 0` along the throat and
   `T_kk^ext ≥ 0` falsifies E3. A scan must look for one rather than assume
   none.
4. A discontinuity of the bracket `1 + 4αμ/f⁴` across the seam falsifies E4.
5. `w ≠ −1` at the critical coupling falsifies E6.

---

## Scope

This settles **global existence** for constant-coupling EGB on *this* geometry:
there is no open set of `α_GB`, and the single surviving value empties the
universe. It does **not** address stability or the graviton kinetic term at that
coupling — the second half of the ledger's open item — and it does not touch
dilatonic `α(φ)L_GB` or `f(R)`, where the scalar's own stress enters. Nor does
it rule out a *different* exterior: the constraint is derived for the round
`S⁴_R` completion this programme assumes.

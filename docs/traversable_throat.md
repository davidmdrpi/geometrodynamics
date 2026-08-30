# The whole-throat S-matrix of a supported traversable 5D throat

*Modules `geometrodynamics/tangherlini/traversable_throat.py` and
`geometrodynamics/transaction/derived_network.py`; probes
`traversable_throat_probe.py` (7/7) and `derived_network_probe.py` (6/6); tests
`tests/test_traversable_throat.py` (31/31) and `tests/test_derived_network.py`
(22/22).*

---

## What this round is for

`transaction/network.py` (PR #216) did something genuinely good: it replaced
`handshake.py`'s advanced confirmation — which its own docstring calls
*"phenomenological … a conjugated ('advanced') phase assigned by fiat"* — with a
Morris–Thorne–Yurtsever mechanism. Every wave propagates forward in local time;
the mouths' clock offset `Δ_BA` carries the return into the exterior past. The
postulate moved from a bare nonlocal rule to a geometric one.

But the geometry was never supplied. `MouthPort.t`, `r_in`, `r_out` are inputs,
`tau_th` is a constant, and

```python
def closure_offset(d_A, d_B, tau_th):
    """The clock offset Delta_BA for which the network path returns
    exactly to the emission event."""
    return -(d_A + d_B + tau_th)
```

*solves* for the offset that makes the loop close. This round supplies the
missing physical input and recomputes the closure without tuning `Δ` to the
answer.

---

## Three results

**1. Traversability has an exact price.**

```
∫ T_ab k^a k^b dλ  =  −3/(16 G₅ a)
```

along a complete radial null geodesic — exact, and scaling as `1/(G₅a)`.

**2. No finite frequency closes carrier and packet together.** The branch-free
residual `C_ℓ(ω) = Arg exp[i(θ_ℓ − ωθ_ℓ′)]`, with `θ = arg(η_topo T_ℓ)` in
`network.py`'s own convention, has **no root** over `[0.2, 12]` at 900 points.
Its decay is analytic: `ωC → −∫V_ℓ ds = −9π/8`, so `C` vanishes only as `1/ω`.

**3. No finite-frequency perfect-transmission point was found** — by direct
search for `R_ℓ(ω) = 0`, stated as a finding and *not* as a theorem, since a
positive barrier can support such resonances — and the deficit closes as
`1 − |T|² ~ e^{−4aω}`, the first-Born prediction with nothing fitted.

---

## Why not Tangherlini — an argument, not a computation

The cross-exterior channel on a maximally extended bridge vanishes for a reason
needing no numerics:

```
supp G_ret(c,s) ⊂ J⁺(s)      supp G_adv(c,d) ⊂ J⁻(d)
```

so `G_ret(c,s)·G_adv(c,d) ≠ 0` requires `c ∈ J⁺(s) ∩ J⁻(d)`. But then
`s → c → d` is a causal chain, so `d ∈ J⁺(s)`. Contrapositive:

> **If `d ∉ J⁺(s)`, the product vanishes for every `c`.**

**An advanced leg by itself does not evade ER non-traversability.** PR #275's
`T_ℓ` is exterior→*horizon*, not exterior→exterior — exactly the zero this
argument predicts for the cross-mouth channel. Computing `G_sym = ½(G_ret+G_adv)`
would produce another zero and arbitrate nothing.

So the network mechanism needs a genuinely traversable throat. That is a choice,
and this round prices it.

---

## The geometry, derived

`ds² = −dt² + ds² + (s²+a²)dΩ₃²`, with `f = √(s²+a²)` the scalar-flat solution
of `f′² = 1 − a²/f²`. With `Φ = e^{−iωt}Y_ℓ(Ω)u(s)` and `ψ = f^{3/2}u`, removing
the first-derivative term contributes `(3/4)(f′/f)² + (3/2)f″/f`, giving

```
ψ″ + [ω² − V_ℓ(s)]ψ = 0 ,   V_ℓ = ℓ(ℓ+2)/f² + (3/4)(s² + 2a²)/f⁴
```

| `ℓ` | `V_ℓ(0)` | exact `[ℓ(ℓ+2)+3/2]/a²` | tail coeff | exact `(ℓ+1)²−¼` |
|--|--|--|--|--|
| 0 | `1.500000` | `1.500000` | `0.750000` | `0.750000` |
| 1 | `4.500000` | `4.500000` | `3.750000` | `3.750000` |
| 2 | `9.500000` | `9.500000` | `8.750000` | `8.750000` |
| 3 | `16.500000` | `16.500000` | `15.750000` | `15.750000` |

Two structural facts, both load-bearing:

**It is one barrier, not two.** `V` is a single smooth symmetric maximum at
`s = 0`, positive everywhere, strictly decreasing outward. The Fabry–Pérot form
`t_A t_B/(1 − r_inA r_inB e^{2iωτ})` that PR #216 assumes is a *modelling
choice the geometry does not supply*, and this module does not impose it — the
whole throat is solved as one scattering problem over `s ∈ (−∞,∞)`. If a later
geometry really does contain two separated interfaces, the factorisation becomes
derived rather than assumed.

**The tail is PR #275's.** `V → [(ℓ+1)² − ¼]/s²` at *both* ends — numerically
the same centrifugal form as the Tangherlini exterior — so that round's
validated Riccati–Hankel Jost condition (order `ν = ℓ+1`, normalised to
`e^{±ix}`) is imported rather than rewritten. Plane waves would be the wrong
basis at low `ω` for exactly the same reason as there.

*A small echo: `V_ℓ(0) = [ℓ(ℓ+2)+3/2]/a²` carries the same `ℓ(ℓ+2)+3/2` that
#275 derived as the exact `∫V dr*` on Tangherlini.*

---

## The GR price

From the 5D Einstein tensor in an orthonormal frame — and `R₅ = 0`, though the
spacetime is not vacuum:

| `s` | `8πG₅ρ` | `8πG₅p_s` | `8πG₅p_Ω` |
|--|--|--|--|
| `0` | `0.000000` | `−3.000000` | `+1.000000` |
| `0.5` | `0.000000` | `−1.920000` | `+0.640000` |
| `1` | `0.000000` | `−0.750000` | `+0.250000` |
| `2` | `0.000000` | `−0.120000` | `+0.040000` |

```
ρ = 0    exactly
p_s = −3a²/(8πG₅ f⁴)      p_Ω = +a²/(8πG₅ f⁴)   (each of three)
```

The radial NEC `ρ + p_s < 0` fails **everywhere**, and the complete null-ray
integral is exactly `−3/(16G₅a)`. Exoticity scales as `1/(G₅a)`: a wider throat
is cheaper, which is the only knob available.

> **This is not the cost of the clock offset.** It is the support required by
> the *static* traversable throat. A frozen final metric can carry a large
> accumulated `Δ_BA` with no local stress tensor proportional to it — the offset
> is produced by the mouths' differential-aging *history*. Costing `Δ` needs
> moving-mouth dynamics, and this round does not attempt it.

---

## The scattering matrix

Piecewise-constant transfer matrix vectorised over `ω`, the method #275
validated, with Jost matching at both ends.

| `ω` | `\|T\|` | `\|R\|` | `\|R\|²+\|T\|²−1` |
|--|--|--|--|
| `0.05` | `0.00098705` | `0.99999951` | `−1.6e-14` |
| `0.20` | `0.01665893` | `0.99986123` | `+1.0e-14` |
| `1.00` | `0.64142490` | `0.76718583` | `−1.4e-14` |
| `2.00` | `0.99657150` | `0.08273606` | `−3.7e-14` |
| `10.0` | `1.00000000` | `0.00000000` | `−6.5e-14` |

Unitarity to `6.6e-14`, imposed nowhere. **Reciprocity is structural**: `V` is
even in `s` to machine precision, so `R_L = R_R` and `T_LR = T_RL` identically —
a property of the geometry, not a numerical coincidence.

---

## The threshold law, and the false regression avoided

The exact static monopole response has `I₃ = ∫ds/f³ = 2/a²` and
`g = 2π²/I₃ = π²a²`. It is tempting to regress `T₀(0)` against `g`. **That would
be wrong.** A flux-normalised transmission amplitude between asymptotically flat
four-dimensional spatial ends vanishes at threshold, from the asymptotic radial
normalisation. What the static solution controls is the *coefficient* of the
low-frequency expansion:

```
|T₀(ω)| → (π/8)(aω)²
```

| `ω` | `\|T₀\|/(π/8)(aω)²` | `arg(ratio)/π` |
|--|--|--|
| `0.010` | `1.000294` | `+0.5000` |
| `0.015` | `1.000618` | `+0.5001` |
| `0.020` | `1.001043` | `+0.5001` |
| `0.030` | `1.002168` | `+0.5002` |
| `0.050` | `1.005398` | `+0.5006` |

The magnitude law is confirmed and converges. The phase of the ratio is a
**constant** `+π/2` to within `5e-3` across the band — i.e. exactly a factor of
`i`, from this module's Riccati–Hankel normalisation to `e^{±ix}`. The oracle
was stated up to the asymptotic phase convention, and that is precisely what the
residual is.

Static conductance and dynamical threshold law are two faces of one calculation
— but they are not the same number, and the module says so.

---

## Wiring it into `network.py`, not beside it

`T_ℓ(ω)` is only worth having if the module that assumed it actually runs on it.
`NetworkThroat` therefore gained a **second backend** — an optional callable
`whole_throat_transfer` — beside the existing `MouthPort` one, which is retained
untouched. `NetworkThroat.transfer(ω)` dispatches to whichever is present, and

```python
def derived_loop_eigenvalue(throat, w, d_A, d_B, delta):
    return throat.topological_factor * throat.transfer(w) \
        * np.exp(1j * w * (d_A + d_B + delta))
```

forms `Λ_ℓ(ω, Δ) = η_topo · T_ℓ(ω) · e^{iω(d_A + d_B + Δ)}` with
`η_topo = NetworkThroat.topological_factor`, the module's own deck orientations
and mouth phases — `(−1)^ℓ`-type parity and `η_M` kept as the separate
operations they are, not collapsed into one sign.

No fake `MouthPort`s are manufactured to fit the old API: the ports on a derived
throat are transparent, and the transfer comes entirely from the metric.

**There is deliberately no `τ_th` phase in this path.** A whole-throat `T`
already carries the transit in `arg T`; adding `τ_glob` on top would
double-count the Wigner delay. The Fabry–Pérot backend needs it only because its
`t_AB` is an *excess* factor over free interior propagation — a different object.

Two consistency checks, both structural rather than tuned. `|η_topo| = 1`, so
`|Λ| = |T|`, confirmed to `< 1e-12` at every probe frequency. And the batched
scan used for the continuous searches below is asserted equal to the scalar
`network.py` path to `1e-12`, so the searches interrogate the object the module
actually exposes.

This matters for more than tidiness: the closure residual defined next **shifts
under a constant rephasing of the transfer**, so it is only well posed inside
the network's own convention, with its own `η_topo`.

---

## No finite frequency closes carrier and packet together

Once the throat is dispersive there is no unique `τ_th`. The geometry supplies
`δ_ℓ(ω) = arg T_ℓ(ω)`, and closure demands

```
phase:   ω(D + Δ) + θ_ℓ = 2πn
group:   D + Δ + θ_ℓ′   = 0
```

**Comparing the two offsets directly is not an invariant test.** An earlier
draft of this round did exactly that at `n = 0` with `φ_topo = 0` and reported
gaps of `6.78, 4.14, …`. Phase branches are `2π/ω` apart, so at `ω = 1` a raw
gap of `4.14` is `2.14` to the nearest branch; and a constant rephasing of `T`
shifts `δ/ω` — hence `Δ_phase` — while leaving the Wigner delay untouched. The
raw dispersion table is retained in the module, flagged as data rather than
verdict.

Eliminating `Δ` removes both problems:

```
θ_ℓ − ω θ_ℓ′ = 2πn   ⟺   C_ℓ(ω) ≡ Arg exp[i(θ_ℓ − ω θ_ℓ′)] = 0
```

`C` searches over `n` automatically. It is still shifted by a constant
rephasing of the transfer, which is why it must be computed **end-to-end** with
`network.py`'s own `η_topo` rather than assembled from a bare `T`.

| `ω` | `C_ℓ(ω)` |
|--|--|
| `0.20` | `+3.10855` |
| `2.17` | `−1.82070` |
| `4.14` | `−0.87519` |
| `6.11` | `−0.58503` |
| `8.08` | `−0.44033` |
| `10.04` | `−0.35325` |

**No roots** over `[0.2, 12]` at 900 points; smallest `|C| = 0.295`.

### And the decay is analytic

At high frequency `θ ≈ −c_ℓ/ω` with `c_ℓ = ½∫V_ℓ ds`, so `θ′ ≈ +c_ℓ/ω²` and

```
C_ℓ = θ − ωθ′ ≈ −2c_ℓ/ω = −∫V_ℓ ds / ω = −(π/a)[ℓ(ℓ+2) + 9/8] / ω
```

| `ω` | `C` | `ωC` |
|--|--|--|
| `4` | `−0.906983` | `−3.627931` |
| `8` | `−0.444534` | `−3.556271` |
| `12` | `−0.295330` | `−3.543958` |
| `20` | `−0.176888` | `−3.537753` |

Predicted `−9π/8 = −3.534292`, reached to `0.1%` by `ω = 20`, monotonically and
with nothing fitted.

> **`C` vanishes as `1/ω` and no faster.** Simultaneous carrier-and-packet
> closure is a UV limit, never attained at finite frequency — and it is *the
> same limit* in which `|T| → 1`. The loop's magnitude and both of its phase
> conditions all close only in the ultraviolet.

---

## Can the transaction complete?

PR #216 carries `Λ_ℓ(ω) = t_net η_topo e^{iωD_loop}` and `G_eff = G₀/(1−Λ)`,
calling `Λ → 1` the completed transaction. Wired to the derived geometry,
`Λ_ℓ(ω,Δ) = η_topo T_ℓ(ω) e^{iω(d_A+d_B+Δ)}` — with **no** `tau_th` phase,
since `arg T` already carries the transit and adding one would double-count the
Wigner delay. `|η_topo| = 1`, so `|Λ| = |T|` to `1e-12`.

**A direct search, not a theorem.** A positive barrier does *not* in general
forbid perfect transmission — positive barriers can have transmission
resonances — so `|T| < 1` does not follow from `V > 0`. Searching for zeros of
`R_ℓ(ω)` over `[0.05, 12]` at 1200 points finds **no interior minimum**: `|R|`
falls monotonically and its smallest value sits at the top of the range. So no
finite-frequency perfect-transmission point was found, and where `|T| < 1` we
have `|Λ| < 1` and `G_eff` has no pole.

**The falloff is Born, not fitted.** For `ℓ = 0`,
`Ṽ₀(q) = (3π/8a)(3 + a|q|)e^{−a|q|}`, so first-Born reflection at momentum
transfer `2ω` gives `1 − |T|² ~ e^{−4aω}`: slope `−4a`, nothing fitted. The
*local* slope descends monotonically toward it —

| `ω` | local slope |
|--|--|
| `2.00` | `−4.7188` |
| `3.00` | `−4.4287` |
| `4.00` | `−4.2678` |
| `5.00` | `−4.1844` |
| `6.00` | `−4.1351` |

— so the `−4.25` an earlier draft quoted from a global fit was the
finite-frequency approach to the analytic asymptote, not a new constant.

| `D_loop` | regime |
|--|--|
| `> 0` | ordinary delayed feedback |
| `= 0` | closed-null / marginal chronology condition |
| `< 0` | return before emission (CTC regime) |

---

## Ledger

| claim | verdict | evidence |
|--|--|--|
| an advanced leg alone evades ER non-traversability | **NO — PROVED, NOT COMPUTED** | support composition forces `d ∈ J⁺(s)` |
| PR #216's throat transfer is supplied by a geometry | **NOW IT IS** | `T_ℓ`, `R_ℓ` from the whole throat instead of port inputs |
| PR #216's loop can be driven by that geometry | **YES, END TO END** | `NetworkThroat.whole_throat_transfer`; `Λ` from `network.derived_loop_eigenvalue` with the module's own `η_topo`; `MouthPort` backend retained |
| the derived loop needs a separate `τ_th` phase | **NO — IT WOULD DOUBLE-COUNT** | `arg T` already carries the transit; `τ_glob` is an excess-factor artefact of the Fabry–Pérot backend |
| comparing `Δ_phase` at `n = 0` to `Δ_group` is the test | **NO — BRANCH DEPENDENT** | branches are `2π/ω` apart; the `4.14` raw gap at `ω = 1` is `2.14` to the nearest branch |
| the residual's `1/ω` decay is a fitted observation | **NO — IT IS ANALYTIC** | `ωC → −(π/a)[ℓ(ℓ+2)+9/8] = −9π/8`, matched to `0.1 %` by `ω = 20` |
| the UV falloff constant `−4.25` is a new constant | **NO — IT IS BORN, AND IT IS `−4a`** | the *local* slope descends `−4.72, −4.43, −4.27, −4.18, −4.13` toward the predicted `−4` |
| the throat factorises as two mouths plus a cavity | **NOT BY THIS GEOMETRY** | one smooth symmetric barrier; not imposed here |
| the static conductance `g` is `T₀(0)` | **NO — IT IS THE THRESHOLD COEFFICIENT** | `\|T₀\| → (π/8)(aω)²`, ratio `1.0003` at `ω = 0.01` |
| the traversable throat is free | **NO** | NEC violated everywhere; null-ray integral exactly `−3/(16G₅a)` |
| that NEC integral is the cost of the clock offset | **NO — DIFFERENT QUESTION** | static support vs aging history |
| one clock offset closes the loop | **NOT AT ANY FINITE `ω`** | branch-free `C` has no root over `[0.2,12]`; `ωC → −9π/8` so `C ~ 1/ω` |
| `Λ = 1` can occur | **NO SUCH POINT FOUND — NOT A THEOREM** | direct search for `R(ω)=0` finds no interior zero; deficit `~ e^{−4aω}` from Born |

---

## What this round establishes, and what it does not

PR #216 relocated the advanced wave from a bare postulate into a
Morris–Thorne–Yurtsever geometry. That was real progress. This round supplies
the geometry it was missing, **drives PR #216's own loop with it**, and finds
three things: the relocation has an exact price, the closure it assumed cannot
be met at any finite frequency on a branch-free test, and its completed
transaction is a high-`Q` UV limit rather than an attainable pole. Because the
wiring is end-to-end, those are statements about the BAM module itself rather
than about a reconstruction beside it.

**Two claims were weakened in the writing.** The `Δ_phase`/`Δ_group` gap was
branch dependent and is withdrawn as a verdict — the invariant statement needed
eliminating `Δ` first. And `|T| < 1` at finite `ω` does **not** follow from
`V > 0`: positive barriers can have perfect-transmission resonances, so what
stands is a search that found none, not a theorem. In both cases the stronger
sentence was available and wrong.

**Scope, kept explicit.** The benchmark has two asymptotically flat ends at
`s → ±∞`, while `network.py` conceptually has two *finite* mouths embedded in
the closed `S³` exterior. `T_ℓ` is a whole-throat oracle, not a literal glued
finite-mouth solution; its UV normalisation `T → 1` is what makes it usable as
an excess transfer factor. Finite matching surfaces to the `S³` exterior, and
their junction stress, are a later construction — and should not be introduced
merely to make the old `MouthPort` API fit.

**What remains postulated:** that the throat is traversable at all, and that the
mouths carry the required clock offset. The first now has an exact price. The
second needs moving-mouth dynamics this round does not do.

**Still open:** the history that produces `Δ_BA`, and whether any BAM exchange
kernel is meant to approximate the whole-throat `T_ℓ` derived here.

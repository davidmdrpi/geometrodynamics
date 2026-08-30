# The whole-throat S-matrix of a supported traversable 5D throat

*Modules `geometrodynamics/tangherlini/traversable_throat.py` and
`geometrodynamics/transaction/derived_network.py`; probes
`traversable_throat_probe.py` (7/7) and `derived_network_probe.py` (8/8); tests
`tests/test_traversable_throat.py` (31/31) and `tests/test_derived_network.py`
(35/35).*

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

**2. On the topology the repository declares, the loop closes at `ω = 1.4617`.**
The branch-free closure function `Ψ_ℓ(ω) = θ_ℓ − ωθ_ℓ′` with
`θ = arg(η_topo T_ℓ)` reaches `2πn` there, stable to `1e-6` under refinement of
the matching edge, the spatial step and the difference step. **But whether it
closes is gauge dependent**: `Ψ` sweeps `3.9676 < 2π`, so a constant rephasing
can create or remove the root. What is invariant is `dΨ/dω = −ωθ″` and that
total variation.

**3. No perfect-transmission point was found on the tested band `[0.05, 12]`** —
by direct search for `R_ℓ(ω) = 0`, stated as a band-limited finding and *not* as
a theorem, since a positive barrier can support such resonances — and the
deficit closes as `1 − |T|² ~ e^{−4aω}`, the first-Born prediction with nothing
fitted.

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
untouched.

**The dispatch lives in `t_AB`,** the one primitive that `traverse_throat`,
`network_confirmation`, `projected_kernel`, `loop_eigenvalue` and
`effective_green` all already call. That placement is the whole point. A first
draft of this round dispatched only inside a *new* `derived_loop_eigenvalue`,
which left every one of those five entry points reading the transparent ports —
so `effective_green(derived_throat, …)` was **not** the `G₀/(1−Λ)` the round was
discussing, and the object's behaviour depended on which API a caller happened
to reach for. Measured agreement now:

| API | disagreement with the derived `T` |
|--|--|
| `t_AB` | `0` |
| `traverse_throat` | `0` |
| `loop_eigenvalue` vs `derived_loop_eigenvalue` | `0` |
| `effective_green` | `5.6e-17` |

`derived_loop_eigenvalue` is a convenience spelling for scans that sweep `Δ`,
not a second code path — with `Δ = throat.delta_BA` it returns exactly what
`loop_eigenvalue` returns, and that is asserted rather than assumed.

**The double count is made unconstructible, not merely avoided.**
`NetworkThroat.__post_init__` rejects a derived backend with `tau_th ≠ 0`, so
the traversal leg's free transit phase `e^{−iωτ_th}` is exactly `1` and `arg T`
is counted once. `loop_expansion` and `r_AA` raise rather than answer from a
transparent port: a smooth single barrier has no echo train, and that
decomposition is a property of the two-interface model, not of the throat.

No fake `MouthPort`s are manufactured to fit the old API: the ports on a derived
throat are transparent, and the transfer comes entirely from the metric. The
mouths themselves come from `make_singlet_pair()` — so `η_topo` is derived, and
orientation, wrap parity, the closure phase and the transfer stay the four
separate operations they are rather than collapsing into one sign.

`|η_topo| = 1`, so `|Λ| = |T|` to `< 1e-12`; and the batched scan used below is
asserted equal to the scalar `network.py` path to `1e-12`.

---

## Where the loop closes — and how much of that is gauge

Once the throat is dispersive there is no unique `τ_th`. The geometry supplies
`δ_ℓ(ω) = arg T_ℓ(ω)`, and closure demands

```
phase:   ω(D + Δ) + θ_ℓ = 2πn
group:   D + Δ + θ_ℓ′   = 0
```

The second is *group-delay closure at the carrier* — the necessary first-order
condition for a finite-band packet to return to the emission event, not exact
packet closure, which would also constrain the amplitude variation and every
higher phase derivative across the band.

**Comparing the two offsets directly is not an invariant test.** An earlier
draft of this round did exactly that at `n = 0` with `φ_topo = 0` and reported
gaps of `6.78, 4.14, …`. Phase branches are `2π/ω` apart, so at `ω = 1` a raw
gap of `4.14` is `2.14` to the nearest branch; and a constant rephasing of `T`
shifts `δ/ω` — hence `Δ_phase` — while leaving the Wigner delay untouched. The
raw dispersion table is retained in the module, flagged as data rather than
verdict.

Eliminating `Δ` removes both problems:

```
Ψ_ℓ(ω) ≡ θ_ℓ − ω θ_ℓ′ = 2πn ,     θ_ℓ = arg(η_topo T_ℓ)
```

`Ψ` searches over `n` automatically.

### The sign of `η_topo` is not free, and it decides the answer

`ConjugatePair.__post_init__` **asserts** that the two mouths of one throat
carry opposite orientations, and `make_singlet_pair` builds `(+1, −1)`. So
`η_topo = orientation_A · orientation_B · e^{i(φ_A+φ_B)}` contains a factor
`−1` for any real BAM throat. A first draft of this round used `(+1, +1)` — a
*chosen* pair — and that choice, not the geometry, produced its headline.

`network_mouth_from_defect` now maps `ThroatDefect → NetworkMouth`, so the
orientations are read off the repository's own construction. And the wrap
parity is kept **separate**: `ThroatDefect.spinor_sign()` is the Hopf-holonomy
sign a *spinor* acquires, which a scalar does not see. The two products both
equal `−1`, so they cancel:

| channel | `η_topo` | roots of `Ψ = 2πn` on `[0.05, 20]` |
|--|--|--|
| scalar (the field solved here) | `−1` | **`ω = 1.4617`** |
| spinor (orientation × wrap) | `+1` | none |

Collapsing those two operations into one sign would have hidden a genuine
disagreement, so both are computed.

Roots are found two ways, because either alone is incomplete: sign changes of
`Ψ − 2πn` are bracketed with Brent, and every interior local minimum of the
smooth objective `|e^{iΨ} − 1|` is refined with a bounded minimiser — the
second is what catches a *tangential* zero that touches a branch without
crossing it.

### The phase needs its own convergence study

`|R|² + |T|² = 1` holds to `1e-13` no matter how badly `arg T` is resolved:
unitarity constrains moduli. `Ψ` differentiates `arg T`, so the root is
re-measured against all three knobs that could move it.

| variant | root |
|--|--|
| baseline (`edge 200`, `60000` steps, `fd 1e-4`) | `1.461703899` |
| `edge 150` / `edge 300` | `1.461703753` / `1.461704316` |
| `30000` / `120000` steps | `1.461704898` / `1.461703649` |
| `fd 1e-3` / `fd 1e-5` | `1.461704151` / `1.461703897` |

Worst shift `1.0e-6`, so the root is quoted as **`ω = 1.4617`**.

### What a rephasing can move, and what it cannot

`Ψ` is *not* invariant under a constant rephasing of the Jost basis: `θ → θ+c`
sends `Ψ → Ψ+c`. And over the whole band `Ψ` sweeps only

```
Ψ ∈ [−1.0029, 2.9647] ,   total variation 3.9676  <  2π = 6.2832
```

Because the swing is **less than one full branch**, that constant can create or
destroy the root. So neither *"it closes"* nor *"it never closes"* is a
property of the geometry by itself.

Two things survive:

```
dΨ/dω = θ′ − θ′ − ωθ″ = −ω θ″        (no constant in it)
```

verified to `2.2e-5` with clean second-order convergence (`5.6e-4 → 1.4e-4 →
2.2e-5`), together with the total variation. A *linear* reference-plane phase
`b·ω` is harmless for a separate reason — it cancels identically from
`θ − ωθ′` — so the entire residual freedom is one constant, not a function.

Part of that constant is now derived: the topological sign, from
`ConjugatePair`. The rest — the Jost basis constant — is still a software
convention, and fixing it physically needs the finite-mouth matching surfaces
this round does not build. **The closure verdict is therefore stated relative
to that basis, and deliberately not promoted to a basis-free claim about BAM.**

### And the decay is analytic — to `arg η_topo`, not to zero

At high frequency the *bare* eikonal phase is `arg T ≈ −c_ℓ/ω` with
`c_ℓ = ½∫V_ℓ ds`, so

```
Ψ_ℓ − arg η_topo ≈ −2c_ℓ/ω = −∫V_ℓ ds / ω = −(π/a)[ℓ(ℓ+2) + 9/8] / ω
```

**The topological constant has to come out first.** `θ = arg(η_topo T)` carries
`arg η_topo` at every frequency, so `ωΨ` would diverge linearly instead of
tending to a limit — an omission that is invisible at `η_topo = +1`, where the
constant is zero, and which the first draft therefore did not notice.

| `ω` | `Ψ` | `Ψ − π` | `ω(Ψ − π)` |
|--|--|--|--|
| `4` | `+2.234610` | `−0.906983` | `−3.627931` |
| `8` | `+2.697059` | `−0.444534` | `−3.556271` |
| `12` | `+2.846263` | `−0.295330` | `−3.543958` |
| `20` | `+2.964705` | `−0.176888` | `−3.537753` |

Predicted `−9π/8 = −3.534292`, reached to `0.10%` by `ω = 20`, monotonically
and with nothing fitted.

> **`Ψ` tends to `π`, not to `0`** — the furthest a phase can sit from a branch
> `2πn`. The loop is *least* closed in the ultraviolet, and the finite root at
> `ω = 1.4617` is a crossing on the way there. An earlier draft read the `1/ω`
> law as proving no finite root existed; an asymptotic law constrains the tail,
> not the interior.

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
| PR #216's loop can be driven by that geometry | **YES — THROUGH THE EXISTING APIs** | dispatch is in `t_AB`, so `traverse_throat`, `network_confirmation`, `projected_kernel`, `loop_eigenvalue`, `effective_green` all see the derived `T` |
| the derived loop needs a separate `τ_th` phase | **NO — AND IT CANNOT CARRY ONE** | `arg T` already carries the transit; `__post_init__` *rejects* a derived throat with `τ_th ≠ 0` |
| `η_topo` may be chosen as `+1` | **NO — IT IS DERIVED, AND IT IS `−1`** | `ConjugatePair` asserts opposite mouth orientations; `make_singlet_pair` builds `(+1,−1)` |
| comparing `Δ_phase` at `n = 0` to `Δ_group` is the test | **NO — BRANCH DEPENDENT** | branches are `2π/ω` apart; the `4.14` raw gap at `ω = 1` is `2.14` to the nearest branch |
| the closure verdict is a property of the geometry | **NO — IT IS GAUGE DEPENDENT** | `Ψ` sweeps `3.9676 < 2π`, so a constant rephasing can create or remove the root |
| nothing about the closure survives a rephasing | **TWO THINGS DO** | `dΨ/dω = −ωθ″` (to `2.2e-5`, second order) and the total variation; a *linear* `b·ω` cancels identically |
| the scalar and spinor channels agree | **NO** | scalar `η_topo = −1` closes; spinor also carries `spinor_sign`, giving `+1`, and does not |
| `Ψ` decays to zero as `1/ω` | **NO — TO `arg η_topo`** | `ω(Ψ − π) → −(π/a)[ℓ(ℓ+2)+9/8] = −9π/8`, matched to `0.10 %` by `ω = 20`; unsubtracted, `ωΨ` diverges |
| that `1/ω` tail law implies no finite root | **NO — TAIL ONLY** | it does cross, at `ω = 1.4617` |
| the UV falloff constant `−4.25` is a new constant | **NO — IT IS BORN, AND IT IS `−4a`** | the *local* slope descends `−4.72, −4.43, −4.27, −4.18, −4.13` toward the predicted `−4` |
| the throat factorises as two mouths plus a cavity | **NOT BY THIS GEOMETRY** | one smooth symmetric barrier; not imposed here |
| the static conductance `g` is `T₀(0)` | **NO — IT IS THE THRESHOLD COEFFICIENT** | `\|T₀\| → (π/8)(aω)²`, ratio `1.0003` at `ω = 0.01` |
| the traversable throat is free | **NO** | NEC violated everywhere; null-ray integral exactly `−3/(16G₅a)` |
| that NEC integral is the cost of the clock offset | **NO — DIFFERENT QUESTION** | static support vs aging history |
| one clock offset closes the loop | **YES, AT `ω = 1.4617`** — *relative to this basis* | `Ψ = 2πn` on the declared topology, stable to `1.0e-6` under `edge`/`steps`/`fd` refinement |
| `Λ = 1` can occur | **NO SUCH POINT FOUND ON THE TESTED BAND — NOT A THEOREM** | direct search for `R(ω)=0` on `[0.05,12]` finds no interior zero; nothing is claimed outside it. Deficit `~ e^{−4aω}` from Born |

---

## What this round establishes, and what it does not

PR #216 relocated the advanced wave from a bare postulate into a
Morris–Thorne–Yurtsever geometry. That was real progress. This round supplies
the geometry it was missing and **drives PR #216's own loop with it, through
the APIs that already existed**. Three findings: the relocation has an exact
price `−3/(16G₅a)`; on the topology the repository declares, one clock offset
*does* serve both the carrier and the packet, at `ω = 1.4617`; and no
perfect-transmission point was found on the tested band.

**The headline reversed during review, and the reason is worth recording.**
The first draft searched with `η_topo = +1` and concluded that simultaneous
closure was a UV limit never reached at finite frequency. That `η` was chosen,
not derived. `ConjugatePair` asserts the mouths of one throat carry *opposite*
orientations, which puts a factor `−1` in `η_topo` and shifts `Ψ` by `π` — and
a root appears. Two further errors were hiding behind the convenient sign: the
UV tail law had a missing `arg η_topo` subtraction that is invisible when the
constant is zero, and the `1/ω` decay was read as proving no finite root
existed, which an asymptotic law cannot do.

**Three claims are weaker than the first draft's.** The `Δ_phase`/`Δ_group` gap
was branch dependent and is withdrawn as a verdict. `|T| < 1` at finite `ω`
does **not** follow from `V > 0` — positive barriers can have
perfect-transmission resonances — so what stands is a band-limited search that
found none. And the closure verdict itself is **gauge dependent**: `Ψ` sweeps
less than `2π`, so a constant rephasing of the Jost basis can create or remove
the root. What is invariant is `dΨ/dω = −ωθ″` and the total variation. The
topological part of that constant is now derived; the basis part is not, and
fixing it needs the finite-mouth matching.

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

**Still open:** the history that produces `Δ_BA`; the finite-mouth junction to
the `S³` exterior, which is also what would fix the Jost basis constant and turn
the closure verdict from basis-relative into basis-free; whether the physical
probe is the scalar or a spinor, since the two channels answer differently; and
whether any BAM exchange kernel is meant to approximate the whole-throat `T_ℓ`
derived here.

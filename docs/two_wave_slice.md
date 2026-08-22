# Two waves, and where they connect inner to outer

**Module:** `geometrodynamics/viz/two_wave_slice.py`
**Probe:** `python -m experiments.closure_ledger.two_wave_slice_probe` (8/8)
**Tests:** `tests/test_viz_two_wave_slice.py` (36)
**Renderer:** `python scripts/geometrodynamics_v66_two_wave_slice.py --still v66.png`
**Superseded reading:** `docs/one_surface.md` (v67)

---

> **The construction here is ambiguous, and v67 says why.** This page draws
> **two** curves and reads their overlap through the seam as a connection. Two
> curves in one frame are two surfaces. If the two contributions are two pieces
> of *one* scalar deformation of *one* surface there is only ever one curve, and
> the question becomes whether it reaches both boundaries.
>
> **Every number below survives**, because the quantity plotted here as a
> "separation between two curves" *is* the one-surface deformation — the same
> array, to two ulps of `R_mid`. But the two configurations **swap names**, which
> inverts the headline: the cheapest-when-co-located result belongs to the
> *like* pair, and the identically-zero bisector is the **node of the opposed
> field**. Read `docs/one_surface.md` alongside this.

## What v46 left open

v46 put **one** scalar wave on the great circle through a source and its
antipode, drew it as a radial height `r = R_mid + εu(|σ|)` in the vacuole, and
glued `R_outer` to `R_inner` so the radial direction is a circle and the picture
lives on a torus. Its result was negative, and sharp:

> the curve is a **graph** `r = f(σ)`, so its radial winding number is
> identically zero. Every outward crossing of the seam is paid for by an inward
> one. **A height field cannot wind**, and one wave running to its antipode
> never meets itself.

That is a statement about *one* wave. Everything since has needed **two** — two
mouths, and an interference term that is bilinear and vanishes unless both are
present. So:

> one wave pulsing **outward** and one pulsing **inward**, both refocusing at
> the antipode — do they connect, inner to outer, and where?

## The construction

Two curves over the same circle, mirror images about the mid-radius:

    r_A(σ) = R_mid + ε u(d_A)      (driven outward)
    r_B(σ) = R_mid − ε u(d_B)      (driven inward)

with `d_A = |σ|`, and `d_B` measured from `B`'s own source — on top of `A`'s
(`offset = 0`) or at the far side (`offset = π`).

On the torus two points at the same `σ` are joined by two radial arcs. The one
that matters is the arc **through the seam**: out past `R_outer`, round, back in
at `R_inner`. Its length is `gap − |δ|` with `δ = r_A − r_B = ε(u_A + u_B)`, and
it closes when `|δ| = gap`.

## The answer

**Yes — at the antipode, on the seam, at the refocus.**

| | |
|--|--|
| where | `σ = ±π`, the antipode |
| at what radius | outward-driven → `0.740` = `R_inner`; inward-driven → `1.260` = `R_outer` |
| distance to the seam | `0.0` |
| contact gain there | `0.4997` |

Glued, `R_inner` and `R_outer` are one point. The two pulses meet *on the seam*.

### And the threshold is not a new number

A single wave crosses the seam when `R_mid + εu > R_outer`, i.e. `εu > gap/2`.
The pair spans `|δ| = 2ε|u|` of the radial circle and touches through the seam
when that reaches `gap` — i.e. `εu = gap/2`. **The same inequality.**

| | |
|--|--|
| single-wave wrap gain (at the run peak) | `0.220059` |
| pair-contact gain (same time, same `σ`) | `0.220059` |
| difference | `0.0` |

v46's *"the wave comes back inside the circle"* and this round's *"the two
pulses connect inner to outer"* are **one event described twice**. That is an
identity, not a coincidence, and it is asserted at zero rather than at a
tolerance.

### The roles are the other way round from the guess

The antipodal refocus is a **rarefaction** (`u < 0`). So it is the
*inward*-driven wave that bulges out to `R_outer`, and the *outward*-driven one
that dips to `R_inner`. Naming a wave "outward" says which way it is driven, not
which way it goes.

## What two waves can do that one cannot

| gain / threshold | covered fraction | connected | arcs |
|--|--|--|--|
| `0.980` | `0.9800` | no | 0 |
| `0.999` | `0.9990` | no | 0 |
| `1.000` | `1.0000` | **yes** | 1 |
| `1.050` | `1.0500` | yes | 1 |
| `1.300` | `1.3000` | yes | 1 |

Below threshold nothing connects. At threshold there is a single tangency at
`σ = ±π`. Above it, that point opens into **one arc**, bounded by two genuine
crossings, on which the band between the two curves covers the *entire* radial
circle — at those `σ` no radius is outside the pair.

A single wave past its own wrap threshold does nothing of the kind. It crosses
the seam, reappears inside, and still leaves every radius outside itself,
because it is a graph — re-checked here at four gains, still winding zero.

> **Two graphs bound a band, and a band can be radially surjective.** That is
> the whole difference between one wave and two, and it is the reason the
> question needed asking again.

## Where else they could have met, and why they do not

With the sources at opposite ends, the two travelling pulses cross at the
quarter points at `t = π/2`. That looks like the natural place for an
inner-to-outer connection. It is the **worst** one — they partially cancel, so
`|u_A + u_B|` is *smaller* there than either pulse alone.

| configuration | cheapest gain | at `t/π` | mid-flight penalty |
|--|--|--|--|
| co-located | `0.2226` | `1.99` | `7.4×` |
| antipodal sources | `0.4506` | `1.99` | `9.0×` |

(A first draft of this page guessed `4×` from a partial scan; the measured
penalties are those.)

The cheapest connection is always at a **refocus**, in both configurations. And
a co-located pair reaches it at `2.02×` less gain than an antipodal one: at a
refocus *both* of a co-located pair are at peak, while only one of an antipodal
pair is, so `|u_A + u_B|` is twice as large.

## What is still put in

Everything v46 put in is still put in, and it is worth restating because this
round's result is more positive than that one's and so easier to overread.

* The crossing rule that glues `R_outer` to `R_inner` is a **representation**
  choice, not a derived boundary condition.
* Nothing makes either wave dynamically aware of the seam, or of the other wave.
  The field is a **linear** scalar on a **fixed** round background, so the two
  waves do **not interact**. They are drawn on the same torus and the question
  is only whether their *images* meet.
* The gain is a **display** amplitude, as it was in v46. Every threshold quoted
  here is in those units and is a statement about the picture.

So this is not a claim that two physical waves reconnect a throat. It is the
statement that the obstruction v46 found — a single height field cannot wind, so
one wave can never meet itself — **does not apply to two**, and that the
amplitude at which it stops applying is exactly the one v46 already reported.

---

# Off the degenerate axis: the offset, and the signs

Everything above is the **co-located** pair, and that is the most degenerate
configuration the construction has. Both wave histories hang off one antipodal
axis, so bringing them together at one pole invites either exact overlap or
exact cancellation and tests neither. The question that needs the sources apart
is:

> can an inner-going branch launched from one axis meet an outer-going branch
> that has crossed the identification and re-entered at the inner boundary on a
> *different* axis — and does such a pair reach places a like-signed pair
> cannot?

## Two knobs

`offset` is the angle `α` between the two sources; `signs` is the radial sense
each wave is driven in.

    r_A(σ) = R_mid + s_A ε u(d_A)        d_A = |σ|
    r_B(σ) = R_mid + s_B ε u(d_B)        d_B = |σ − α|   (wrapped)

    δ = r_A − r_B = ε (s_A u_A − s_B u_B)

**Opposed** signs give the *sum* of the two fields; **like** signs give the
*difference*. That single line is the whole asymmetry.

## First, a correction to the framing

The expectation going in was that inner–inner and outer–outer would behave
differently. They do not — they are **the same case**, exactly:

| | |
|--|--|
| `(out,out)` vs `(in,in)`, as a difference of fields | `0.0` — the same bits |
| the same, through the drawn radii | `1.0` ulp of `R_mid` |
| `(out,in)` vs `(in,out)` | likewise |

Flipping both signs is a reflection about `R_mid`, which is an isometry of the
glued radial circle. So there are **two** configurations here, not four, and the
reason is worth being blunt about:

> the radial direction in this picture carries the field's **amplitude**, not its
> direction of propagation. A picture in which an inner-going branch differs
> from an outer-going one would have to encode propagation in the curve. This one
> does not.

That is a limitation of the representation, stated rather than worked around.
The two cases it *can* distinguish are **opposed** and **like-signed**, and that
distinction turns out to carry the answer.

(The 1-ulp row is not a hedge. A first draft asserted the drawn-radius residue
at zero and it failed; the identity was never in doubt, only which quantity it
was being asserted about. `(R_mid + εu_A) − (R_mid + εu_B)` and
`(R_mid − εu_A) − (R_mid − εu_B)` round differently.)

## The bisector

`σ = α/2` is equidistant from both sources, so `u_A = u_B` there **identically**
— at every time, at every amplitude. Then:

* **like-signed:** `δ ≡ 0`. The two curves are *the same curve* on that axis.
  They are never separated by so much as a hair, so no gain however large
  carries them through the seam there.
* **opposed:** `δ = 2ε u(α/2)` — as large as it can be.

The bisector is the axis where one pair is maximally connected and the other is
identically not. There are **two** of them, `α/2` and `α/2 − π`.

| offset `α/π` | bisector `σ/π` | like-signed `\|δ\|` | opposed threshold | far bisector | its threshold |
|--|--|--|--|--|--|
| `0.00` | `+0.0000` | `0.0e+00` | `0.2226` | `-1.0000` | `0.2786` |
| `0.15` | `+0.0750` | `0.0e+00` | `0.8605` | `-0.9250` | `0.7223` |
| `0.30` | `+0.1500` | `0.0e+00` | `1.1735` | `-0.8500` | `1.0623` |
| `0.50` | `+0.2500` | `0.0e+00` | `1.4367` | `-0.7500` | `1.3562` |
| `0.70` | `+0.3500` | `0.0e+00` | `1.5925` | `-0.6500` | `1.5430` |
| `0.85` | `+0.4250` | `6.6e-16` | `1.6506` | `-0.5750` | `1.6255` |
| `1.00` | `+0.5000` | `0.0e+00` | `1.6609` | `-0.5000` | `1.6609` |

The far bisector is the **cheaper** of the two at every genuine offset, because
it sits nearer the antipodal caustic where the field is amplified. The residues
in the third column are the floating-point floor of the wrapped distance — the
`mod`, not the field; worst relative residue over the whole sweep is `4.5e-15`.

## The answer

**Yes.** Drive each configuration to `1.15×` the opposed pair's own bisector
threshold, at the time that threshold is reached:

| offset `α/π` | bisector `σ/π` | arc (rad) | centre − bisector | like-signed samples on it | to nearest source | to nearest antipode |
|--|--|--|--|--|--|--|
| `0.15` | `+0.0750` | `0.0960` | `0.0e+00` | **0** | `0.075π` | `0.925π` |
| `0.30` | `+0.1500` | `0.0960` | `0.0e+00` | **0** | `0.150π` | `0.850π` |
| `0.50` | `+0.2500` | `0.0960` | `0.0e+00` | **0** | `0.250π` | `0.750π` |
| `0.70` | `+0.3500` | `0.0960` | `0.0e+00` | **0** | `0.350π` | `0.650π` |
| `1.00` | `+0.5000` | `0.0960` | `0.0e+00` | **0** | `0.500π` | `0.500π` |

> An opposed (inner–outer) pair connects through the seam on an arc **centred on
> the bisector to machine zero**, off both the sources and their antipodes, on
> which a like-signed pair connects at **no amplitude at all**.

At `α = 0` the bisector collapses onto the source axis and there is nothing
off-axis to find. That is the degeneracy — recovered by measurement, as a
coordinate fact, rather than assumed.

## What the slider does

| | |
|--|--|
| threshold at `α = 0` | `0.2201` |
| threshold at `α = π` | `1.6639` |
| over the sweep | `7.56×` |
| timing against `t = α/2` | `0.0031π` |
| price of exclusivity | `1.68–3.74×` |
| cheapest point pins to an axis from | `α = 0.125π` |

Turning `α` up slides the exclusive connection continuously from the source axis
to the quarter point. The **timing** is the pulse-crossing time: the two
outgoing pulses reach the bisector together at `t = α/2`, and the measured peak
leads that by a constant `0.107` rad — which is where a single pulse's own
extremum sits relative to its geodesic arrival, checked against one wave rather
than fitted.

And **exclusive is not cheap.** From `α = 0.125π` up, the globally cheapest
connection sits exactly on one of the four axes — either source or either
antipode — is available to *both* pairs, and costs `1.7–3.7×` less. Below that
offset it drifts off axis by at most `0.03π`, while the two return focuses are
still close enough to blend.

### One number with a history

The threshold rises monotonically in `α` but for a turn-over of about `0.02%`
into the symmetric endpoint `α = π`, visible only on a fine enough sweep.

A first pass put that dip at `0.08%` — four times larger — and "confirmed" it by
refining the **time** grid fourfold and watching it fail to shrink. It failed to
shrink because time was never the problem: that pass evaluated the bisector at
the nearest point of the `σ` grid, and at `α = 0.958π` the bisector falls exactly
halfway between two samples.

> Refining the axis a discrepancy does not live on will always leave it standing,
> and leaving it standing is not evidence.

The bisector is evaluated off-grid now, at the angle it actually has. What
survives that is one fifth the size and still real.

## Scope, again

Unchanged, and it still binds. The crossing rule is a representation choice; the
field is a linear scalar on a fixed round background, so the two waves do not
interact; the gain is a display amplitude. Nothing here says two physical waves
reconnect a throat on a bisector. It says that in this construction the opposed
pair has an off-antipodal contact locus that the like-signed pair provably lacks,
that the locus is pinned by the source geometry rather than by the amplitude, and
that the offset is what brings it into existence.

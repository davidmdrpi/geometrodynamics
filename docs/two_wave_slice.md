# Two waves, and where they connect inner to outer

**Module:** `geometrodynamics/viz/two_wave_slice.py`
**Probe:** `python -m experiments.closure_ledger.two_wave_slice_probe` (4/4)
**Tests:** `tests/test_viz_two_wave_slice.py` (19)
**Renderer:** `python scripts/geometrodynamics_v66_two_wave_slice.py --still v66.png`

---

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

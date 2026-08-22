# Two-wave slice probe — do the antipodal pulses connect?

**Question.** one wave pulsing outward and one pulsing inward, both refocusing at the antipode -- do they connect, inner to outer, and where?

**Answer.** yes: at the antipode, on the seam, at the refocus -- and at exactly the amplitude where one wave would have wrapped

| | |
|--|--|
| where | `σ = -1π` — the antipode |
| at what radius | outward-driven → `0.740` = `R_inner`, inward-driven → `1.260` = `R_outer` |
| distance to the seam | `0.0e+00` |
| contact gain there | `0.499716` |
| single-wave wrap gain | `0.220059` (at the run peak, `t = 6.236`) |
| the two thresholds differ by | `0.0e+00` |

**8/8 checks pass.**

| id | check | result |
|----|-------|--------|
| T1 | *** one wave wrapping IS two touching *** | PASS |
| T2 | *** where: antipode, seam, refocus *** | PASS |
| T3 | *** the band covers what one graph cannot *** | PASS |
| T4 | meeting mid-flight is harder | PASS |
| T5 | two cases, not four (a limitation) | PASS |
| T6 | *** the bisector is degenerate for like signs *** | PASS |
| T7 | *** an off-antipodal connection only one pair has *** | PASS |
| T8 | the slider: it moves, and it costs | PASS |

## Below, at, and above threshold

| gain / threshold | covered fraction | connected | arcs |
|--|--|--|--|
| `0.980` | `0.9800` | no | 0 |
| `0.999` | `0.9990` | no | 0 |
| `1.000` | `1.0000` | yes | 1 |
| `1.001` | `1.0010` | yes | 1 |
| `1.050` | `1.0500` | yes | 1 |
| `1.300` | `1.3000` | yes | 1 |

> a single wave is a graph, so its radial winding is zero and every radius stays outside it however hard it is driven; two graphs bound a band, and past threshold that band leaves no radius out

## Where else they could have met

| configuration | cheapest gain | at `t/π` | mid-flight penalty |
|--|--|--|--|
| co located | `0.2226` | `1.99` | `7.4×` |
| antipodal sources | `0.4506` | `1.99` | `9.0×` |

The cheapest connection is always at a **refocus**. A co-located pair reaches it at `2.02×` less gain than an antipodal one.

## Off the degenerate axis

**Question.** does an inner-outer pair possess off-antipodal intersections that neither inner-inner nor outer-outer pairs possess?

**Answer.** yes -- an opposed (inner-outer) pair connects through the seam on an arc centred on the bisector between the two source axes, which is off both the sources and their antipodes, and on which a like-signed pair is identically coincident and so connects at no amplitude at all

### There are two configurations, not four

| | |
|--|--|
| `(out,out)` vs `(in,in)`, as fields | `0.0e+00` |
| the same, through the drawn radii | `1.0` ulp of `R_mid` |

> the radial direction carries field amplitude, not direction of propagation, so this representation cannot distinguish inner-inner from outer-outer

### The bisector

| offset `α/π` | bisector `σ/π` | like-signed `|δ|` | opposed threshold | far bisector | its threshold |
|--|--|--|--|--|--|
| `0.00` | `+0.0000` | `0.0e+00` | `0.2226` | `-1.0000` | `0.2786` |
| `0.15` | `+0.0750` | `0.0e+00` | `0.8605` | `-0.9250` | `0.7223` |
| `0.30` | `+0.1500` | `0.0e+00` | `1.1735` | `-0.8500` | `1.0623` |
| `0.50` | `+0.2500` | `0.0e+00` | `1.4367` | `-0.7500` | `1.3562` |
| `0.70` | `+0.3500` | `0.0e+00` | `1.5925` | `-0.6500` | `1.5430` |
| `0.85` | `+0.4250` | `6.6e-16` | `1.6506` | `-0.5750` | `1.6255` |
| `1.00` | `+0.5000` | `0.0e+00` | `1.6609` | `-0.5000` | `1.6609` |

> sigma = alpha/2 is equidistant from both sources, so u_A = u_B there identically; like signs subtract that to zero and opposed signs add it to twice one field

### The exclusive arc

Driven to `1.15×` the opposed pair's own bisector threshold:

| offset `α/π` | bisector `σ/π` | arc (rad) | centre − bisector | like-signed samples on it | to nearest source | to nearest antipode |
|--|--|--|--|--|--|--|
| `0.15` | `+0.0750` | `0.0960` | `0.0e+00` | **0** | `0.075π` | `0.925π` |
| `0.30` | `+0.1500` | `0.0960` | `0.0e+00` | **0** | `0.150π` | `0.850π` |
| `0.50` | `+0.2500` | `0.0960` | `0.0e+00` | **0** | `0.250π` | `0.750π` |
| `0.70` | `+0.3500` | `0.0960` | `0.0e+00` | **0** | `0.350π` | `0.650π` |
| `1.00` | `+0.5000` | `0.0960` | `0.0e+00` | **0** | `0.500π` | `0.500π` |

### What the slider does

| | |
|--|--|
| threshold at `α = 0` | `0.2201` |
| threshold at `α = π` | `1.6639` |
| over the sweep | `7.56×` |
| timing vs `t = α/2` | `0.0031π` |
| price of exclusivity | `1.68–3.74×` |
| cheapest point pins to an axis from | `α = 0.125π` |

It rises monotonically but for a turn-over of `0.017%` into the symmetric endpoint `α = 1.00π`. A first pass put that at `0.08%` and confirmed it by refining the *time* grid — the wrong axis. The bisector is evaluated off-grid here.

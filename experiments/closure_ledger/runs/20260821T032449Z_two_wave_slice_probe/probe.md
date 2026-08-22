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

**4/4 checks pass.**

| id | check | result |
|----|-------|--------|
| T1 | *** one wave wrapping IS two touching *** | PASS |
| T2 | *** where: antipode, seam, refocus *** | PASS |
| T3 | *** the band covers what one graph cannot *** | PASS |
| T4 | meeting mid-flight is harder | PASS |

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

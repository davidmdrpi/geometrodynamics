# Physical-throat probe — which throat, and the sign it gives

**Question.** the throat's area and length were free parameters. which values are physical?

**Answer.** the vacuum throat, glued with no surface layer -- an Einstein-Rosen bridge with no free parameter at all: f0 = sin^3 a is forced, and f0 = 2M, so M = sin^3(a)/2

| | this round (vacuum throat) | PR #264 (product tube) |
|--|--|--|
| neck / area | `f₀ = 1.24844e-04` | `𝒜 = 12.5664` |
| length | `0.10075` (forced) | `0.9` (chosen) |
| shunt | `0.0` — zero by identity | `6.0702` |
| mass | `M = 6.24219e-05` — **derived** | not defined |
| `ΔA/A` | `(+6.6375e+00, +8.5813e+00)` | `(-2.0575e-03, -1.8818e-03)` |
| verdict | **opens / opens** | closes / closes |

In units of 2 pi G, with the source being PR #263's interference dT00.  Mass law: `M = sin^3(a) / 2, in units of the ESU curvature radius`.

**7/7 checks pass.**

| id | check | result |
|----|-------|--------|
| T1 | the curvature, derived from the metric | PASS |
| T2 | the gluing forces the neck radius, L and I | PASS |
| T3 | what the other throats would need | PASS |
| T4 | *** there is no cavity *** | PASS |
| T5 | *** the shunt decides the sign *** | PASS |
| T6 | *** it is an Einstein-Rosen bridge, and M is derived *** | PASS |
| T7 | *** the answer, and it reverses *** | PASS |

## The mechanism

| conductance | shunt | `ΔA/A` mouth 1 | verdict |
|--|--|--|--|
| `3.9270e-04` | `0.0000` | `+6.63750e+00` | opens/opens |
| `1.6042e+01` | `0.0000` | `+7.60933e+00` | opens/opens |
| `3.9270e-04` | `6.0702` | `-1.59264e-03` | closes/closes |
| `1.6042e+01` | `6.0702` | `-2.05839e-03` | closes/closes |

The conductance is scanned over eight orders and never changes the sign. The shunt passes through a pole near `2e-03` and flips it. A vacuum tube has zero shunt **by identity** — `(f²u')' = 0` — so there is nowhere for monopole flux to go.

**What it costs.** the response is 3000x larger and scales as a^-3, because a vacuum throat has zero shunt by identity and so barely lifts the constraint's degeneracy.

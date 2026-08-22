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

**9/9 checks pass.**

| id | check | result |
|----|-------|--------|
| T1 | the curvature, derived from the metric | PASS |
| T2 | the gluing forces the neck radius, L and I | PASS |
| T3 | what the other throats would need | PASS |
| T4 | *** there is no cavity *** | PASS |
| T5 | *** the shunt decides the sign *** | PASS |
| T6 | *** it is an Einstein-Rosen bridge, and M is derived *** | PASS |
| T7 | *** the answer, and it reverses *** | PASS |
| T8 | *** where the Riccati solve stops resolving (PR #267) *** | PASS |
| T9 | *** the mouth-to-mouth hierarchy C_l ~ a^(4l+3) *** | PASS |

## The mechanism

| conductance | shunt | `ΔA/A` mouth 1 | verdict |
|--|--|--|--|
| `3.9270e-04` | `0.0000` | `+6.63750e+00` | opens/opens |
| `1.6042e+01` | `0.0000` | `+7.60933e+00` | opens/opens |
| `3.9270e-04` | `6.0702` | `-1.59264e-03` | closes/closes |
| `1.6042e+01` | `6.0702` | `-2.05839e-03` | closes/closes |

The conductance is scanned over eight orders and never changes the sign. The shunt passes through a pole near `2e-03` and flips it. A vacuum tube has zero shunt **by identity** — `(f²u')' = 0` — so there is nowhere for monopole flux to go.

**What it costs.** the response is 3000x larger and scales as a^-3, because a vacuum throat has zero shunt by identity and so barely lifts the constraint's degeneracy.

## The two-port in closed form, and where the solve stopped resolving

`f = f₀cosh²x` with `ds = 2f dx` turns `(f²u')' = ℓ(ℓ+1)u` into `y'' = (2ℓ+1)²y`, and `e^{−X} = tan(a/2)` exactly.  With `k = 2ℓ+1` and `q = tan^{2k}(a/2)`:

> `D_ℓ = −2π sin a [ k(1+q²)/(1−q²) − cos a ]`  ,  `C_ℓ = +4π sin a · kq/(1−q²)`

That is now the production admittance at `a = 0.05`.  The Riccati solve is kept as an independent validator, not deleted.

| `ℓ` | `C_ℓ` closed form | `C_ℓ` Riccati | rel. error | signs | diagonal rel. error |
|--|--|--|--|--|--|
| 0 | `+3.92699e-04` | `+3.92699e-04` | `7.48e-14` | agree | `7.5e-14` |
| 1 | `+4.60578e-10` | `+4.60517e-10` | `1.32e-04` | agree | `9.7e-14` |
| 2 | `+3.00105e-16` | `-1.16573e-14` | `3.98e+01` | **opposite** | `9.4e-15` |
| 3 | `+1.64257e-22` | `-7.54952e-15` | `4.60e+07` | **opposite** | `4.0e-15` |

The cross term is a *difference* of two eigenchannel values in the solve and a *product* of small factors in the closed form.  The floor is around `1e-12`; the diagonal was never affected.

The `39×` factor is deliberately **not** pinned: it is one solver's step sequence in one build of SciPy and would move under any of them.  What is pinned is the boundary — that the sign is wrong and the magnitude is more than an order of magnitude out, in a channel whose honest size is below the solver's floor.

## The hierarchy

`C_l ~ 4 pi (2l+1) sin(a) tan^(4l+2)(a/2) ~ a^(4l+3)` — fitted exponents `3.000000`, `7.000004`, `11.000009`, `15.000013` against `3, 7, 11, 15`.

| `a` | `C₀` | `C₁` | `C₀/C₁` | `C₁/D₁` |
|--|--|--|--|--|
| `0.05` | `3.926992e-04` | `4.605780e-10` | `8.5262e+05` | `7.329e-10` |
| `0.10` | `3.141614e-03` | `5.910168e-08` | `5.3156e+04` | `4.699e-08` |
| `0.20` | `2.513546e-02` | `7.641316e-06` | `3.2894e+03` | `3.031e-06` |
| `0.30` | `8.487017e-02` | `1.327738e-04` | `6.3921e+02` | `3.497e-05` |

the four n=1 harmonics x^A split locally as 1 (X^0 = cos chi) + 3 (X^i = sin chi n^i), so the ESU kernel splits `1 ⊕ 3` and the two pieces cross four powers of `a` apart.

**The statement this supports:** the static scalar Laplacian on this scalar-flat spatial throat suppresses the local l=1 mouth-to-mouth channel by ~1e-09 at a = 0.05, while preserving a much stronger monopole channel.

**What is not claimed:** that orientation information cannot cross the throat -- one operator, one slice, no lapse chosen, and the l=1 channel is small rather than zero.

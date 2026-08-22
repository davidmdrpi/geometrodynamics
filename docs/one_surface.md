# One field on one surface

**Module:** `geometrodynamics/viz/one_surface.py`
**Probe:** `python -m experiments.closure_ledger.one_surface_probe` (9/9)
**Tests:** `tests/test_viz_one_surface.py` (35)
**Renderers:** `geometrodynamics_v67_one_surface.py`, `geometrodynamics_v68_two_fronts.py`
**Supersedes the reading of:** `docs/two_wave_slice.md` (v66)

---

## The objection, and it lands

v66 drew **two** curves — `r_A = R_mid + s_A ε u_A` and `r_B = R_mid + s_B ε u_B`
— and asked whether their images meet through the glued seam. Its docstring was
careful that the two waves do not interact, but that caveat does not repair the
construction; it only labels the problem. **Two curves in one frame are two
surfaces**, and reading their overlap as a connection is a statement about a
picture, not about a field.

If the two contributions are two pieces of one scalar deformation of one
surface, there is only ever one curve:

    u(θ, t) = s_A u_A + s_B u_B          ← one field
    r(θ, t) = R_mid + ε u(θ, t)          ← one curve

and the question is not whether two images meet, but whether that single curve
reaches `R_outer` at some `θ` and `R_inner` at another — so that the surface
itself passes through the identification.

## The repair costs nothing, and that is the surprise

`δ = r_A − r_B = ε(s_A u_A − s_B u_B)` — the quantity v66 plotted as a
*separation between two curves* — **is** the one-surface deformation with the
second sign flipped. The same array:

| | |
|--|--|
| v66 separation vs one-surface field, worst over 5 offsets × a full run | `2.0` ulp of `R_mid` |

Two ulps and not zero because v66's `separation` forms
`(R_mid + a) − (R_mid + b)`, so it carries the mid-radius's rounding whatever the
fields do. The identity is about the fields. (Second time this exact species has
appeared in this arc — see §25.)

So every number v66 reported survives the change of construction, with the two
configurations **swapping names**:

    v66 "like-signed"  δ = ε(u_A − u_B)   ==   one surface, OPPOSED
    v66 "opposed"      δ = ε(u_A + u_B)   ==   one surface, LIKE-signed

That inversion is not cosmetic. v66's headline — *the opposed pair connects most
cheaply when co-located* — is, correctly read, a statement about the **like**
pair. And v66's *the like pair is identically zero on the bisector* is the
**node of the opposed field**, which is where it belongs and where it makes
sense: `u_A = u_B` on the equidistant axis, so the opposed combination vanishes
there identically.

## Coincident foci cancel, exactly

| | |
|--|--|
| the opposed field at `α = 0`, worst over a run | `0.0` — exactly zero |
| its span threshold | `∞` |
| the like pair's, unaffected | `0.4908` |

Two identical contributions with opposite orientation subtract to nothing. No
amplitude connects them, and the divergence of the required gain is not an
artefact — it is the statement that there is no deformation at all.

## The monochromatic law

For one spatial mode `m`, equal amplitude, opposite orientation, foci at `∓α/2`:

    u = A cos[m(θ+α/2) − ωt] − A cos[m(θ−α/2) − ωt]
      = −2A sin(mα/2) sin(mθ − ωt)

verified symbolically and on a grid. So the deformation amplitude is
`B = 2A|sin(mα/2)|` and the peak-to-trough span is `Δr = 4A|sin(mα/2)|`.

**The optimum is half a wavelength**, `α* = π/m` — and `(2j+1)π/m` after it. The
antipode is simply the `m = 1` member of that family.

| mode `m` | `α*` = half a wavelength | at the antipode | verdict |
|--|--|--|--|
| `1` | `1.0000π` | `2.0000` | **maximal** |
| `2` | `0.5000π` | `0.0000` | **cancels** |
| `3` | `0.3333π` | `2.0000` | **maximal** |
| `4` | `0.2500π` | `0.0000` | **cancels** |
| `5` | `0.2000π` | `2.0000` | **maximal** |
| `6` | `0.1667π` | `0.0000` | **cancels** |
| `7` | `0.1429π` | `2.0000` | **maximal** |
| `8` | `0.1250π` | `0.0000` | **cancels** |

> **The antipode is parity-dependent.** `sin(mπ/2)` is `±1` for odd `m` and `0`
> for even `m`. "Opposite poles give the greatest separation" is a statement
> about the *lowest* mode, not a general one — for even modes opposite poles
> cancel exactly.

With `D = R_outer − R_inner = 0.52`, the curve spans the gap when
`4A|sin(mα/2)| ≥ D`, i.e. `A_req = D/(4|sin(mα/2)|)` — properly singular as
`α → 0`. Hand-checked against a bisection on the measured amplitude:

| `α` | `A_req` |
|--|--|
| `180°` | `0.130000` |
| `90°` | `0.183848` |
| `60°` | `0.260000` |
| `30°` | `0.502281` |

## Where the two models part company

The parity above is not a plane-wave artefact — it is in the real spectrum. The
ESU eigenfunctions here are `Z_n(χ) = sin[(n+1)χ]/[(n+1) sin χ]` with
`ω_n = n + 1`, the same `ω_n = n+1` this repository derives for the conformally
coupled ESU. `Z_n(π) = (−1)ⁿ`, checked, so the antipodal difference field
doubles for odd `n` and cancels for even `n` — **exactly**, ratio `2.0000` and
`0.0000`.

But the *location of the optimum* does not carry across, and this is the one
place the derivation and the measurement disagree.

A zonal harmonic is **centred**: `Z_n(0) = 1` is a global maximum and
`|Z_n| ≤ 1` everywhere. So the opposed pair obeys `|Z_A − Z_B| ≤ 2`, with
equality only where one focus sees `+1` and the other `−1` **at the same
point** — which happens exactly when the two poles are antipodal and `n` is odd.

| zonal order `n` | `ω = n+1` | measured `α*` | peak strength | at the antipode | at half a wavelength |
|--|--|--|--|--|--|
| `1` | `2` | `1.0000π` | `2.0000` | `2.00e+00` | `1.4142` |
| `2` | `3` | `0.5002π` | `1.3333` | `1.51e-13` | `1.1547` |
| `3` | `4` | `1.0000π` | `2.0000` | `2.00e+00` | `1.1217` |
| `4` | `5` | `0.2902π` † | `1.2500` | `9.05e-14` | `1.1097` |
| `5` | `6` | `1.0000π` | `2.0000` | `2.00e+00` | `1.1039` |
| `6` | `7` | `0.2053π` † | `1.2330` | `6.46e-13` | `1.1006` |
| `8` | `9` | `0.8401π` † | `1.2266` | `5.03e-13` | `1.0971` |
| `10` | `11` | `0.1303π` † | `1.2234` | `4.11e-13` | `1.0954` |

> For zonal modes `α* = π` for **every** odd `n`, and it **saturates** the
> bound. Half a wavelength reaches only a fraction of it — `1.41` at `n = 1`,
> falling to `1.10` by `n = 5`. For even `n` the antipode cancels exactly and
> nothing reaches the bound at all: the best available is `1.22–1.33`.

† **The location is not unique for even `n ≥ 4`.** Two separations reach within
`1e-3` of the same peak, so which one `argmax` returns moves with the grid — a
1501-point sweep puts `n = 6` at `0.794π` and a 601-point sweep at `0.205π`,
agreeing on the height to seven digits. The **strength** is the result; the
location is reported with the degeneracy counted rather than hidden behind
whichever grid ran last. For odd `n` the maximiser is unique and is `π`.

A plane wave has no distinguished centre, so nothing picks out a separation
other than the wavelength, and `α* = π/m` follows. A zonal mode has one, and the
antipode is where two centres coincide with opposite sign. **The parity survives
the change of model exactly; the location of the optimum does not.**

The kernel this programme actually cares about is `n = 1`, which is odd — so for
it the antipode is both optimal *and* saturating.

## A pulse is not a mode

v46's field is a launched pulse. Measured against the `S²` zonal basis at the
refocus it carries a power-weighted mean of `n ≈ 10`, with **fifteen** modes
holding 90% of the power. So the single-mode laws apply to it mode by mode, not
in bulk — and one of them visibly fails.

Pulse width `0.18` rad. The one-surface threshold plateaus at `0.2163`
(spread `4.2e-04`):

| offset `α/π` | in pulse widths | monochromatic `A_req` | pulse threshold | ratio |
|--|--|--|--|--|
| `0.02` | `0.35` | `4.1387` | `0.6833` | `6.06×` |
| `0.04` | `0.70` | `2.0704` | `0.3655` | `5.66×` |
| `0.06` | `1.05` | `1.3814` | `0.2725` | `5.07×` |
| `0.10` | `1.75` | `0.8310` | `0.2217` | `3.75×` |
| `0.20` | `3.49` | `0.4207` | `0.2156` | `1.95×` |
| `0.40` | `6.98` | `0.2212` | `0.2161` | `1.02×` |
| `0.60` | `10.47` | `0.1607` | `0.2164` | `0.74×` |
| `1.00` | `17.45` | `0.1300` | `0.2166` | `0.60×` |

> Two localized pulses cancel only **while they overlap**. Past roughly one
> pulse width nothing cancels, and the threshold **saturates** at the
> single-pulse value instead of continuing down like `1/|sin|`.

The coincident-foci cancellation is real and exact for a pulse too. The law
governing its *approach* is not.

## The chord, and what it costs

At the optimum the outward and inward extrema sit `π/m` apart on the circle, so
the straight chord between them is

    L = √(D² + 4 R_out R_in sin²(π/2m))

— identical to the plain law of cosines to `2e-16`, and grouped so the purely
radial gap `D` is visible as the `Δθ → 0` limit.

| mode `m` | separation | span / `A` | bulk chord | `L/D` | span at fixed energy | spans the gap? |
|--|--|--|--|--|--|--|
| `1` | `1.0000π` | `2.0000` | `2.0000` | `3.846` | `4.0000` | yes |
| `2` | `0.5000π` | `2.0000` | `1.4612` | `2.810` | `2.0000` | yes |
| `3` | `0.3333π` | `2.0000` | `1.0967` | `2.109` | `1.3333` | yes |
| `4` | `0.2500π` | `2.0000` | `0.9037` | `1.738` | `1.0000` | yes |
| `6` | `0.1667π` | `2.0000` | `0.7213` | `1.387` | `0.6667` | yes |
| `8` | `0.1250π` | `2.0000` | `0.6421` | `1.235` | `0.5000` | **no** |
| `16` | `0.0625π` | `2.0000` | `0.5534` | `1.064` | `0.2500` | **no** |

At fixed **display amplitude** the span is flat at `2.0000` across the whole
half-wavelength family while the chord falls from `2.000` to `0.553` — the same
deformation on a progressively shorter connection, tending to the purely radial
gap `0.520`.

## But frequency is not free

`E ∝ ω²A²`, so at fixed **energy** `A ∝ 1/ω` and the attainable span falls as
fast as the chord does:

> low frequency buys a large deformation on a long connection; high frequency
> buys a short connection with a small deformation.

At fixed energy the highest mode that still spans the gap at all is `m = 7`.
That is a real bound and it is worth having, but it is *not* an optimum:

**Not claimed.** No favourable frequency, because that needs an energy
normalisation and a packet focusing law this model does not contain. The point
that matters for the visualisation is narrower and sharper: **a frequency slider
cannot hold displacement fixed and then be read as constant-energy physics.**

## Where each front sits on the one surface

The offset question, asked *inside* the one-surface object. As the two axes move
apart, where do `A` and `B` individually sit on that same surface, what are
their signs, where do they interfere, and how does the surface answer?

**Renderer:** `python scripts/geometrodynamics_v68_two_fronts.py --still v68.png`

### The rule the figure obeys

Only one thing is ever drawn as a closed curve in the annulus: **the surface**.
`c_A = s_A u_A` and `c_B = s_B u_B` are *components* of its deformation, so they
appear only on graphs of field against `σ`, where nothing invites reading them
as separate objects. Drawing them as closed curves is exactly the v66 error and
it is not repeated. What the annulus panels do instead is **colour the single
surface** by which front is live along each arc — and because an inward dent is
a negative contribution and an outward one positive, each front's sign is read
straight off the surface.

### What the sweep shows

| offset `α/π` | at `t/π` | peak `c_A` | peak `u` | amplification | overlap arc | `σ` of `c_A` | `σ` of `c_B` | `σ` of `u` |
|--|--|--|--|--|--|--|--|--|
| `0.00` | `0.01` | `0.6710` | `0.0000` | `0.0000` | `0.794` | `−0.01π` | `−0.01π` | `−1.00π` |
| `0.15` | `1.99` | `1.1796` | `1.1994` | `1.0168` | `0.148` | `+0.00π` | `+0.15π` | `+0.00π` |
| `0.25` | `1.99` | `1.1796` | `1.1978` | `1.0155` | `0.000` | `+0.00π` | `+0.25π` | `+0.00π` |
| `0.50` | `1.99` | `1.1796` | `1.1950` | `1.0131` | `0.000` | `+0.00π` | `+0.50π` | `+0.00π` |
| `0.75` | `1.99` | `1.1796` | `1.1938` | `1.0121` | `0.000` | `+0.00π` | `+0.75π` | `+0.00π` |
| `1.00` | `1.99` | `1.1796` | `1.1935` | `1.0118` | `0.000` | `+0.00π` | `−1.00π` | `−1.00π` |

Two things fall out, and the second was not what I expected to be writing.

**At `α = 0` the surface is a perfect circle.** Not a small deformation — a
circle. The two contributions are the same function of the same distance with
opposite sign, so they cancel identically, and there is nothing for the bulk to
connect at any amplitude.

**Past one pulse width the two contributions stop overlapping at all.** The
overlap arc goes to `0.000` by `α = 0.25π`, and the amplification
`max|u| / max|c_A|` sits at `1.012–1.017` across the whole rest of the sweep.
The total is one contribution plus almost nothing, and it peaks exactly where a
single contribution peaks.

> So the offset does **not** turn interference on. It turns the *cancellation*
> off, by pulling the two fronts apart — and what is left is two nearly
> independent dents in one surface, not a reinforced one.

That is a sharper and less flattering statement than "the offset produces an
off-antipodal connection", and it is what the decomposition says. Interference
in this construction is confined to the narrow arc where both fronts are
actually present, and for a localized pulse that arc closes as soon as the foci
clear each other. The monochromatic case is the opposite extreme — a mode fills
the whole circle, so it cancels everywhere and the `1/|sin(mα/2)|` law holds all
the way out.

### And in time

The right-hand panel of v68 is `σ` against `t` at `α = π/2`: each source
launches **two** arms, they cross at the quarter points, cross again, and the
surface answers wherever a front is — most strongly at a refocus. The crossings
are visible as the bright cells; they are where the overlap arc is briefly
nonzero even at large offset.

## Scope

Unchanged from v46 and v66, and it still binds. The crossing rule that glues
`R_outer` to `R_inner` is a **representation** choice, not a derived boundary
condition. The field is a **linear** scalar on a **fixed** round background. The
gain is a **display** amplitude — and the fixed-energy section is the one place
that distinction is made to do work.

What the round changes is not the physics but the object: v66 asked whether two
drawn curves meet, which was never a well-posed question about a field. This
asks whether one surface reaches both boundaries, which is.

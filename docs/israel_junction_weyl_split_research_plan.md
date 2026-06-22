# Israel junction audit for the non-orientable throat gluing — the braneworld Weyl split (PR #167)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The question, sharpened

The naive Israel audit is already settled: a thin-shell throat needs
NEC/WEC-violating **exotic** surface matter (`σ < 0`), and the
non-orientable (antipodal Z₂ / C-swap) self-gluing does **not** rescue the
sign. That is not the interesting probe. The honest, **not-predetermined**
question is the **braneworld split**:

> Does BAM's specific 5D Tangherlini bulk supply the throat's negative
> effective-4D `σ` through the *projected bulk Weyl term* `E_μν` — with
> **ordinary 5D matter** — or does it still need a genuine exotic source?

The deliverable is the split: how much of the required exotic stress comes
from the bulk Weyl projection versus how much remains irreducible. If the
bulk covers it, the "geometrically enforced, no exotic matter" claim
survives and the citation is **Bronnikov–Kim** (brane wormholes) /
**Dadhich–Maartens–Papadopoulos–Rezania** (tidal-charge brane black holes),
via the **Shiromizu–Maeda–Sasaki** effective brane equations.

## The eight deliverables

| # | deliverable | result |
|---|---|---|
| 1 | thin-layer stress tensor `S_ab` | `S^a_b = diag(−σ, p_t, p_t)`, Lanczos |
| 2 | surface energy density `σ` | `σ = −√f(a)/(2πa) < 0` (all `a>r_s`) |
| 3 | tangential pressure `p_t` | `p_t = (1/4π)[f'/(2√f) + √f/a]` |
| 4 | `S_ab k^a k^b` (null) | `∝ σ + p_t`; bulk NEC violated tangentially |
| 5 | NEC / WEC status | WEC violated on the shell (`σ<0`) — exotic |
| 6 | right sign of `σ`? | `σ<0`: right for a wormhole, wrong for ordinary matter |
| 7 | right scale of `σ`? | `\|σ\| ≈ 1/(2πa)`, the inverse throat scale |
| 8 | discrete `P` as thickness → 0 | a tanh wall of width `δ` → the discrete Lanczos `σ` |

## The decisive fact (computed)

The throat metric `f(r) = 1 − (r_s/r)²` (Tangherlini / tidal-charge form)
is **Ricci-flat**, `R = 0` (to machine precision). Its effective 4D stress
is **traceless** with the `r⁻⁴` "radiation" form

```
ρ_eff = −r_s²/(8πG r⁴) < 0,   p_r = −ρ_eff,   p_t = +ρ_eff,
−ρ + p_r + 2p_t = 0,           r⁴·ρ_eff = −r_s² (constant).
```

A traceless, `r⁻⁴` stress is exactly the form a **projected bulk Weyl
tensor** `E_μν` (traceless by construction) takes.

## The split

By Shiromizu–Maeda–Sasaki, a **vacuum brane** (`T_brane = 0`) obeys
`G_μν = −E_μν`, which forces `R = 0` — **satisfied here**. So the entire
effective-exotic stress is the projected bulk Weyl term, carrying tidal
charge `Q = −r_s²` (Dadhich et al.), sourced by the **ordinary 5D
Tangherlini vacuum** (Ricci-flat in 5D). The split is

```
~100% bulk Weyl projection   /   ~0% irreducible brane exotic matter.
```

The exotic-looking 4D stress is the **geometric shadow of the ordinary 5D
bulk**, not real exotic matter on the brane (Bronnikov–Kim). The surgical
thin-shell surface term itself **vanishes** as the gluing approaches the
throat (`f → 0`).

## Honest caveats

- `f = 0` at `r = r_s` is a **horizon** (null surface): the throat sits at
  a degenerate locus, and the vanishing surface term is a vanishing at a
  null surface, not a generic traversable neck.
- `R = 0` + traceless is the **necessary** 4D signature of a vacuum brane;
  the full 5D embedding (a 5D Tangherlini bulk whose Weyl projection
  reproduces this induced metric and `E_μν`) is the Dadhich/Bronnikov–Kim
  construction, **cited**, not re-solved as a 5D boundary-value problem.
- The tidal charge is **negative**; the 4D WEC is genuinely violated — the
  claim is not that nothing is exotic in 4D, but that the 4D exotic stress
  is the geometric projection of an ordinary 5D bulk, requiring no exotic
  brane matter.
- This is a **static** junction audit; the dynamical version (a
  self-gravitating focused wave crossing a real threshold, PR #166) is the
  motivated frontier, not addressed here.

## Reproduce

```bash
python -m experiments.closure_ledger.israel_junction_weyl_split_probe
# Verdict: BULK_WEYL_SUPPLIES_EFFECTIVE_EXOTIC_SIGMA_RICCI_FLAT_VACUUM_BRANE
```

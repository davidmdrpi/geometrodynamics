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

## The split — what is, and is not, established

The effective stress is **100% of the bulk-Weyl FORM** (this is computed).
But the *attribution* — that this form **is** a bulk-Weyl projection with
no on-brane exotic matter — is **consistent-with, not proven**. The split
is honestly a split of necessary vs pending:

**Necessary conditions — met:**

- `R = 0` (Ricci-flat). By Shiromizu–Maeda–Sasaki, a vacuum brane
  (`T_brane = 0`) obeys `G_μν = −E_μν`, which *forces* `R = 0`. Satisfied.
- A **negative** tidal charge (`Q = −r_s²`, `ρ_eff < 0`). The traceless
  `r⁻⁴` form is shared by (a) a projected bulk Weyl tensor and (b) a real
  on-brane Maxwell field (Reissner–Nordström) — and **only the sign
  distinguishes them**: a real brane gauge field gives `ρ_EM > 0`, whereas
  here `ρ_eff < 0`. So the stress is **not** a real on-brane gauge field —
  provided BAM carries **no fundamental brane gauge field** that would
  force reading (b).

**Sufficient step — pending:**

- The explicit 5D embedding (a 5D Tangherlini bulk whose Weyl projection
  sources exactly this `E_μν`). The Dadhich/Bronnikov–Kim construction
  exists for tidal-charge metrics and is **cited**, not re-solved here as a
  5D boundary-value problem.

**The honest line.** Throat stress is the tidal-charge / bulk-Weyl form;
on-brane exotic matter is **avoidable if** the 5D embedding sources `E_μν`
and BAM has no fundamental brane gauge field — **necessary conditions met,
5D derivation pending**. This is a *consistent-with*, the program's
narrowest and most closable gap, not a proof.

## What this does NOT establish

- It does **not** prove the stress *is* a bulk-Weyl projection — only that
  it is consistent with one and the necessary conditions hold.
- It does **not** evade the **horizon at `f = 0`**: that locus is null and
  degenerate, and the surgical surface term merely **vanishes** there —
  relocating σ, not removing it.
- The 4D WEC is genuinely violated; the claim is only that the 4D exotic
  stress *could be* the projection of an ordinary 5D bulk, not that nothing
  is exotic in 4D.
- This is a **static** junction audit; the dynamical version (a
  self-gravitating focused wave crossing a real threshold, PR #166) is the
  motivated frontier, not addressed here.

## Reproduce

```bash
python -m experiments.closure_ledger.israel_junction_weyl_split_probe
# Verdict: THROAT_IS_TIDAL_CHARGE_BULK_WEYL_FORM_NECESSARY_CONDITIONS_MET_5D_PENDING
```

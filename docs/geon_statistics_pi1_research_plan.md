# π₁ of the two-mouth configuration space and FR-homotopy survival for the Pin⁻ throat (PR #171)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The gap this closes

PR #170 derived the throat's Fermi statistics but **cited** the
Finkelstein–Rubinstein (FR) homotopy — and the standard FR result is for
**orientable** defects. The BAM throat mouth is the non-orientable `RP²`
with a **Pin⁻** structure, so the orientable FR result does not transfer
for free. The correct framework is **geon statistics** (Friedman–Sorkin;
Sorkin; Aneziris–Balachandran–Bourdeau–Jo–Ramadas–Sorkin; Dowker–Sorkin),
where a topological lump's statistics is a representation of `π₁` of the
configuration space (equivalently the mapping class group), and the
spin–statistics correlation is a **theorem with hypotheses — known to fail
for some geons**. This probe computes `π₁` of the two-mouth configuration
space and checks whether the −1 survives a Pin⁻ mouth, against that
literature.

## What is computed

- **`π₁` of the two-mouth configuration space** (model): the exchange `σ`
  with `σ² = e` (in ≥3 spatial dimensions the configuration-space
  fundamental group is the symmetric group — **no braiding**, so the only
  statistics are the ±1 representations); the per-geon 2π rotation `R_i`;
  and, because each mouth is the non-orientable `RP²`, an
  **orientation-reversing loop `τ_i`** absent for orientable mouths.
- **The single geon is spinorial**: the 2π rotation acts on its Pin⁻ spinor
  as `−I` (4π → `+I`) — the Pin⁻ holonomy of #170 and Friedman–Sorkin's
  "spin-½ from gravity."
- **The new ingredient orientable FR lacks**: the non-orientable exchange
  carries an **orientation reversal — a reflection** — and `RP²` admits
  **Pin⁻ only** (`w₂ ≠ 0` kills Pin⁺), in which a reflection **squares to
  −1**:

  | structure | reflection² | RP² admits? |
  |---|---:|---|
  | Pin⁻ | **−1** | yes (`w₂+w₁²=0`) |
  | Pin⁺ | +1 | no (`w₂≠0`) |

- **Achirality**: a non-orientable geon is its own mirror image, so the
  geon/mirror-geon handedness obstruction the spin–statistics theorem must
  control is trivial — non-orientability *helps*.

## The result

Spinorial (`2π = −1`) + the non-orientable exchange's reflection squaring
to −1 in Pin⁻ + achirality ⟹ in the geon-statistics framework the exchange
is the **−1 irrep of `σ`: a fermion**. The FR homotopy **survives** the
Pin⁻ mouth — now on the correct non-orientable footing, with the orientable
rotation argument replaced by the Pin⁻ reflection² = −1 plus the achiral
geon-statistics theorem.

## Honest scope

- **Conditional** on the **Dowker–Sorkin exchangeability ("slide")
  hypothesis** — that two identical throats can be exchanged by an ambient
  diffeomorphism — which holds for identical asymptotically-flat throats and
  is **cited, not derived** from the full BAM field theory.
- The geon literature contains spin–statistics **violation** examples, so
  this is a genuine check that BAM's `Pin⁻ + achiral + exchangeable` mouth
  **passes**, not an automatic result.
- The remaining honest gap is that hypothesis (and the BAM field-theory
  mapping class group, modeled topologically here), **not** the spinor sign
  or the reflection algebra.

## What changes vs #170

The orientable Finkelstein–Rubinstein citation is replaced by the correct
non-orientable geon-statistics framework; the Fermi conclusion survives,
now with the Pin⁻ reflection² = −1 as the explicit new ingredient and the
exchangeability hypothesis made explicit.

## Citations

- J. L. Friedman & R. D. Sorkin, *Spin ½ from gravity*, Phys. Rev. Lett.
  **44**, 1100 (1980).
- C. Aneziris, A. P. Balachandran, M. Bourdeau, S. Jo, T. R. Ramadas,
  R. D. Sorkin, *Aspects of spin and statistics for gravitational geons*,
  Int. J. Mod. Phys. A **4**, 5459 (1989).
- H. F. Dowker & R. D. Sorkin, *A spin-statistics theorem for certain
  topological geons*, Class. Quantum Grav. **15**, 1153 (1998).

## Reproduce

```bash
python -m experiments.closure_ledger.geon_statistics_pi1_probe
# Verdict: FR_HOMOTOPY_SURVIVES_PIN_MINUS_MOUTH_GEON_FERMI_CONDITIONAL_ON_EXCHANGEABILITY
```

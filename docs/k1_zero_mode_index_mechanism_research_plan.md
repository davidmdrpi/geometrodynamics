# The index mechanism: a Pin/Dirac zero mode for the k=1 sector (PR #195)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity.

## The question #194 posed

#194 proved the surrogate's electron near-zero is **dialed** — an
unprotected cancellation, with no symmetry, index, seesaw, or attractor —
and posed the mechanism question sharply: is there structure in the
throat geometry that could pin a k=1 zero mode (an index, a chiral
grading of the winding sectors, a geometric identity tying the base
action to the transport repulsion)?

**Answer: yes — and with no new ingredients.**

## The mechanism

BAM throats are **Pin⁻** (the #183/#188 grading: `T = iσ_y`, `T² = −I`) —
the throat's internal mode is a *spinor*, not a scalar. On the #193
sector reduction, the winding-k mode lives on the base S² coupled to a
monopole of charge `q = k/2`. For the spinor the operator is the **Dirac
operator with monopole charge**, and the Atiyah–Singer index theorem
gives exactly `2q = k` chiral zero modes in sector k:

> sectors {1, 3, 5} → {1, 3, 5} zero modes, all one chirality,
> **pinned at E = 0 by topology — not by tuning**.

The SUSY decomposition makes it computable with the validated #193
monopole solver:

```
D²₊ = L_{q−½} − (q − ½)     (kernel: 2q modes at l = q−½)
D²₋ = L_{q+½} + (q + ½)     (gapped: lowest 2q + 1)
```

The "geometric identity tying the diagonal to the repulsion" that #194
asked for **is** the factorization `D² = A†A`: for the spinor, the
would-be diagonal action and the repulsion are the two pieces of a
perfect square, guaranteed to cancel on the kernel. The scalar has no
factorization — hence the surrogate's dialing.

## The index, verified

| k | zero modes | predicted | worst residual | 1st excited (+) | lowest (−) | exact gap 2q+1 |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1.5e-10 | 2.000000 | 2.000000 | 2 |
| 3 | 3 | 3 | 2.9e-08 | 4.000000 | 4.000000 | 4 |
| 5 | 5 | 5 | 9.3e-08 | 5.999999 | 6.000000 | 6 |

The nonzero towers match the exact Dirac spectrum `λ² = (j+½)² − q²`;
the chirality asymmetry (no `j = q−½` level on the − side) *is* the
index. The count coincidence — {1,3,5} zero modes in the {1,3,5}
generation depths — is noted, not built upon.

## The protection, certified (the discriminator vs #194)

Because `D²₊ ⪰ 0`, a variational residual gives a **two-sided pin** on
the lowest eigenvalue — no eigensolve trust required:

- **Gauge wobble** (flux-preserving, ε = 0.05): the deformed zero mode is
  constructed in closed form (`v₀·e^{−ε∫f}`) and certified — the zero
  energy is pinned in **[0, 6×10⁻¹⁰]** while the wavefunction deforms.
  Energy-pinned, wavefunction-mobile: the defining signature of index
  protection.
- **Metric deformation** of the base (in 2D every metric deformation is
  conformal): `ker D` is conformally rigid — the same zero mode is
  annihilated on the deformed metric, certificate **9×10⁻¹⁶**,
  Ω-independent.
- **The scalar contrast**: on the same deformed metric the scalar sector
  ground moves 1.5000 → 1.4012 (O(ε)) — **8 orders of magnitude**
  between the pinned spinor zero and the moving scalar level. The #192
  Berger squash moved the surrogate's electron because the surrogate is
  a scalar-type cancellation; the spinor zero mode would not have moved.
- **Flux-change control**: the count jumps only with a flux quantum
  (1 → 3 → 5) — topology moves the index, smooth deformation cannot.
  The #183 algebra dichotomy, realized at the level of the spectrum.

## The natural mass (two-mouth pairing)

- **Within one mouth, a first-order mass lift is forbidden** by angular
  momentum: the zero modes sit at `j = q−½` and the opposite-chirality
  tower starts at `j = q+½` — no partner exists. A Weyl-like protection.
- **The lift requires the throat's two mouths** (opposite winding ±k ⟹
  opposite monopole charge ⟹ opposite chirality — a genuine Dirac pair;
  the BAM wormhole supplies exactly this). For k=1 both mouths' zero
  modes are the l=0 constants; the antipodal pairing overlap is
  `o = 1.000`, and the lifted level is

  ```
  |E_e| = ε · o     (linear, multiplicative, sign-stable)
  ```

  't Hooft-natural: ε → 0 restores the two independent chiral zero
  modes, so E_e is radiatively protected. **Contrast #194:** the
  surrogate's E_e was the *difference* of two O(7) numbers (sign-flips
  under ±2%, Δ = 74.7); here it is a *product* — small because the
  mouths couple weakly, not because two big numbers cancel.

## Honest scope

- **Established:** the throat geometry, with the field content BAM
  already requires, contains an index mechanism that pins a k=1 zero
  mode, energy-rigid under smooth deformations, with a natural mass from
  the two-mouth pairing.
- **Not established:** the lepton mass ladder is *not* re-derived — the
  locked surrogate's {1,3,5} energies and observed ratios are untouched.
  Identifying the surrogate's k=1 state with the spinor zero mode,
  computing the mouth coupling ε from the throat overlap machinery
  (#185/#190), and rebuilding the ladder on the Dirac tower is the
  follow-up.
- The base is the round (and conformally deformed) S²; the full Berger
  3-geometry Dirac spectrum is a further check.

## Reproduce

```bash
python -m experiments.closure_ledger.k1_zero_mode_index_mechanism_probe
# Verdict: K1_ZERO_MODE_IS_INDEX_PROTECTED_THE_PIN_DIRAC_STRUCTURE_SUPPLIES
#          _THE_MECHANISM_THE_SCALAR_SURROGATE_LACKS
```

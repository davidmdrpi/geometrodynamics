# Probe I — what caps BAM's nonlocality at Tsirelson? (PR #236)

_Run 2026-08-05T00:26:13.876618+00:00 · 8/8 PASS_

**Q. The bridge escapes Bell. What stops it at 2√2 rather than the algebraic 4?**

**A. Not the geometry. `B² = I`.**

## The ceiling is real

| quantity | value |
|---|---:|
| max CHSH over all 8 holonomies and all preparations | 2.828427125 |
| Tsirelson | 2.828427125 |
| excess | +8.9e-16 |

## Where each ingredient comes from

| ingredient | supplied by |
|---|---|
| the pair state (the singlet) | GEOMETRY |
| the nonlocality (violating Bell at all) | GEOMETRY |
| the measurement settings | GEOMETRY |
| THE CEILING (why not 4?) | NOT THE GEOMETRY |

Locality is **not** the binding constraint: the mouth projectors commute to 4.7e-113, but dropping commutativity entirely still gives 2.828427095, because `(B₀+B₁)² + (B₀−B₁)² = 4I` for any dichotomic observables (residual 5.3e-15).

## No-signaling is equal-time only

| `gb` | `dt` | Bob's marginal shift |
|---:|---:|---:|
| 0.8 | 0.0 | 0.000e+00 |
| 0.8 | 0.5 | 2.273e-02 |
| 0.8 | 2.0 | 1.019e-01 |
| 0.8 | 8.0 | 2.315e-01 |
| 0.2 | 0.0 | 0.000e+00 |
| 0.2 | 0.5 | 5.262e-03 |
| 0.2 | 2.0 | 1.344e-02 |
| 0.2 | 8.0 | 3.600e-02 |
| 0.0 | 0.0 | 0.000e+00 |
| 0.0 | 0.5 | 1.499e-15 |
| 0.0 | 2.0 | 1.166e-15 |
| 0.0 | 8.0 | 1.499e-15 |

## Verdict

**BAM_REACHES_TSIRELSON_EXACTLY_AND_NEVER_EXCEEDS_IT_BUT_THE_CEILING_IS_NOT_GEOMETRIC_IT_FOLLOWS_FROM_DICHOTOMIC_READOUT_ALONE_SO_THE_IMPORTED_QUANTIZATION_REDUCES_TO_B_SQUARED_EQUALS_I_AND_NO_SIGNALING_IS_EQUAL_TIME**

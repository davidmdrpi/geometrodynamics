# The v3+CP joint re-lock: the v4 candidate lock (PR #163)

> **Framing (to avoid a category error).** QFT on the *fixed classical*
> throat geometry — geometry → fields, **not** quantum gravity. The re-lock
> expresses the #161 target state through a structured parameter set; the
> library migration is staged separately.

The flagged #161 successor: realize the tabulated target state in the
model's own kind of parameters — and the realization comes with a sharp
structural finding.

## The minimal-law no-go (exact)

The v3 off-diagonal law (shared transport `−t·e^{−α·dk}`, `dk = max`, plus
the one existing targeted coupling `η_35^−`) enforces two equalities the
targets **break**:

1. **Partition symmetry of the transport magnitude**: `H₊[12] = H₋[12]` by
   construction (machine-exact at the lock); the targets need ×1.287 vs
   ×1.832 — a partition split of ratio **1.424**.
2. **The dk = max degeneracy**: `H[13] = H[23]` within each block; the
   minus-block d–b element must be enhanced **×1.996**.

Meanwhile, in the up block — where the data permits — the law survives
**exactly** (deviations 5×10⁻⁶). The breaking pattern is the partition
asymmetry, concentrated on the minus block's d-row: precisely where #155
located the physical mixing, and the sector where the v3 lock already
carries targeted couplings (χ_k3, η_35^−) for the s/t outliers.

## The v4 candidate lock

**Element level**: the twelve tabulated target numbers + the derived phases
`e^{±i(π/k₅)·dk}` reproduce the v3 masses exactly (eigenvalues to 1e-15) and
all nine flavor-CP observables at ≤ 1% — the first complete flavor state of
the program in one parameter set.

**Structured level**: the v3 law + exactly three new targeted couplings +
one retune:

| coupling | value | element | factor vs v3 |
|---|---|---|---|
| η_12^+ (new) | −0.1018 | H₊[12] | ×1.287 |
| η_12^− (new) | −0.2953 | H₋[12] | ×1.832 |
| η_13^− (new) | −0.2671 | H₋[13] | ×1.996 |
| η_35^− (retune) | 5.0 → 5.586 | H₋[23] | ×1.111 |

plus diagonal retunes within the existing diagonal law's reach (up ±0.002;
down +0.113/−0.065/−0.048). Three of the four touched couplings live in the
minus block — the extension *continues* the lock's own pattern.

## The counting (honest)

+3 parameters buy +5 independent observables (V_us, V_cb, V_ub, β, γ; the
other four of the nine follow from unitarity and the derived phase): **net
predictive surplus +2** relative to v3 — and the entire CP sector costs
**zero parameters** (φ_h = π/k₅ derived, #158–#160). The six masses are
inherited from the v3 calibration exactly. The #150 budget is unchanged.

## The staged migration (the follow-up PR)

1. Three new `QuarkParams` fields (`eta_12_plus`, `eta_12_minus`,
   `eta_13_minus`), following the `eta_k3k5_minus` pattern.
2. The complexified same-partition transport with φ_h = π/k₅ as the
   structural default (the #158 relocation, in code).
3. The `LOCKED_QUARK_PARAMS` update (η_35^− retune + diagonal-law retune).
4. Regression re-baseline of the #155–#162 probes against the new lock.

This probe deliberately leaves `geometrodynamics/qcd` untouched.

## Reproduce

```bash
python -m experiments.closure_ledger.v4_relock_realization_probe
# Verdict: V4_LOCK_REALIZED_THREE_TARGETED_COUPLINGS_MINIMAL_LAW_NO_GO_PARTITION
```

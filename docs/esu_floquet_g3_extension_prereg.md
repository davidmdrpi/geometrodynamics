# Prospective extension of gate G3 (integrator convergence) for #310

Date: 2026-09-26. Freeze [`4e65c3e`](esu_floquet_refocusing_prereg.md) is
unchanged. The results commit `3251f72` retains its outcome: V and S are
**UNRESOLVED** under the frozen G3, and that record is not relabelled.
Publish this addendum before implementing or running the extension.

## Why

The frozen RK4 ratio test compared 2^16- and 2^17-step traces against the
primary DOP853 map. At high degree, both differences reach the primary's own
error, about 1e-11 to 1e-12, so their ratio measures reference noise rather
than RK4 truncation. The agreement conditions passed everywhere, with
maxima of 3.4e-12 (RK4) and 9.1e-10 (secondary DOP853). The extension asks
the question G3 intended to ask: do independent fourth-order integrations
converge to the primary maps at the fourth-order rate, in a regime where
truncation error dominates?

## Test (all sectors T, V, S; n = 2..80; archived primary maps unchanged)

1. Compute RK4 fundamental matrices over [0, pi] with 2^10, 2^11 and 2^12
   fixed steps. Use the same vectorised routine and right-hand sides, with
   no change to any equation.
2. Let e_k = |tr M_RK4(2^k) - tr M_primary| for k = 10, 11, 12.
3. Evaluate each successive ratio e_k/e_(k+1) only when e_k > 1e-9. That
   floor is about 100x the observed reference error. Each evaluated ratio
   must lie in [8, 32].
4. Non-vacuity: every sector must have at least one evaluated ratio at some
   n >= 20.
5. The frozen agreement conditions remain in force: 2^17-step RK4 and
   secondary DOP853 each within 1e-7 of the primary trace, as already
   measured.

If a sector passes, compute extended verdicts from the unchanged archived
maps with the frozen classifier and thresholds, and name them
`<X>_STABILITY_EXT`, `<X>_REFOCUSING_EXT` and `V_WKB_PREDICTION_EXT`. The
frozen classifier is applied as is, including its known weakness for S
(section 6 of the results). If a sector fails, it stays UNRESOLVED and the
failure is reported. No further tolerance or step-count search follows.

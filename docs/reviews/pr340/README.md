# Review assets for PR #340: two moving antipodal mouth pairs

All scripts use PR #340's own `MovingMouthPairs` (`geometrodynamics/waves/moving_mouth_pairs.py`
on that branch) with its default parameters (M = 0.1). Run them from a checkout of #340:

    PYTHONPATH=. xvfb-run -a python docs/reviews/pr340/render_long_horizon_pairs.py --output long_horizon
    PYTHONPATH=. python docs/reviews/pr340/review_controls.py        # writes controls.json, long_*.npz
    python docs/reviews/pr340/long_horizon_figure.py                  # reads long_*.npz
    python docs/reviews/pr340/window_impulses.py                      # reads long_directional.npz
    PYTHONPATH=. python docs/reviews/pr340/long_horizon_checks.py    # dt/2, cutoff 40, displaced pair
    PYTHONPATH=. python docs/reviews/pr340/displaced_pair_windows.py
    PYTHONPATH=. python docs/reviews/pr340/amplitude_scaling.py      # follow-up review (650e825)
    PYTHONPATH=. python docs/reviews/pr340/symmetric_stop.py         # follow-up review (650e825)

## The movie

`long_horizon_pairs.gif` (final frame `long_horizon_end.png`) runs the same model and
parameters to 4π + 1.2 instead of π + 1.2. Only the horizon changes:

- **Directional packets:** the pairs travel 29° and 49°.
- **Symmetric packets:** the pairs converge. The panel freezes at t = 12.30, the last
  output before mouths of different handles reach the PR protocol's 2·0.3 + 0.05 rad
  disjointness limit. This is a display stop: the script integrates the full horizon, and
  states past the limit are not valid finite-mouth dynamics. `symmetric_stop.py` stops the
  integration itself; the last admissible state is at t = 12.312.

Energy is conserved to 7e-6 relative and total angular momentum to 1e-14.

The drawing differs from the PR renderer in two ways, both needed to show large motion:

- **A sphere fixed in space.** The mouths are cut out of it each frame, and the field is
  sampled at its fixed points. The PR's collar warp folds the sphere mesh once a mouth
  moves more than about 12°: the minimum orientation ratio is 0.05 at 11° and −0.04 at 12°.
- **Field drawn as colour only,** with no radial display. The only geometric motion on
  screen is therefore the solved motion of the mouths.

The necks are #338's neutral drawings, rotated rigidly with their pair axis. A Klein handle
crosses the sphere wall only because it cannot be drawn in 3D without crossing.

## Follow-up review scripts

- `amplitude_scaling.py` scales the prepared field by λ = 1, ½ and ¼. It separates the bilinear
  part of the single-packet nonadditivity residual from its trajectory and backreaction part.
- `symmetric_stop.py` stops the symmetric run at the first integration step past the
  disjoint-mouth limit, rather than freezing the display as the movie does.

## The controls

`long_horizon_controls.png` plots pair momenta and true angles to 4π + 1.2 for both
preparations. The scripts above reproduce every number in the review.

This is a conservative scalar comparator on S² with operator −Δ + ¼. It is not a GR
calculation and makes no quantization claim.

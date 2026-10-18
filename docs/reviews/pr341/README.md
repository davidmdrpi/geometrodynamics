# Review assets for PR #341: leading mouth exchange

All scripts use PR #341's own `MovingMouthPairs` model and preparation, with default
parameters (M = 0.1, cutoff 32). Run them from this folder with a checkout of #341 on
`PYTHONPATH`:

    PYTHONPATH=<pr341> python fixed_axis_mechanism.py   # coefficient + split by term, mouth and window
    PYTHONPATH=<pr341> python mechanism_figure.py       # mechanism.png (after the line above)
    PYTHONPATH=<pr341> python small_lambda.py           # moving model at lambda = 1/32, 1/64
    PYTHONPATH=<pr341> python start_offset.py           # coefficient vs packet start time delta
    PYTHONPATH=<pr341> python start_offset_moving.py    # moving-model cross-check at delta = 1.2

- **`fixed_axis_mechanism.py`** is an independent fixed-axis linear response. It is a Strang
  integrator with exact free flows, separate from the PR's DOP853 solver.
  - At dt = 0.002 and 0.001, extrapolated to dt → 0, it reproduces the PR's leading
    coefficients to 3e-10 relative.
  - It splits pair 1's mixed y torque by term (e^z·D_yU^x and e^x·D_yU^z), by mouth (±z)
    and by window.
- **`small_lambda.py`** evaluates the PR's moving model directly at λ = 1/32 and 1/64, with
  no extrapolation.
- **`start_offset.py`** backpropagates the prepared packets by π + δ, so that the mouths
  rest before a complete source focus. It reports the leading coefficient for
  δ ∈ {0, 0.3, 0.6, 0.9, 1.2}; δ = 0 is the PR's preparation.

This is a conservative scalar comparator on S² with operator −Δ + ¼. It is not a GR
calculation and makes no quantization claim.

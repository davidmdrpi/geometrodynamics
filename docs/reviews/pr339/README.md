# Review asset for PR #339: visible momentum transfer

`momentum_transfer.gif` is rendered by `render_momentum_transfer.py` from PR
#339's own model (`geometrodynamics/waves/antipodal_kick.py` on that branch),
run from a checkout of #339:

    PYTHONPATH=. python docs/reviews/pr339/render_momentum_transfer.py --output momentum_transfer.gif

The equations, integrator, waves and coupling are unchanged. The only change is
the attachment mass, which is M = 0.1 instead of the default 0.3, so the solved
motion is visible without magnifying positions. The horizon is extended to
3π + 1.2, so the movie covers the first focus at π, the coasting that follows,
and the next focus at 3π.

Bead positions, trails and momentum arrows come from the solved X(t) and L(t).
The ledger chart plots the field's and the attachment's angular momentum about
z; their sum is conserved to about 6e-15. Field height and colour use a fixed
display gain.

This is a conservative scalar comparator on S² with operator −Δ + ¼. It is not
a GR calculation and makes no quantization claim.

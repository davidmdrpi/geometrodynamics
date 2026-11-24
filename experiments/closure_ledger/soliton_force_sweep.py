"""Exploratory response to PR #361 review, preserving the registered probe."""
import json

from .soliton_current import PLUS_WEIGHT, run


def measure(c, force, fine=False):
    grid = dict(n=2048, dt=.0005) if fine else {}
    rows = run(c=c, force=force, measure_components=True, **grid)
    last = rows[-1]
    result = dict(c=c, force=force, fine=fine, peak=last["spectral_peak"],
                  median=last["number_median"], centroid=last["charge_centroid"],
                  component_centroids=last["component_centroids"],
                  outside_at_start=rows[0]["outside_three"], outside_at_end=last["outside_three"],
                  max_peak_median_gap=max(abs(r["spectral_peak"]-r["number_median"]) for r in rows))
    if c is None:
        # The uniform-gradient, initially resting NLS solution. The drive here
        # is tapered only in the tails, so retain the measured discrepancy.
        result["uniform_force_prediction"] = force*(2*PLUS_WEIGHT-1)
        result["centroid_residual"] = last["charge_centroid"]-result["uniform_force_prediction"]
    return result


def probe():
    return dict(component_order=["+1", "-1"], outside_radius=3.,
                primary=[measure(c, force) for c in (8., None)
                         for force in (0., .125, .25, .5, 1.)],
                refinement=[measure(c, force, True) for c in (8., None) for force in (.125, 1.)])


if __name__ == "__main__":
    print(json.dumps(probe(), indent=2, allow_nan=False))

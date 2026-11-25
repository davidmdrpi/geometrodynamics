"""Exploratory transfer of #363's balanced fluence sweep to its shared-ring successor."""
import json

import numpy as np

from .closed_loop_records import base, evolve


def probe():
    q, p = base.preparation(512, 363, .02)
    rows = []
    for ratio in np.arange(1., 2.26, .125):
        pair = [ratio/2, ratio/2]
        for closed, delay in ((False, 4.), (True, 3.), (True, 4.), (True, 6.)):
            row, _, _ = evolve(pair, q, p, closed=closed, delay=delay)
            rows.append(dict(energy_ratio=float(ratio), closed=closed, delay=delay,
                             persistent=row["persistent_counts"], endpoint=row["endpoint_counts"],
                             exclusive=row["exclusive"], budget=row["max_energy_residual"]))
    return rows


if __name__ == "__main__":
    print(json.dumps(probe(), indent=2, allow_nan=False))

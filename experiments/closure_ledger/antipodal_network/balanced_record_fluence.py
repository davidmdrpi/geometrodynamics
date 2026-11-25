"""Post-registration control: increase total fluence while holding the split balanced."""
import json

import numpy as np

from .two_detector_records import evolve, preparation


def probe():
    q, p = preparation(512, 363, .02)
    gains = np.arange(1., 2.26, .125)
    fractions = np.column_stack((gains/2, gains/2))
    results = []
    for mu in (0., .25):
        rows, labels, residual = evolve(fractions, q, p, mu, .025)
        results.append(dict(mu=mu, max_budget_residual=residual,
                            rows=[dict(energy_ratio=float(gain), counts=row["counts"],
                                       exclusive=float(np.mean((codes == 1) | (codes == 2))))
                                  for gain, row, codes in zip(gains, rows, labels)]))
    return results


if __name__ == "__main__":
    print(json.dumps(probe(), indent=2, allow_nan=False))

"""Render the wave-mouth-network visualization and its controls.

Every panel is drawn from an integrated history of one Hamiltonian; see
``geometrodynamics/viz/wave_mouth_network.py`` and
``docs/wave_mouth_network_visualization.md``.
"""

import argparse
import json
from pathlib import Path

from geometrodynamics.viz import wave_mouth_network as w

RUNS = {
    "wave_mouth_network": dict(
        amplitude=w.AMPLITUDE, stored=(4.0, 0.0, 0.0),
        title="A wave moves mouths, collapses a loaded throat into a larger "
              "wave, and the ledger closes"),
    "wave_mouth_network_subthreshold": dict(
        amplitude=w.SUBTHRESHOLD, stored=(4.0, 0.0, 0.0),
        title="Control: a weaker wave still moves the mouths; the throat "
              "stays open"),
    "wave_mouth_network_unloaded": dict(
        amplitude=w.AMPLITUDE, stored=(0.0, 0.0, 0.0),
        title="Control: same wave, no stored energy; nothing to release"),
    "wave_mouth_network_no_wave": dict(
        amplitude=0.0, stored=(4.0, 0.0, 0.0),
        title="Control: no wave; nothing moves and nothing collapses"),
}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        default=Path("docs/images/wave_mouth_network"))
    args = parser.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = {}
    for name, cfg in RUNS.items():
        hist = w.scenario(cfg["amplitude"], cfg["stored"])
        report[name] = dict(cfg, **w.render(hist, args.output_dir / f"{name}.png",
                                            cfg["title"]))
        print(name, json.dumps({k: report[name][k] for k in (
            "energy_max_relative_error", "collapse_times", "wave_energy_initial",
            "wave_energy_final")}))
    (args.output_dir / "summary.json").write_text(
        json.dumps(report, indent=2, default=float) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

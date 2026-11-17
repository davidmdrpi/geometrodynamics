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
        title="A wave moves mouths and releases a loaded throat's stored "
              "energy into the field; the ledger closes"),
    "wave_mouth_network_subthreshold": dict(
        amplitude=w.SUBTHRESHOLD, stored=(4.0, 0.0, 0.0),
        title="Control: a weaker wave still moves the mouths; the throat "
              "stays open"),
    "wave_mouth_network_unloaded": dict(
        amplitude=w.AMPLITUDE, stored=(0.0, 0.0, 0.0),
        title="Control: same wave, no stored energy; nothing to release"),
    "wave_mouth_network_no_wave": dict(
        amplitude=0.0, stored=(4.0, 0.0, 0.0),
        title="Control: no wave; nothing moves and no transition"),
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
            "energy_max_relative_error", "transition_times",
            "threshold_crossings", "field_energy_initial",
            "field_energy_final", "far_field_energy_final")}))
    # Control for the detector: a stable throat that swings below 0.5.
    swing = w.summary(w.scenario(1.0, (0.0, 0.0, 0.0), duration=40.0))
    report["detector_control"] = {
        "amplitude": 1.0, "stored": [0.0, 0.0, 0.0], "duration": 40.0,
        **{k: swing[k] for k in ("transition_times", "threshold_crossings",
                                 "min_coordinate", "final_coordinates")}}
    (args.output_dir / "summary.json").write_text(
        json.dumps(report, indent=2, default=float) + "\n")
    # Energy flux through fixed boundaries around mouth 1, every run.
    scalars = ("box", "step", "duration", "net_outward_final",
               "net_outward_left_final", "net_outward_right_final",
               "gross_outward", "gross_inward", "inside_initial",
               "inside_final", "box_field_energy_final",
               "balance_max_residual", "mouth_margin")
    fluxes = {}
    for name, cfg in RUNS.items():
        flux = w.flux_ledger(cfg["amplitude"], cfg["stored"])
        fluxes[name] = {k: flux[k] for k in scalars}
        if name == "wave_mouth_network":
            fluxes[name]["net_outward_minimum"] = float(flux["net_outward"].min())
            w.render_flux(flux, args.output_dir / "wave_mouth_network_flux.png",
                          report[name]["transition_times"][0])
    half = w.flux_ledger(h=w.STEP / 2, sample=40)
    fluxes["wave_mouth_network_half_step"] = {k: half[k] for k in scalars}
    (args.output_dir / "flux.json").write_text(
        json.dumps(fluxes, indent=2, default=float) + "\n")
    print("flux", json.dumps({k: v["net_outward_final"]
                              for k, v in fluxes.items()}))
    validation = w.validate()
    (args.output_dir / "validation.json").write_text(
        json.dumps(validation, indent=2) + "\n")
    print("validation", json.dumps(validation))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

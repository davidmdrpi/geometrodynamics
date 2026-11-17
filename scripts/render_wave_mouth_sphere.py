#!/usr/bin/env python
"""Render the sphere-with-Alice-handles figures from the solved model."""
import argparse
import json
from pathlib import Path

from geometrodynamics.viz import wave_mouth_sphere as v
from geometrodynamics.waves import wave_mouth_sphere as m

RUNS = {
    "wave_mouth_sphere": dict(
        amplitude=m.AMPLITUDE, stored=((2.0, 0.0), (0.0, 0.0)), title=None),
    "wave_mouth_sphere_subthreshold": dict(
        amplitude=m.SUBTHRESHOLD, stored=((2.0, 0.0), (0.0, 0.0)),
        title="Control: a weaker wave moves the mouths; the throat stays open"),
    "wave_mouth_sphere_unloaded": dict(
        amplitude=m.AMPLITUDE, stored=((0.0, 0.0), (0.0, 0.0)),
        title="Control: same wave, no stored energy; no transition"),
    "wave_mouth_sphere_no_wave": dict(
        amplitude=0.0, stored=((2.0, 0.0), (0.0, 0.0)),
        title="Control: no wave; nothing moves and no transition"),
}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        default=Path("docs/images/wave_mouth_sphere"))
    parser.add_argument("--no-gif", action="store_true")
    args = parser.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report, scale = {}, None           # controls share the documented scale
    for name, cfg in RUNS.items():
        model, data = m.scenario(cfg["amplitude"], cfg["stored"])
        report[name] = dict(
            amplitude=cfg["amplitude"], stored=cfg["stored"],
            **v.render(model, data, args.output_dir / f"{name}.png", cfg["title"],
                       scale=scale))
        scale = scale or report[name]["field_scale"]
        if name == "wave_mouth_sphere" and not args.no_gif:
            v.animate(model, data, args.output_dir / f"{name}.gif")
        print(name, json.dumps({k: report[name][k] for k in (
            "energy_max_relative_error", "transition_times",
            "field_energy_initial", "field_energy_final",
            "pair_displacements")}))
    _, half = m.scenario(h=m.STEP / 2, sample=40)
    model, _ = m.scenario(duration=1.0)
    report["half_step"] = {k: v.summary(model, half)[k] for k in (
        "energy_max_relative_error", "field_energy_final", "transition_times",
        "pair_displacements")}
    # Energy flux through a fixed circle around the loaded mouth, every run.
    scalars = ("centre", "radius", "duration", "step", "net_outward_final",
               "net_outward_minimum", "gross_outward", "gross_inward",
               "inside_initial", "inside_final", "cap_field_energy_final",
               "rotor_work_final", "balance_max_residual",
               "balance_max_residual_with_tails", "tail_work_final",
               "net_flux_rise", "net_flux_fall", "edge_points",
               "loaded_mouth_margin", "other_mouth_margin")
    fluxes = {}
    for name, cfg in RUNS.items():
        flux = m.cap_flux_ledger(cfg["amplitude"], cfg["stored"])
        fluxes[name] = {k: flux[k] for k in scalars}
        if name == "wave_mouth_sphere":
            fluxes[name]["net_outward_maximum"] = float(flux["net_outward"].max())
            v.render_flux(flux, args.output_dir / "wave_mouth_sphere_flux.png",
                          report[name]["transition_times"][0])
    half_flux = m.cap_flux_ledger(h=m.STEP / 2, sample=40)
    fluxes["half_step"] = {k: half_flux[k] for k in scalars}
    quarter = m.cap_flux_ledger(h=m.STEP / 4, sample=80)
    fluxes["quarter_step"] = {k: quarter[k] for k in scalars}
    fine_edge = m.cap_flux_ledger(edge_points=416)
    fluxes["edge_416"] = {k: fine_edge[k] for k in scalars}
    (args.output_dir / "flux.json").write_text(
        json.dumps(fluxes, indent=2, default=float) + "\n")
    print("flux", json.dumps({k: f["net_outward_final"]
                              for k, f in fluxes.items()}))
    (args.output_dir / "summary.json").write_text(
        json.dumps(report, indent=2, default=float) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

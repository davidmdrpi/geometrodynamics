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
    (args.output_dir / "summary.json").write_text(
        json.dumps(report, indent=2, default=float) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

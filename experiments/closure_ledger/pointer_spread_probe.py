"""Publicly preregistered finite-pointer-spread experiment; nonzero on failure."""

import argparse
from dataclasses import asdict
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from geometrodynamics.bulk import pointer_spread as spread
from geometrodynamics.bulk.closure_equilibrium import rotor_potential
from geometrodynamics.bulk.mouth_spin_frame import qmul
from geometrodynamics.bulk.source_readout import source_sphere, record_statistics


def ode_reference(x, p, P, model=spread.PointerModel()):
    """Independent embedded Hamilton ODE, transport ODE and external work."""
    m = np.asarray(model.axis) / np.linalg.norm(model.axis)

    def rhs(t, y):
        x, p, pointer_p, C = y[:3], y[3:6], y[7], y[8:12]
        if 0 <= t <= model.duration:
            h = 2 * model.strength / model.duration * math.sin(math.pi * t / model.duration) ** 2
            dh = 2 * model.strength * math.pi / model.duration ** 2 * math.sin(2 * math.pi * t / model.duration)
        else:
            h = dh = 0.
        f = m @ x
        return np.r_[p, -(p @ p) * x - h * pointer_p * (m - f * x),
                     pointer_p / model.mass + h * f, 0.,
                     .5 * qmul(np.r_[0., np.cross(x, p)], C), dh * pointer_p * f]

    initial = np.r_[x, p, 0., P, 1., 0., 0., 0., 0.]
    options = dict(method="DOP853", rtol=2e-12, atol=2e-14)
    pulse = solve_ivp(rhs, (0, model.duration), initial, **options)
    if not pulse.success:
        raise RuntimeError(pulse.message)
    read = pulse.y[:, -1]
    late = solve_ivp(rhs, (model.duration, model.future_time), read, **options)
    if not late.success:
        raise RuntimeError(late.message)
    final = late.y[:, -1]
    energy_change = .5 * (read[3:6] @ read[3:6] - p @ p)
    return {"x_read": read[:3], "p_read": read[3:6], "record_shift": read[6],
            "x_final": final[:3], "p_final": final[3:6], "transport": final[8:12],
            "work_residual": float(abs(energy_change - read[12])),
            "P_change": float(abs(final[7] - P))}


def mechanics_controls():
    rows = []
    rng = np.random.default_rng(607)
    for P in (-.3, 0., .2):
        x = rng.normal(size=3)
        x /= np.linalg.norm(x)
        p = rng.normal(size=3) * .2
        p -= (p @ x) * x
        ref = ode_reference(x, p, P)
        errors = []
        for steps in (32, 64, 128):
            result = spread.evolve(x[None], p[None], [P], steps=steps)
            error = max(float(np.max(np.abs(result[k][0] - ref[k]))) for k in (
                "x_read", "p_read", "record_shift", "x_final", "p_final", "transport"))
            errors.append(error)
        rows.append({"P": P, "errors_32_64_128": errors,
                     "ratios": [errors[0] / errors[1], errors[1] / errors[2]],
                     "work_residual": ref["work_residual"], "P_change": ref["P_change"]})
    return rows


def geometry_controls():
    rng = np.random.default_rng(20260906)
    product_error = max(float(np.max(np.abs(spread.qmultiply(a, b) - qmul(a, b))))
                        for a, b in rng.normal(size=(16, 2, 4)))
    x, _ = spread.source_nodes(7, 20260906, 0.)
    paths = spread.evolve(x, np.zeros_like(x), np.zeros(len(x)))
    error = 0.
    for b in spread.FUTURE_ANALYZERS:
        weights = spread.history_weights(paths, spread.FIXED_ANALYZER, b)
        expected = np.column_stack([rotor_potential(x, sa * spread.FIXED_ANALYZER, -sb * b)
                                    for sa in (-1, 1) for sb in (-1, 1)])
        error = max(error, float(np.max(np.abs(weights["energies"] - expected))))
    return {"quaternion_product_residual": product_error, "stationary_triangle_energy_residual": error,
            "stationary_record_residual": float(np.max(np.abs(paths["record_shift"] - x @ spread.READOUT_AXIS)))}


def preparation_controls():
    """Constant kernel and instrument-preserving rotation, including dynamics."""
    x, p = spread.source_nodes(5)
    P = np.linspace(-.2, .2, len(x))
    R = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    model = spread.PointerModel(axis=(0., 0., 1.))
    first = spread.evolve(x, p, P, model)
    second = spread.evolve(x @ R.T, p @ R.T, P, model)
    w0 = spread.history_weights(first, spread.FIXED_ANALYZER, spread.FUTURE_ANALYZERS[0])["weight"]
    w1 = spread.history_weights(second, spread.FIXED_ANALYZER, spread.FUTURE_ANALYZERS[1])["weight"]
    _, prior = spread.pointer_nodes(.1, 8)
    measures = spread.normalized_measures(np.tile(w0, (8, 1)), prior, w1)
    return {"record_rotation_residual": float(np.max(abs(first["record_shift"] - second["record_shift"]))),
            "weight_rotation_residual": float(np.max(abs(w0 - w1))),
            "record_variance": float(np.var(first["record_shift"])),
            "constant_kernel_residual": float(max(abs(np.sum(measures[k] * .37) - .37)
                                                   for k in ("fixed", "joint", "frozen")))}


def run_probe(progress=None):
    rows, checks, refinements = [], {}, []

    def note(message):
        if progress is not None:
            progress(message)

    def check(name, value):
        checks[name] = bool(value)

    mechanics = mechanics_controls()
    geometry = geometry_controls()
    preparation = preparation_controls()
    check("P4 rotation-covariant blind instrument", preparation["record_rotation_residual"] < 1e-11
          and preparation["weight_rotation_residual"] < 1e-11 and preparation["record_variance"] > .1)
    check("P4 constant response erases information", preparation["constant_kernel_residual"] < 1e-12)
    check("P1 quaternion product", geometry["quaternion_product_residual"] < 1e-11)
    check("P1 stationary triangle bridge", geometry["stationary_triangle_energy_residual"] < 1e-11)
    check("P1 stationary pointer bridge", geometry["stationary_record_residual"] < 1e-11)
    for row in mechanics:
        check(f"P5 ODE at P={row['P']}", row["errors_32_64_128"][-1] < 1e-3)
        check(f"P5 work at P={row['P']}", row["work_residual"] < 1e-10 and row["P_change"] < 1e-12)
        if row["P"] != 0:
            check(f"P5 second-order at P={row['P']}", all(3 < r < 5 for r in row["ratios"]))
    stationary = spread.simulate_spread(0., model=spread.PointerModel(source_sigma=0.))
    baseline = []
    for b in spread.FUTURE_ANALYZERS:
        x, w = source_sphere(spread.FIXED_ANALYZER, b)
        baseline.append(record_statistics(x, w, noise=.15)["tail_probability"])
    check("P1 independent sphere probability bridge", max(abs(stationary["choices"][i]["fixed"]["tail"] - baseline[i])
                                                         for i in range(2)) < 3e-3)
    for sigma in spread.SPREADS:
        for seed in spread.SEEDS:
            row = spread.simulate_spread(sigma, seed=seed)
            rows.append(row)
            note(f"sigma={sigma:g}, seed={seed}: fixed contrast={row['contrasts']['fixed']:.8f}")
    primary = [r for r in rows if r["sigma_pointer"] == .1]
    check("P2 primary contrast exceeds 0.1 in every scramble", all(r["contrasts"]["fixed"] > .1 for r in primary))
    invariant_names = ("sphere_residual", "tangency_residual", "transport_norm_residual",
                       "axial_momentum_residual", "loop_axis_residual")
    check("P5 flow and transport invariants", max(r["diagnostics"][k] for r in rows for k in invariant_names) < 1e-9)
    check("P4 weight and record parity", max(r["diagnostics"][k] for r in rows
          for k in ("weight_parity_residual", "record_antipodal_residual")) < 1e-9)
    check("P4 fixed Gaussian pointer moments", max(max(abs(c["fixed"]["P_mean"]),
          abs(c["fixed"]["P_variance"] - r["sigma_pointer"] ** 2)) for r in rows for c in r["choices"]) < 1e-10)
    check("P4 odd means and smoothed signs", max(max(abs(c[k]["mean"]),
          abs(c[k]["positive_probability"] - .5)) for r in rows for c in r["choices"]
          for k in ("fixed", "joint", "frozen")) < 1e-9)
    check("P5 nonzero primary source recoil", min(r["diagnostics"]["reference_momentum_recoil_rms"] for r in primary) > 1e-3)
    check("P5 zero-P source motion is free", max(r["diagnostics"][k] for r in rows if r["sigma_pointer"] == 0
          for k in ("reference_momentum_recoil_rms", "reference_position_recoil_rms")) < 1e-10)
    for name, kwargs, tolerance in (("steps", {"steps": 128}, 1e-3),
                                    ("Hermite", {"hermite": 32}, 1e-2),
                                    ("source", {"power": 12}, 1e-2)):
        changes = []
        for base in primary:
            refined = spread.simulate_spread(.1, seed=base["seed"], **kwargs)
            changes.append(max(abs(refined["choices"][i][k]["tail"] - base["choices"][i][k]["tail"])
                               for i in range(2) for k in ("fixed", "joint", "frozen")))
            refinements.append({"kind": name, "max_probability_change": changes[-1], "result": refined})
        check("resolution primary " + name, max(changes) < tolerance)
        note(f"primary {name} refinement: max probability change={max(changes):.4g}")
    wide_changes = []
    for base in (r for r in rows if r["sigma_pointer"] == .5):
        refined = spread.simulate_spread(.5, seed=base["seed"], hermite=32)
        change = max(abs(refined["choices"][i][k]["tail"] - base["choices"][i][k]["tail"])
                     for i in range(2) for k in ("fixed", "joint", "frozen"))
        wide_changes.append(change)
        refinements.append({"kind": "wide Hermite", "max_probability_change": change, "result": refined})
    check("resolution wide Hermite", max(wide_changes) < 1e-2)
    summaries = []
    for sigma in spread.SPREADS:
        selected = [r for r in rows if r["sigma_pointer"] == sigma]
        summary = {"sigma_pointer": sigma}
        for kind in ("fixed", "joint", "frozen"):
            values = np.array([r["contrasts"][kind] for r in selected])
            summary[kind] = {"contrast": float(values.mean()), "scramble_standard_error": float(values.std(ddof=1) / 2),
                             "min_contrast": float(values.min()), "max_contrast": float(values.max()),
                             "tails": [float(np.mean([r["choices"][i][kind]["tail"] for r in selected])) for i in range(2)]}
        summary["momentum_recoil_rms"] = float(np.mean([r["diagnostics"]["reference_momentum_recoil_rms"] for r in selected]))
        summaries.append(summary)
    return {"public_preregistration": spread.PUBLIC_PREREG, "model": asdict(spread.PointerModel()),
            "checks": checks, "passed": all(checks.values()), "summaries": summaries,
            "replicates": rows, "refinements": refinements, "mechanics": mechanics,
            "geometry": geometry, "preparation": preparation, "stationary_bridge": stationary,
            "independent_stationary_tails": baseline,
            "scope": "Specified coupled rotor history laws; no operational BAM field or intervention derivation."}


def render(report):
    lines = ["# Finite pointer-spread experiment", "", report["scope"], "",
             f"Public freeze: `{report['public_preregistration']}`.", "",
             "| sigma P | fixed-marginal contrast | scramble SE | joint-conditioned contrast | frozen-posterior contrast | momentum recoil RMS |",
             "|---:|---:|---:|---:|---:|---:|"]
    for row in report["summaries"]:
        lines.append(f"| {row['sigma_pointer']:.3g} | {row['fixed']['contrast']:.8f} | {row['fixed']['scramble_standard_error']:.2g} | {row['joint']['contrast']:.8f} | {row['frozen']['contrast']:.8f} | {row['momentum_recoil_rms']:.6g} |")
    lines += ["", "Scramble SE is a numerical diagnostic, not a rigorous confidence interval.", "",
              "| criterion | pass |", "|---|---|"]
    lines += [f"| {name} | {passed} |" for name, passed in report["checks"].items()]
    lines += ["", f"Passed {sum(report['checks'].values())}/{len(report['checks'])} criteria.", "",
              "The JSON retains all preparation-specific records, posterior P moments, path diagnostics and independent refinements.", ""]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    report = run_probe(progress=lambda s: print(s, flush=True))
    summary = render(report)
    if args.output_dir is not None:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "probe.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        (args.output_dir / "probe.md").write_text(summary)
    print(summary)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

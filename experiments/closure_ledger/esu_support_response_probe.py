"""Check the c07891d ESU support and cubic proper-time reaction freeze."""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from geometrodynamics.waves import reciprocal_scalar_tt as rt
from geometrodynamics.waves import esu_support_response as esu


def evolution_checks(model, times, spatial):
    reduced = model.integrate(times)
    fluid = model.integrate(times, independent_fluid=True)
    coarse = model.integrate(times, independent_fluid=True, rtol=1e-10, atol=1e-12)
    reconstructed = np.array([np.r_[z, model.fields(t, z)["rho_f"], model.fields(t, z)["J_f"]]
                              for t, z in zip(times, reduced)])
    ham, mom, curvature = [], [], []
    for t, z in zip(times, fluid):
        H, M = model.constraints(t, z)
        ham.append(H)
        mom.append(M)
        f = model.fields(t, z)
        curvature.append(f["curvature"]-f["trace_curvature"])
    # Unused constraints of the independent fluid solve, not identities made
    # zero by reconstructing rho and J from the constraints.
    alpha = np.array([model.fields(t, z)["alpha"] for t, z in zip(times, fluid)])
    psi_values = spatial.H @ fluid[:, :4].T
    alpha_values = spatial.H @ alpha.T
    bounds = model.continuous_potential_bounds(float(times[-1]))
    return {"sound_speed_squared": model.sound_speed_squared,
            "independent_fluid_scaled_error": esu.scaled_error(fluid, reconstructed),
            "ODE_refinement_scaled_error": esu.scaled_error(coarse, fluid),
            "unused_Hamiltonian_residual": float(np.max(np.abs(ham))),
            "unused_momentum_residual": float(np.max(np.abs(mom))),
            "curvature_trace_residual": float(np.max(np.abs(curvature))),
            "sampled_max_abs_psi": float(np.max(np.abs(psi_values))),
            "sampled_max_abs_alpha": float(np.max(np.abs(alpha_values))),
            "continuous_all_space_potential_bounds": bounds,
            "initial_proper_coefficient": model.initial_proper_force()/model.cubic_scale,
            "initial_coordinate_coefficient": model.coordinate_force(0., model.initial())/model.cubic_scale,
            "times": times.tolist(), "psi_group_amplitudes": fluid[:, :4].tolist(),
            "support_density_group_amplitudes": fluid[:, 8:12].tolist(),
            "coordinate_scalar_force_over_cubic_scale": [model.coordinate_force(t, z)/model.cubic_scale for t, z in zip(times, fluid)],
            "TT_force_over_cubic_scale": [model.induced_tt_force(t)/model.cubic_scale for t in times]}, reduced


def initial_force_checks(model, spatial):
    forces = spatial.force_fields(0., model.initial())
    projected = {k: spatial.project_force(v)/model.cubic_scale for k, v in forces.items()}
    return {"proper_quadrature_coefficient": projected["initial_proper_geometry"],
            "coordinate_quadrature_coefficient": projected["coordinate"],
            "proper_conversion_coefficient": projected["initial_proper_conversion"],
            "proper_coefficient_error": abs(projected["initial_proper_geometry"]+7976/875),
            "coordinate_coefficient_error": abs(projected["coordinate"]-55096/875),
            "clock_conversion_scaled_error": esu.scaled_error(
                forces["initial_proper_conversion"]/model.cubic_scale,
                forces["initial_proper_geometry"]/model.cubic_scale)}


def metric_variation_checks(model, spatial, times, states):
    rows = []
    for index in (0, 37, 100, 200):
        t, z = times[index], states[index]
        target = spatial.force_fields(t, z)["coordinate"]
        derivatives = []
        for epsilon in (.04, .02, .01):
            plus = spatial.exact_coordinate_acceleration(t, z, epsilon)
            minus = spatial.exact_coordinate_acceleration(t, z, -epsilon)
            derivatives.append((plus-minus)/(2*epsilon))
        richardson = (4*derivatives[-1]-derivatives[-2])/3
        rows.append({"time": float(t),
                     "centered_scaled_errors": [esu.scaled_error(d/model.cubic_scale, target/model.cubic_scale) for d in derivatives],
                     "Richardson_scaled_error": esu.scaled_error(richardson/model.cubic_scale, target/model.cubic_scale)})
    return rows


def induced_tt_check(model, times, direction):
    rm = rt.ReciprocalModel(radius=model.radius, kappa=model.kappa)
    initial = rm.pack(np.zeros(5), np.zeros(5), model.amplitude*direction, np.zeros(16))
    # Independent old scalar/TT implementation in its prescribed-source mode.
    states = rm.integrate(initial, times, reciprocal=False, rtol=1e-12, atol=1e-14)
    reference = []
    for t, z in zip(times, states):
        b, _, q, _ = rm.unpack(z)
        reference.append(float(direction @ (2*np.einsum("a,aij,j->i", b, rm.F, q))))
    predicted = np.array([model.induced_tt_force(t) for t in times])
    return esu.scaled_error(predicted/model.cubic_scale, np.array(reference)/model.cubic_scale)


def run_probe(progress=lambda message: None):
    progress("metric/connection derivation and exact initial coefficients")
    geometry = esu.metric_einstein_derivation()
    exact = esu.initial_coefficient_certificate()
    model = esu.SupportResponse()
    times = np.linspace(0., 2., 201)
    progress("improved stress, anisotropic scalar projection and proper-time coefficient")
    grids = [esu.SpatialFields(model, *orders) for orders in ((8, 16), (12, 24))]
    stress = [{"points": len(grid.grid.points), "time": t, **grid.stress_check(t)}
              for grid in grids for t in (0., math.pi/16, .37, 1.)]
    initial = [initial_force_checks(model, grid) for grid in grids]
    spatial_refinement = max(abs(initial[0][key]-initial[1][key]) for key in (
        "proper_quadrature_coefficient", "coordinate_quadrature_coefficient"))
    progress("independent support conservation and unused Einstein constraints")
    controls = []
    primary_states = None
    for cs2 in (1/3, 0., .2, 1.):
        candidate = esu.SupportResponse(sound_speed_squared=cs2)
        result, states = evolution_checks(candidate, times, grids[0])
        controls.append(result)
        if cs2 == 1/3:
            primary_states = states
    progress("nonlinear metric control, amplitude/radius scaling and TT comparison")
    variation = metric_variation_checks(model, grids[0], times, primary_states)
    amplitudes = []
    for amplitude in (.02, .01, .005):
        m = esu.SupportResponse(amplitude=amplitude)
        amplitudes.append({"amplitude": amplitude,
                           "initial_proper_force": m.initial_proper_force(),
                           "initial_metric_group_norm": float(np.linalg.norm(m.initial()[:4]))})
    force_scaling = max(abs(amplitudes[i]["initial_proper_force"]/amplitudes[i+1]["initial_proper_force"]-8) for i in (0, 1))
    metric_scaling = max(abs(amplitudes[i]["initial_metric_group_norm"]/amplitudes[i+1]["initial_metric_group_norm"]-4) for i in (0, 1))
    radii = []
    for radius in (.7, 2.):
        m = esu.SupportResponse(radius=radius, kappa=.4)
        fields = esu.SpatialFields(m)
        radii.append({"radius": radius, "kappa": m.kappa, **initial_force_checks(m, fields)})
    tt_error = induced_tt_check(model, times, grids[0].wave.direction)
    checks = {
        "metric_Einstein_derivation": geometry["all_exact_zero"],
        "improved_stress_and_scalar_anisotropy": max(r[k] for r in stress for k in (
            "rho", "pressure", "momentum", "anisotropic_scalar_projection")) < 1e-9,
        "two_spatial_rules": spatial_refinement < 1e-9,
        "independent_fluid_evolution": max(r["independent_fluid_scaled_error"] for r in controls) < 1e-8,
        "propagated_constraints": max(max(r["unused_Hamiltonian_residual"], r["unused_momentum_residual"]) for r in controls) < 1e-8,
        "ODE_refinement": max(r["ODE_refinement_scaled_error"] for r in controls) < 1e-8,
        "trace_curvature": max(r["curvature_trace_residual"] for r in controls) < 1e-9,
        "fluid_only_frequency": geometry["ell2_threshold"] == "1/5" and geometry["ell0_frequency_times_a2"] == "-3*cs2 - 1",
        "nonlinear_metric_variation": max(r["Richardson_scaled_error"] for r in variation) < 1e-7,
        "initial_proper_and_coordinate_coefficients": max(max(r["proper_coefficient_error"], r["coordinate_coefficient_error"]) for r in initial+radii) < 1e-9,
        "clock_conversion": max(r["clock_conversion_scaled_error"] for r in initial+radii) < 1e-9,
        "amplitude_and_radius_scaling": force_scaling < 1e-8 and metric_scaling < 1e-8,
        "induced_TT_force": tt_error < 1e-8 and model.induced_tt_force(0.) == 0.,
        "small_metric_regime": max(max(r["continuous_all_space_potential_bounds"].values()) for r in controls) < .05
             and all(r["sampled_max_abs_psi"] <= r["continuous_all_space_potential_bounds"]["psi"]
                     and r["sampled_max_abs_alpha"] <= r["continuous_all_space_potential_bounds"]["alpha"] for r in controls),
        "exact_preparation_obstruction": exact["proper_expected"] and exact["coordinate_expected"]
             and exact["proper_total"] != "0" and all(abs(r["initial_proper_coefficient"]+7976/875) < 1e-9 for r in controls),
    }
    checks = {key: bool(value) for key, value in checks.items()}
    return {"baseline": esu.BASELINE, "public_preregistration": esu.PUBLIC_PREREG, "seed": esu.SEED,
            "primary_model": {"radius": model.radius, "kappa": model.kappa, "amplitude": model.amplitude,
                              "sound_speed_squared": model.sound_speed_squared, "time_interval": [0., 2.],
                              "coefficient_unit": "kappa s^3/(V a^2)"},
            "scope": "parameterized perfect-fluid scalar response, with initial zero support perturbations and no free TT data",
            "metric_derivation": geometry, "exact_initial_certificate": exact,
            "stress_checks": stress, "initial_force_checks": initial,
            "spatial_refinement_error": spatial_refinement, "support_controls": controls,
            "metric_variation_checks": variation, "amplitude_checks": amplitudes,
            "force_scaling_error": force_scaling, "metric_scaling_error": metric_scaling,
            "radius_checks": radii, "induced_TT_scaled_error": tt_error,
            "checks": checks, "checks_passed": not esu.failed_checks(checks), "verdict": esu.verdict(checks)}


def render(report):
    lines = ["# ESU support response and first cubic reaction", "",
             "The support class is assumed; its sound speed is a parameter. The proper-time coefficient is preparation-scoped.",
             "The continuous metric bounds apply to t in [0,2], with a=kappa=1 and s=0.02.", ""]
    lines += [f"- {key}: **{value}**" for key, value in report["verdict"].items()]
    exact = report["exact_initial_certificate"]
    lines += ["", "| Initial projection, units kappa s^3/(V a^2) | Exact coefficient |", "|---|---:|",
              f"| Fluid proper-time scalar force | {exact['proper_total']} |",
              f"| Newtonian coordinate-time scalar force | {exact['coordinate_total']} |",
              "| Induced homogeneous TT force | 0 |", "",
              "| c_s^2 | Unused Hamiltonian residual | Unused momentum residual | All-space/all-time psi upper bound |",
              "|---:|---:|---:|---:|"]
    for r in report["support_controls"]:
        lines.append(f"| {r['sound_speed_squared']:.6g} | {r['unused_Hamiltonian_residual']:.3g} | "
                     f"{r['unused_momentum_residual']:.3g} | {r['continuous_all_space_potential_bounds']['psi']:.6g} |")
    lines += ["", "| Required check | Pass |", "|---|---|"]
    lines += [f"| {key} | {ok} |" for key, ok in report["checks"].items()]
    lines += ["", f"Passed {sum(report['checks'].values())}/{len(esu.REQUIRED_CHECKS)} required checks.", "",
              "No full corrected scalar evolution, microscopic BAM support, or source-local readout is derived."]
    return "\n".join(lines)+"\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    report = run_probe(progress=lambda message: print(message, flush=True))
    summary = render(report)
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir/"probe.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
        (args.output_dir/"probe.md").write_text(summary)
    print(summary)
    return 0 if report["checks_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

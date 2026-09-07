"""Verify the e1706b1 conditional constraint-response bounds on PR #289."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from geometrodynamics.bulk.tt_triangle_rotor import nearest_uniaxial
from geometrodynamics.waves import reciprocal_scalar_tt as rt
from geometrodynamics.waves import scalar_tt_constraints as sc


def scaled_error(actual, expected):
    """Euclidean error divided by max(1, Euclidean norm of reference)."""
    return float(np.linalg.norm(actual-expected)/max(1., np.linalg.norm(expected)))


def physical_rms(grid, values):
    fields = np.asarray(values).reshape(len(grid.points), -1)
    return float(np.sqrt(grid.weights @ np.sum(fields*fields, axis=1)/grid.volume))


def grid_checks(grid, wave):
    checks = {"square_relative_L2_reconstruction": physical_rms(grid, grid.evaluate(wave.h)-wave.mode**2)
              /physical_rms(grid, wave.mode**2)}
    rows = []
    for t in (0., math.pi/16, .3, 1., 2., 4.):
        r = wave.response(t)
        rho, j = sc.round_sources(wave.model, r["A"]*wave.direction,
                                  r["dA"]*wave.direction, grid.points)
        rho_fit, j_fit = grid.evaluate(r["rho"]), grid.gradient(r["J"])
        Acoeff = grid.longitudinal_coefficients(r["w"])
        Afield = np.einsum("pn,ijn->pij", grid.Y, Acoeff)
        div = grid.Y @ grid.divergence_coefficients(Acoeff).T
        ham = grid.evaluate(grid.laplacian_coefficients(r["u"])+3*r["u"])
        row = {"time": t, "rho_formula_scaled": scaled_error(rho_fit, rho),
               "j_formula_scaled": scaled_error(j_fit, j),
               "Hamiltonian_residual_scaled": scaled_error(ham, -wave.model.kappa*rho/4),
               "momentum_residual_scaled": scaled_error(div, wave.model.kappa*j),
               "K_norm_identity_error": abs(physical_rms(grid, Afield)-wave.norms(t)["K_longitudinal_rms"]),
               "momentum_compatibility": float(np.max(np.abs(grid.require_compatible_momentum(j)))),
               "energy_dipole_overlap": float(np.max(np.abs(np.einsum("p,p,pa->a", grid.weights, rho, grid.points))))}
        rows.append(row)
    checks["times"] = rows
    return checks


def arbitrary_energy_checks(grid):
    rng = np.random.default_rng(sc.SEED)
    model = rt.ReciprocalModel()
    cases = []
    for _ in range(3):
        q, p = rng.normal(size=(2, 16))
        cases.append((.2*q/np.linalg.norm(q), .3*p/np.linalg.norm(p)))
    direction = model.multiplet.coherent(np.array([2., -3., 1.])/math.sqrt(14))
    cases.append((.2*direction, -.4*direction))
    rows = []
    for q, p in cases:
        rho, _ = sc.round_sources(model, q, p, grid.points)
        coeff = grid.project(rho)
        u = sc.solve_hamiltonian(coeff, grid.degrees)
        lap_u = grid.evaluate(grid.laplacian_coefficients(u)+3*u)
        rows.append({"energy_relative_L2_reconstruction": physical_rms(grid, grid.evaluate(coeff)-rho)/physical_rms(grid, rho),
                     "Hamiltonian_residual_scaled": scaled_error(lap_u, -rho/4),
                     "energy_dipole_overlap": float(np.max(np.abs(np.einsum("p,p,pa->a", grid.weights, rho, grid.points)))),
                     "upper_bound_margin": grid.rms(coeff)/12-grid.rms(u)})
    return rows


def compatibility_controls(grid, wave):
    try:
        sc.solve_hamiltonian([1.], [1])
    except ValueError:
        dipole_rejected = True
    else:
        dipole_rejected = False
    q = .2*wave.direction
    p = wave.model.multiplet.D[0] @ q
    _, j = sc.round_sources(wave.model, q, p, grid.points)
    try:
        grid.require_compatible_momentum(j)
    except ValueError:
        charged_rejected = True
    else:
        charged_rejected = False
    charge = float(grid.weights @ j[:, 0])
    predicted = float(-p @ wave.model.multiplet.D[0] @ q)
    return {"dipole_rejected": dipole_rejected, "charged_momentum_rejected": charged_rejected,
            "left_Killing_charge": charge, "predicted_charge": predicted,
            "charge_error": abs(charge-predicted)}


def bounds_hold(row, tol=1e-10):
    return bool(row["u_rms"] <= row["u_rms_upper"]+tol
                and row["u_inhomogeneous_lower"]-tol <= row["u_inhomogeneous_rms"]
                <= row["u_inhomogeneous_upper"]+tol
                and row["K_longitudinal_lower"]-tol <= row["K_longitudinal_rms"]
                <= row["K_longitudinal_upper"]+tol)


def sharp_bounds(grid):
    zero = np.zeros(84)
    errors = []
    for index, divisor in ((0, 12), (1, 20)):
        rho = zero.copy()
        rho[index] = 1.
        u = sc.solve_hamiltonian(rho, grid.degrees)
        errors.append(abs(np.linalg.norm(u)-1/divisor))
    J = zero.copy()
    J[1] = 1.
    w = zero.copy()
    w[1] = -3/20
    A = grid.longitudinal_coefficients(w)
    errors.append(abs(np.linalg.norm(A)-math.sqrt(3/10)*math.sqrt(8)))
    return errors


def source_shape(model, history, times):
    rows = []
    for t, state in zip(times, history):
        _, _, q, _ = model.unpack(state)
        S = rt.tensor(model.source(q))
        norm = float(np.linalg.norm(S))
        distance = nearest_uniaxial(S)["distance"]
        rows.append({"time": float(t), "source_norm": norm, "uniaxial_distance": distance,
                     "relative_distance": distance/norm if norm else None})
    return {"initial": rows[0], "at_time_2": rows[200],
            "minimum_norm": min(r["source_norm"] for r in rows),
            "maximum_norm": max(r["source_norm"] for r in rows),
            "maximum_relative_distance": max(r["relative_distance"] for r in rows if r["relative_distance"] is not None),
            "scope": "old reciprocal projected history; no constraint evolution inferred"}


def run_probe(progress=lambda message: None):
    progress("even harmonic reconstruction and both pointwise constraints")
    grids = [sc.EvenHarmonicGrid(*order) for order in ((8, 16), (12, 24))]
    waves = [sc.StandingWaveConstraints(grid) for grid in grids]
    results = [grid_checks(grid, wave) for grid, wave in zip(grids, waves)]
    grid, wave = grids[0], waves[0]
    arbitrary = arbitrary_energy_checks(grid)
    negative = compatibility_controls(grid, wave)
    sharp = sharp_bounds(grid)
    times = np.linspace(0., 4., 401)
    history = [{"time": float(t), **wave.norms(t)} for t in times]
    envelope = wave.envelopes()
    fine_envelope = waves[1].envelopes()
    refinement = max(abs(envelope[k][n]-fine_envelope[k][n])
                     for k in ("u", "u_inhomogeneous") for n in ("minimum_rms", "maximum_rms"))
    refinement = max(refinement, scaled_error(wave.h, waves[1].h),
                     abs(envelope["K_longitudinal_maximum_rms"]-fine_envelope["K_longitudinal_maximum_rms"]))
    progress("continuous bounds, amplitude scaling and independent tensor response")
    scaling = []
    for amplitude in (.2, .1, .05):
        candidate = sc.StandingWaveConstraints(grid, amplitude=amplitude)
        r = candidate.norms(math.pi/16)
        scaling.append({"amplitude": amplitude, **{k: r[k] for k in (
            "u_rms", "u_inhomogeneous_rms", "K_longitudinal_rms", "induced_TT_metric_frobenius")}})
    scaling_error = max(abs(scaling[i][key]/scaling[i+1][key]-4)
                        for i in (0, 1) for key in scaling[i] if key != "amplitude")
    model = wave.model
    omega = math.sqrt(model.omega_scalar2)
    f0 = model.source(.2*wave.direction)/model.C
    def forced(t, y):
        return np.r_[y[5:], -model.omega_tensor2*y[:5]+f0*math.cos(omega*t)**2]
    sol = solve_ivp(forced, (0., 4.), np.zeros(10), t_eval=times,
                    method="DOP853", rtol=1e-12, atol=1e-14)
    if not sol.success:
        raise RuntimeError(sol.message)
    tensor_error = float(np.max(np.abs(sol.y[:5].T-np.array([wave.induced_tt(t) for t in times]))))
    progress("computed TT constraint certificates and source-shape review follow-up")
    certificate = sc.homogeneous_tt_constraint_certificate()
    square_certificate = sc.coherent_square_certificate()
    measured_powers = [float(grid.volume*np.sum(wave.h[grid.degrees == l]**2)) for l in (0, 2, 4, 6)]
    square_power_error = scaled_error(np.array(measured_powers), np.array(square_certificate["numeric_powers_times_volume"]))
    primary = model.integrate(rt.primary_data(model), times, rtol=1e-12, atol=1e-14)
    primary_tt = 2*np.linalg.norm(primary[:, :5], axis=1)
    shape = source_shape(model, primary, times)
    h1 = rt.harmonic_multiplet(1)
    anticommutator = max(np.linalg.norm(h1.D[i] @ h1.D[j]+h1.D[j] @ h1.D[i]
                                       + 2*(i == j)*np.eye(4)) for i in range(3) for j in range(3))
    per_time = [r for result in results for r in result["times"]]
    envelopes_hold = all(
        envelope["u"]["minimum_rms"]-1e-10 <= r["u_rms"] <= envelope["u"]["maximum_rms"]+1e-10
        and envelope["u_inhomogeneous"]["minimum_rms"]-1e-10 <= r["u_inhomogeneous_rms"]
        <= envelope["u_inhomogeneous"]["maximum_rms"]+1e-10
        and r["K_longitudinal_rms"] <= envelope["K_longitudinal_maximum_rms"]+1e-10
        and r["induced_TT_metric_frobenius"] <= envelope["induced_TT_metric_all_time_upper"]+1e-10
        and r["scalar_metric_inhomogeneous_rms"] >= envelope["scalar_metric_inhomogeneous_all_time_lower"]-1e-10
        for r in history)
    checks = {
        "full_even_reconstruction": max(r["square_relative_L2_reconstruction"] for r in results) < 1e-9
             and max(r["energy_relative_L2_reconstruction"] for r in arbitrary) < 1e-9,
        "independent_improved_stress": max(max(r["rho_formula_scaled"], r["j_formula_scaled"]) for r in per_time) < 1e-10,
        "Hamiltonian_and_momentum_residuals": max(max(r["Hamiltonian_residual_scaled"], r["momentum_residual_scaled"]) for r in per_time) < 1e-9
             and max(r["Hamiltonian_residual_scaled"] for r in arbitrary) < 1e-9,
        "dipole_and_CKV_compatibility": max(max(r["momentum_compatibility"], r["energy_dipole_overlap"]) for r in per_time) < 1e-10
             and max(r["energy_dipole_overlap"] for r in arbitrary) < 1e-10,
        "K_norm_identity": max(r["K_norm_identity_error"] for r in per_time) < 1e-9,
        "two_quadrature_rules": refinement < 1e-9,
        "spectral_bounds_and_sharp_constants": all(bounds_hold(r) for r in history)
             and max(sharp) < 1e-10 and min(r["upper_bound_margin"] for r in arbitrary) > -1e-10,
        "incompatible_sources_rejected": negative["dipole_rejected"] and negative["charged_momentum_rejected"]
             and abs(negative["left_Killing_charge"]) > 1e-3 and negative["charge_error"] < 1e-10,
        "quadratic_amplitude_scaling": scaling_error < 1e-8,
        "independent_forced_TT_solution": tensor_error < 1e-9,
        "continuous_envelopes": envelopes_hold,
        "computed_homogeneous_TT_zeros": certificate["exact"],
        "degree_one_anticommutator": anticommutator < 1e-12,
        "exact_coherent_power_and_size_certificate": square_certificate["exact"] and square_power_error < 1e-10,
    }
    checks = {key: bool(value) for key, value in checks.items()}
    return {"baseline": sc.BASELINE, "public_preregistration": sc.PUBLIC_PREREG,
            "seed": sc.SEED, "norm_convention": "RMS over physical S3; tensor Frobenius in background orthonormal frame",
            "residual_convention": "scaled means Euclidean error/max(1, reference Euclidean norm)",
            "model": {"radius": 1., "kappa": 1., "amplitude": .2, "degrees": [0, 2, 4, 6],
                      "support_delta_rho_and_delta_j": "chosen zero on independent CMC slices",
                      "evolution_completion": "unspecified"},
            "grid_checks": results, "arbitrary_energy_checks": arbitrary,
            "compatibility_controls": negative, "sharp_bound_errors": sharp,
            "refinement_scaled_error": refinement, "amplitude_scaling": scaling,
            "amplitude_scaling_error": scaling_error, "forced_TT_error": tensor_error,
            "envelopes": envelope, "square_harmonic_coefficients": wave.h.tolist(),
            "square_harmonic_degrees": grid.degrees.tolist(),
            "homogeneous_TT_certificate": certificate, "degree_one_anticommutator_error": float(anticommutator),
            "coherent_square_certificate": square_certificate, "square_power_error": square_power_error,
            "primary_TT_metric_sampled_range": [float(min(primary_tt)), float(max(primary_tt))],
            "primary_source_shape": shape,
            "history": history, "checks": checks, "checks_passed": all(checks.values()),
            "verdict": sc.verdict(checks)}


def render(report):
    lines = ["# Conditional scalar–TT constraint response", "",
             "Each time is an independent linear CMC slice with zero support perturbations.",
             "These are metric bounds; the complete scalar backreaction is not bounded here.", ""]
    lines += [f"- {k}: **{v}**" for k, v in report["verdict"].items()]
    env = report["envelopes"]
    lines += ["", "| Continuous quantity | Value |", "|---|---:|"]
    for key in ("u", "u_inhomogeneous"):
        for extreme in ("minimum_rms", "maximum_rms"):
            lines.append(f"| {key} {extreme} | {env[key][extreme]:.12g} |")
    for key in ("K_longitudinal_maximum_rms", "j_maximum_rms", "scalar_metric_inhomogeneous_all_time_lower", "induced_TT_metric_all_time_upper"):
        lines.append(f"| {key} | {env[key]:.12g} |")
    lines += ["", "| Time | u mean | u RMS without mean | K longitudinal RMS | scalar metric RMS without mean | induced TT metric norm |",
              "|---:|---:|---:|---:|---:|---:|"]
    for r in report["history"][::100]:
        lines.append("| "+" | ".join(f"{r[k]:.10g}" for k in (
            "time", "u_mean", "u_inhomogeneous_rms", "K_longitudinal_rms",
            "scalar_metric_inhomogeneous_rms", "induced_TT_metric_frobenius"))+" |")
    lines += ["", "| Required check | Pass |", "|---|---|"]
    lines += [f"| {key} | {value} |" for key, value in report["checks"].items()]
    lines += ["", f"Passed {sum(report['checks'].values())}/{len(report['checks'])} checks.",
              "", "A nonzero constraint response establishes neither an operational readout nor a completed Einstein history."]
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

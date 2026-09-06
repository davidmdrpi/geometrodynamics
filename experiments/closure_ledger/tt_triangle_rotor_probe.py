"""Preregistered TT-to-rotor test; archive evidence and fail on failed checks."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm, expm_frechet

from geometrodynamics.bulk import tt_triangle_rotor as rotor
from geometrodynamics.waves import backreaction


def ode_flow(beta0, velocity0, times, source=None, model=rotor.TensorModel()):
    """Independent unreduced matrix ODE: no director extraction/projection."""
    times = np.asarray(times)

    def rhs(t, y):
        forcing = np.zeros((3, 3)) if source is None else source(t)
        return np.r_[y[9:], np.asarray(forcing).ravel() - model.omega2 * y[:9]]

    result = solve_ivp(rhs, (0., float(times.max())),
                      np.r_[beta0.ravel(), velocity0.ravel()], t_eval=times,
                      method="DOP853", rtol=1e-11, atol=1e-13)
    if not result.success:
        raise RuntimeError(result.message)
    return result.y[:9].T.reshape(-1, 3, 3), result.y[9:].T.reshape(-1, 3, 3)


def algebra_controls():
    rng = np.random.default_rng(rotor.SEED)
    model = rotor.TensorModel()
    worst = {k: 0. for k in ("kinetic", "lagrangian", "normal_formula",
        "normal_basis", "normal_norm", "radial_equation", "angular_equation",
        "covariance", "antipodal", "metric_finite_difference", "adm_kinetic")}

    def record(name, value):
        worst[name] = max(worst[name], float(value))

    for _ in range(32):
        R, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        R[:, 0] *= np.linalg.det(R)
        e, f, n = R.T
        v = .4 * (rng.normal() * e + rng.normal() * f)
        amplitude = rng.choice([-1., 1.]) * rng.uniform(.01, .03)
        adot = .01 * rng.normal()
        source = .03 * rotor.stf(rng.normal(size=(3, 3)))
        B, Bd = rotor.embedding(amplitude, n), rotor.field_velocity(amplitude, adot, n, v)
        record("kinetic", abs(np.sum(Bd * Bd) - (2 * adot ** 2 / 3 + 2 * amplitude ** 2 * (v @ v))))
        record("lagrangian", abs(rotor.lagrangian(B, Bd, source)
                                  - rotor.restricted_lagrangian(amplitude, adot, n, v, source)))
        addot, nddot = rotor.restricted_acceleration(amplitude, adot, n, v, source)
        Bdd = rotor.field_acceleration(amplitude, adot, addot, n, v, nddot)
        residual = Bdd + model.omega2 * B - source
        parts = rotor.tensor_components(residual, n)
        predicted = rotor.omitted_equation(amplitude, n, v, source)
        record("normal_formula", np.linalg.norm(residual - predicted))
        record("radial_equation", np.linalg.norm(parts["radial"]))
        record("angular_equation", np.linalg.norm(parts["angular"]))
        basis = [(np.outer(e, e) - np.outer(f, f)) / math.sqrt(2),
                 (np.outer(e, f) + np.outer(f, e)) / math.sqrt(2)]
        explicit = sum(np.sum(source * N) * N for N in basis)
        record("normal_basis", np.linalg.norm(explicit - rotor.normal_projection(source, n)))
        record("normal_norm", abs(np.linalg.norm(rotor.omitted_equation(amplitude, n, v))
                                   - math.sqrt(2) * abs(amplitude) * (v @ v)))
        record("covariance", np.linalg.norm(rotor.omitted_equation(amplitude, R @ n, R @ v, R @ source @ R.T)
                                              - R @ predicted @ R.T))
        record("antipodal", max(np.linalg.norm(rotor.embedding(amplitude, -n) - B),
                                 np.linalg.norm(rotor.field_velocity(amplitude, adot, -n, -v) - Bd)))

        def chart(z):
            A, theta, phi = z
            axis = np.array([math.sin(theta) * math.cos(phi),
                             math.sin(theta) * math.sin(phi), math.cos(theta)])
            return rotor.embedding(A, axis).ravel()

        point = np.array([amplitude, rng.uniform(.3, 2.8), rng.uniform(-2, 2)])
        h = 1e-6
        derivative = np.column_stack([(chart(point + h * d) - chart(point - h * d)) / (2 * h)
                                      for d in np.eye(3)])
        metric = model.normalization * derivative.T @ derivative
        expected_metric = rotor.pullback_metric(amplitude, point[1])
        # Relative in each orthonormalized coordinate, so small A cannot hide an error.
        scaling = np.sqrt(np.outer(np.diag(expected_metric), np.diag(expected_metric)))
        record("metric_finite_difference", np.max(abs(metric - expected_metric) / scaling))

        # Noncommuting matrix velocity: ADM K=(1/2)g^{-1}gdot, not a diagonal-only test.
        eps = 1e-3
        g = expm(2 * eps * B)
        gd = expm_frechet(2 * eps * B, 2 * eps * Bd, compute_expm=False)
        K = .5 * np.linalg.solve(g, gd)
        kinetic = (np.trace(K @ K) - np.trace(K) ** 2) / eps ** 2
        record("adm_kinetic", abs(kinetic - np.sum(Bd ** 2)))
    return worst


def free_motion(amplitude=.01, speed=.4):
    n, v = np.array([0., 0., 1.]), np.array([speed, 0., 0.])
    B0, Bd0 = rotor.embedding(amplitude, n), rotor.field_velocity(amplitude, 0., n, v)
    times = np.array([.04, .02, .01, .005])
    B, Bd = rotor.exact_free_flow(B0, Bd0, times)
    nearest = [rotor.nearest_uniaxial(b) for b in B]
    distances = [r["distance"] for r in nearest]
    coefficient = abs(amplitude) * speed ** 2 / math.sqrt(2)
    return {"amplitude": amplitude, "speed": speed,
            "times": times.tolist(), "distances": distances,
            "distances_over_t2": (np.array(distances) / times ** 2).tolist(),
            "predicted_coefficient": coefficient,
            "halving_ratios": [distances[i] / distances[i+1] for i in range(3)] if speed else [],
            "eigenvalue_distance_error": max(abs(r["distance"] - r["eigenvalue_distance"]) for r in nearest),
            "free_normal_residual": float(np.linalg.norm(rotor.omitted_equation(amplitude, n, v))),
            "predicted_normal_residual": math.sqrt(2) * abs(amplitude) * speed ** 2,
            "energy_residual": float(np.max(abs(rotor.free_energy(B, Bd) - rotor.free_energy(B0, Bd0)))),
            "charge_residual": float(np.max(abs(rotor.rotational_charge(B, Bd) - rotor.rotational_charge(B0, Bd0))))}


def dynamics_controls():
    n, v = np.array([0., 0., 1.]), np.array([.4, 0., 0.])
    B0, Bd0 = rotor.embedding(.01, n), rotor.field_velocity(.01, 0., n, v)
    times = np.linspace(0, .3, 121)
    B, Bd = rotor.exact_free_flow(B0, Bd0, times)
    independent, independent_dot = ode_flow(B0, Bd0, times)
    ode_error = max(float(np.max(abs(B - independent))), float(np.max(abs(Bd - independent_dot))))
    field = rotor.uniform_rotation(0.)
    manufactured, _ = ode_flow(field["beta"], field["beta_dot"], times,
                               source=lambda t: rotor.uniform_rotation(t)["source"])
    fields = [rotor.uniform_rotation(t) for t in times]
    manufactured_error = float(np.max(abs(manufactured - np.array([f["beta"] for f in fields]))))
    radial, normal, equation = [], [], []
    for f in fields:
        parts = rotor.tensor_components(f["source"], f["n"])
        radial.append(float(np.linalg.norm(parts["radial"])))
        normal.append(float(np.linalg.norm(parts["normal"])))
        addot, nddot = rotor.restricted_acceleration(.01, 0., f["n"], f["v"], f["source"])
        full = rotor.field_acceleration(.01, 0., addot, f["n"], f["v"], nddot)
        equation.append(float(np.linalg.norm(full + 8 * f["beta"] - f["source"])))
    # A constant-amplitude stationary tensor also requires radial support.
    fixed_amplitude_source = 8 * rotor.embedding(.01, n)
    return {"free_ode_error": ode_error, "manufactured_ode_error": manufactured_error,
            "manufactured_full_equation_residual": max(equation),
            "manufactured_min_radial_force": min(radial), "manufactured_min_normal_force": min(normal),
            "frozen_amplitude_stationary_radial_force": float(np.linalg.norm(fixed_amplitude_source)),
            "max_tensor_norm": float(np.max(np.linalg.norm(B, axis=(1, 2)))),
            "late_distance_t_0_3": rotor.nearest_uniaxial(B[-1])["distance"]}


def inherited_response_control():
    """Check normalization and the actual repository response kernel against ODE."""
    model = rotor.TensorModel()
    integrated = rotor.stf(np.array([[.02, .01, 0.], [.01, -.015, .005], [0., .005, -.005]]))
    S = rotor.normalized_source(integrated)
    errors = []
    for count in (513, 1025):
        times = np.linspace(0, .4, count)
        envelope = np.sin(math.pi * times / .4) ** 2
        field = backreaction.shear_response(envelope[:, None, None] * S, times[1])
        reference, _ = ode_flow(np.zeros((3, 3)), np.zeros((3, 3)), times,
                               source=lambda t: math.sin(math.pi * t / .4) ** 2 * S)
        errors.append(float(np.max(abs(field - reference))))
    return {"errors_513_1025": errors,
            "integrated_source_factor": 1 / model.normalization,
            "source_variation_identity_residual": float(np.linalg.norm(model.normalization * S - integrated))}


def post_review_controls():
    """User-review additions, not part of the original frozen 22 checks."""
    omega = math.sqrt(8)
    times = np.linspace(0, 2 * math.pi / omega, 257)
    spectrum_error = distance_error = return_error = eigenline_error = 0.
    for speed in (.2, .4, .8):
        for amplitude in (-.01, .01):
            B0 = rotor.embedding(amplitude, [0., 0., 1.])
            Bd0 = rotor.field_velocity(amplitude, 0., [0., 0., 1.], [speed, 0., 0.])
            B, _ = rotor.exact_free_flow(B0, Bd0, times)
            closed = rotor.closed_form_free_spectrum(times, amplitude, speed)
            spectrum_error = max(spectrum_error, float(np.max(abs(np.linalg.eigvalsh(B) - closed["eigenvalues"]))))
            distances = np.array([rotor.nearest_uniaxial(b)["distance"] for b in B])
            distance_error = max(distance_error, float(np.max(abs(distances - closed["distance"]))))
            returns, _ = rotor.exact_free_flow(B0, Bd0, np.arange(4) * math.pi / omega)
            return_error = max(return_error, max(rotor.nearest_uniaxial(b)["distance"] for b in returns))
        B0 = rotor.embedding(.01, [0., 0., 1.])
        Bd0 = rotor.field_velocity(.01, 0., [0., 0., 1.], [speed, 0., 0.])
        B, _ = rotor.exact_free_flow(B0, Bd0, times)
        line = rotor.continuous_eigenline(times, speed)
        axes = line["axis"]
        eigenvalues = np.einsum("ni,nij,nj->n", axes, B, axes)
        eigenline_error = max(eigenline_error, float(np.max(abs(np.einsum("nij,nj->ni", B, axes) - eigenvalues[:, None] * axes))))

    fractions = np.array([0., .25, .375, .5, .625, .75, 1.])
    line = rotor.continuous_eigenline(fractions * math.pi / omega)
    full_line = rotor.continuous_eigenline(times)
    quarter = math.pi / (2 * omega)
    B0 = rotor.embedding(.01, [0., 0., 1.])
    Bd0 = rotor.field_velocity(.01, 0., [0., 0., 1.], [.4, 0., 0.])
    midpoint, _ = rotor.exact_free_flow(B0, Bd0, quarter)
    a1, a2 = np.array([1., 0., 1.]) / math.sqrt(2), np.array([-1., 0., 1.]) / math.sqrt(2)
    fits = [rotor.embedding(1.5 * (a @ midpoint @ a), a) for a in (a1, a2)]
    minimum = rotor.nearest_uniaxial(midpoint)
    fit_errors = [abs(np.linalg.norm(midpoint - fit) - minimum["distance"]) for fit in fits]
    before, _ = rotor.exact_free_flow(B0, Bd0, quarter - 1e-7)
    after, _ = rotor.exact_free_flow(B0, Bd0, quarter + 1e-7)
    switch = rotor.nearest_uniaxial(before)["amplitude"] * rotor.nearest_uniaxial(after)["amplitude"] < 0
    checks = {
        "R1 exact full-period spectrum and distance": max(spectrum_error, distance_error) < 1e-12,
        "R2 periodic returns do not imply interval invariance": return_error < 1e-12 and minimum["distance"] > 1e-4,
        "R3 a continuously advancing eigenline exists": eigenline_error < 1e-12
            and np.min(full_line["angular_speed"]) > 0 and np.min(np.diff(full_line["angle"])) > 0,
        "R4 nearest-uniaxial axis switches between distinct optima": max(fit_errors) < 1e-12
            and not minimum["axis_defined"] and abs(a1 @ a2) < 1e-12 and switch,
    }
    return {"status": "POST_REVIEW_ADDITION", "checks": {k: bool(v) for k, v in checks.items()},
            "spectrum_error": spectrum_error, "distance_error": distance_error,
            "periodic_return_error": float(return_error), "eigenline_residual": eigenline_error,
            "phase_over_pi": fractions.tolist(), "continuous_angle_degrees": np.degrees(line["angle"]).tolist(),
            "nearest_axis_defined_at_quarter_period": minimum["axis_defined"],
            "two_optima_distance_errors": [float(e) for e in fit_errors],
            "nearest_amplitude_changes_sign": bool(switch)}


def run_probe(progress=None):
    def note(message):
        if progress:
            progress(message)

    adm = rotor.adm_quadratic_derivation()
    note("Independent ADM curvature/action expansion completed.")
    old = backreaction.derive_the_tensor_mode_equation()
    cartan = {"omega2_radius2": float(old["omega_squared_times_a_squared"]),
              "anchors_pass": bool(old["the_validations_pass"]),
              "momentum_constraint": bool(old["the_momentum_constraint_holds"])}
    note("Existing Cartan/Einstein frequency and constraints checked.")
    algebra = algebra_controls()
    primary = free_motion()
    stationary = free_motion(speed=0.)
    amplitudes = [free_motion(amplitude=a) for a in (.01, .005, .0025)]
    speeds = [free_motion(speed=v) for v in (.2, .4, .8)]
    dynamics = dynamics_controls()
    inherited = inherited_response_control()
    relative_amplitude_error = max(float(np.max(abs(
        np.array(r["distances"]) / abs(r["amplitude"])
        - np.array(primary["distances"]) / abs(primary["amplitude"])))) for r in amplitudes)
    coefficient_error = abs(primary["distances_over_t2"][-1] / primary["predicted_coefficient"] - 1)
    checks = {
        "Q1 independent ADM round and linear anchors": adm["round_anchor"] and adm["linear_vanishes"],
        "Q1 ADM kinetic and Cartan frequency agree": adm["kinetic_is_frobenius"] and cartan["anchors_pass"]
            and cartan["momentum_constraint"] and abs(adm["omega2_radius2"] - cartan["omega2_radius2"]) < 1e-12
            and abs(adm["omega2_radius2"] - backreaction.TENSOR_MODE_FREQUENCY ** 2) < 1e-12,
        "Q1 matrix velocity and restricted action": max(algebra[k] for k in ("kinetic", "lagrangian", "adm_kinetic")) < 1e-10,
        "Q1 finite-difference pullback metric": algebra["metric_finite_difference"] < 1e-7,
        "Q1 relational rotation covariance and antipodal identity": max(algebra[k] for k in ("covariance", "antipodal")) < 1e-10,
        "Q2 omitted equation matches full tensor residual": algebra["normal_formula"] < 1e-10,
        "Q2 normal projector matches explicit biaxial basis": algebra["normal_basis"] < 1e-10,
        "Q2 free normal norm identity": algebra["normal_norm"] < 1e-10,
        "Q2 restricted radial and angular equations vanish": max(algebra[k] for k in ("radial_equation", "angular_equation")) < 1e-10,
        "Q2 primary full-field trajectory leaves the cone": min(primary["distances"]) > 1e-10,
        "Q2 second-order departure in time": all(3.8 < r < 4.2 for r in primary["halving_ratios"]),
        "Q2 predicted small-time coefficient": coefficient_error < .01,
        "Q2 minimization includes signed amplitude and all axes": max(r["eigenvalue_distance_error"] for r in amplitudes + speeds) < 1e-12,
        "Q2 stationary-director invariant control": max(stationary["distances"]) < 1e-12,
        "Q2 departure persists at linear field order": relative_amplitude_error < 1e-9,
        "Q2 every nonzero test speed has nonzero normal residual": all(r["free_normal_residual"] > 1e-5 for r in speeds),
        "Q2 independent unconstrained DOP853 flow": dynamics["free_ode_error"] < 1e-10,
        "Q2 conserved free energy and rotation charge": max(max(r["energy_residual"], r["charge_residual"]) for r in amplitudes + speeds) < 1e-10,
        "Q2 manufactured forced rotation solves full equation": dynamics["manufactured_ode_error"] < 1e-10
            and dynamics["manufactured_full_equation_residual"] < 1e-10,
        "Q2 manufactured rotation requires radial and normal stress": min(dynamics["manufactured_min_radial_force"], dynamics["manufactured_min_normal_force"]) > 1e-4,
        "Q2 frozen stationary amplitude also needs support": dynamics["frozen_amplitude_stationary_radial_force"] > .01,
        "Q3 repository source normalization and response": inherited["source_variation_identity_residual"] < 1e-12
            and inherited["errors_513_1025"][-1] < 1e-8
            and inherited["errors_513_1025"][-1] < inherited["errors_513_1025"][0],
    }
    checks = {k: bool(v) for k, v in checks.items()}
    frozen_checks = dict(checks)
    review = post_review_controls()
    checks.update(review["checks"])
    return {"public_preregistration": rotor.PUBLIC_PREREG,
            "base_commit": "cd7f3001581a37a2fe887ff4452340c45cb9f00a",
            "model": {"radius": 1., "kappa": 1., "volume": rotor.TensorModel().volume,
                      "action_normalization": rotor.TensorModel().normalization},
            "adm": adm, "cartan_reference": cartan, "algebra": algebra,
            "primary": primary, "stationary_director": stationary,
            "amplitude_controls": amplitudes, "speed_controls": speeds,
            "dynamics": dynamics, "inherited_response": inherited,
            "relative_amplitude_error": relative_amplitude_error,
            "primary_coefficient_relative_error": coefficient_error,
            "frozen_checks": frozen_checks, "post_review": review,
            "checks": checks, "passed": all(checks.values()),
            "physical_verdict": rotor.physical_verdict(checks, adm["kinetic_is_frobenius"] and algebra["normal_formula"] < 1e-10),
            "scope": "Unconstrained linear homogeneous ESU TT mode and its proposed uniaxial cone; not a no-go for every BAM field or driven solution.",
            "weight_selection": "NOT_DERIVED", "operational_source_readout": "NOT_DERIVED"}


def render(report):
    lines = ["# TT field-to-triangle-rotor test", "", report["scope"], "",
             "Public freeze: `" + report["public_preregistration"] + "`.", "",
             "Physical verdict: **" + report["physical_verdict"] + "**.", "",
             "| time | distance to full uniaxial cone | distance / time squared |",
             "|---:|---:|---:|"]
    p = report["primary"]
    for t, d, scaled in zip(p["times"], p["distances"], p["distances_over_t2"]):
        lines.append(f"| {t:g} | {d:.12g} | {scaled:.12g} |")
    lines += ["", f"Predicted coefficient: {p['predicted_coefficient']:.12g}.", "",
              "| criterion | pass |", "|---|---|"]
    lines += [f"| {key} | {value} |" for key, value in report["checks"].items()]
    lines += ["", f"Passed {sum(report['checks'].values())}/{len(report['checks'])} checks "
              f"({len(report['frozen_checks'])} frozen, {len(report['post_review']['checks'])} post-review).", "",
              "A passing obstruction check is not a successful rotor reduction. Phi selection and operational source readout remain unproved.", ""]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    report = run_probe(progress=lambda s: print(s, flush=True))
    summary = render(report)
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "probe.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        (args.output_dir / "probe.md").write_text(summary)
    print(summary)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

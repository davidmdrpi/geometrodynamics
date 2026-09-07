"""Derive and test a reciprocal scalar–TT projection; public freeze d8dc90d."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from geometrodynamics.bulk.tt_triangle_rotor import nearest_uniaxial, stf
from geometrodynamics.waves import reciprocal_scalar_tt as rt


def finite_difference_checks(model, state):
    # H has degree at most two in each individual coordinate, and the force
    # has degree at most two. Central differences are exact in real arithmetic
    # here; a larger step avoids cancellation without a Taylor-error tradeoff.
    h = 1e-4
    eye = np.eye(model.size)
    grad = np.array([(model.hamiltonian(state+h*e)-model.hamiltonian(state-h*e))/(2*h)
                     for e in eye])
    b, P, q, p = model.unpack(grad)
    expected = model.pack(P, -b, p, -q)
    got = model.rhs(0., state)
    jac = np.column_stack([(model.rhs(0., state+h*e)-model.rhs(0., state-h*e))/(2*h)
                            for e in eye])
    d = model.multiplet.dimension
    mixed = jac[5:10, 10:10+d]
    reciprocal = jac[10+d:, :5].T
    scale = max(1., np.linalg.norm(mixed))
    one_way_jac = np.column_stack([
        (model.rhs(0., state+h*e, reciprocal=False)
         - model.rhs(0., state-h*e, reciprocal=False))/(2*h) for e in eye])
    one_way_defect = one_way_jac[5:10, 10:10+d]-one_way_jac[10+d:, :5].T
    return {"hamilton_equations_scaled_error": float(np.linalg.norm(got-expected)/max(1., np.linalg.norm(got))),
            "mixed_derivative_scaled_error": float(np.linalg.norm(mixed-reciprocal)/scale),
            "one_way_mixed_derivative_defect": float(np.linalg.norm(one_way_defect)/scale)}


def action_expansion(model, state):
    b, P, q, p = model.unpack(state)
    direction = stf(np.array([[.5, .2, -.1], [.2, -.3, .15], [-.1, .15, -.2]]))
    direction /= np.linalg.norm(direction)
    scalar0 = .5*(p @ p-model.omega_scalar2*(q @ q))
    rows = []
    for scale in (.02, .01, .005):
        beta = scale*direction
        exact = model.static_scalar_lagrangian(beta, q, p)
        linear = scalar0+rt.components(beta) @ model.source(q)
        rows.append({"beta_scale": scale, "exact_static_scalar_L": float(exact),
                     "linear_scalar_L": float(linear), "remainder": float(abs(exact-linear))})
    ratios = [a["remainder"]/b["remainder"] for a, b in zip(rows, rows[1:])]
    return {"rows": rows, "halving_ratios": ratios,
            "scope": "static homogeneous anisotropic scalar action, including xi R phi^2"}


def covariance_checks(model, state):
    axis = np.array([2., -1., 3.])/math.sqrt(14)
    angle = .61
    R = expm(angle*rt.cross_matrix(axis))
    b, P, q, p = model.unpack(state)
    qr = model.multiplet.rotate(q, axis, angle)
    pr = model.multiplet.rotate(p, axis, angle)
    br = rt.components(R @ rt.tensor(b) @ R.T)
    Pr = rt.components(R @ rt.tensor(P) @ R.T)
    rotated = model.pack(br, Pr, qr, pr)
    m = np.array([1., 2., 3.])/math.sqrt(14)
    return {"hamiltonian_error": float(abs(model.hamiltonian(rotated)-model.hamiltonian(state))),
            "source_covariance_error": float(np.linalg.norm(rt.tensor(model.source(qr))
                                          - R @ rt.tensor(model.source(q)) @ R.T)),
            "Q_covariance_error": float(abs((R @ m) @ rt.tensor(br) @ (R @ m)-m @ rt.tensor(b) @ m)),
            "coherent_rotation_error": float(np.linalg.norm(model.multiplet.rotate(
                 model.multiplet.coherent(m), axis, angle)-model.multiplet.coherent(R @ m)))}


def free_control(model, scalar_zero=False):
    times = np.linspace(0, 4, 101)
    initial = rt.primary_data(model)
    if scalar_zero:
        initial[10:] = 0.
    b, P, q, p = model.unpack(initial)
    w, v = math.sqrt(model.omega_tensor2), math.sqrt(model.omega_scalar2)
    rows = [model.pack(b*math.cos(w*t)+P*math.sin(w*t)/(model.C*w),
                       P*math.cos(w*t)-model.C*w*b*math.sin(w*t),
                       q*math.cos(v*t)+p*math.sin(v*t)/v,
                       p*math.cos(v*t)-v*q*math.sin(v*t)) for t in times]
    computed = model.integrate(initial, times, rtol=1e-12, atol=1e-14)
    return float(np.max(np.abs(computed-np.array(rows))))


def constraint_certificate(model):
    """Nonzero omitted constraint source, with an independent analytic anchor.

    For q=A Re(z^3) normalized by sqrt(8/V), p=0: at x0=1 the
    improved T00 is 44 A^2/V. At x0=m.x_vec=0, phi and its first
    derivatives vanish, giving T00=0. A homogeneous background adjustment
    cannot cancel this contrast.
    """
    from geometrodynamics.waves.backreaction import stress_series
    state = rt.primary_data(model)
    b, P, q, p = model.unpack(state)
    m = np.array([1., 2., 3.])/math.sqrt(14)
    transverse = np.cross(m, [0., 0., 1.])
    transverse /= np.linalg.norm(transverse)
    points = np.array([[1., 0., 0., 0.], np.r_[0., transverse]])
    jets, _ = rt.scalar_jets(model, q, p, -model.omega_scalar2*q, points)
    rho = stress_series(jets)[:, 0, 0, 0]
    predicted = 44*.2**2/model.volume
    return {"north_energy_density": float(rho[0]), "transverse_energy_density": float(rho[1]),
            "predicted_contrast_44_amplitude2_over_volume": predicted,
            "certificate_error": float(abs(rho[0]-rho[1]-predicted)),
            "linear_homogeneous_TT_delta_G00": 0.,
            "linear_homogeneous_TT_delta_G0i": 0.,
            "scope": "omitted constraints of this TT-only ansatz; additional metric/support response is open"}


def run_probe(progress=lambda s: None):
    rng = np.random.default_rng(rt.SEED)
    model = rt.ReciprocalModel()
    initial = rt.primary_data(model)
    progress("harmonic representation and degree selection")
    algebra = {str(n): rt.harmonic_multiplet(n).algebra_checks() for n in (1, 3, 5)}
    m = np.array([1., 2., 3.])/math.sqrt(14)
    coherent = []
    for n in (1, 3):
        M = rt.ReciprocalModel(n)
        q = M.multiplet.coherent(m)
        expected = n*(n-1)*(np.outer(m, m)-np.eye(3)/3)
        coherent.append({"degree": n, "interaction_matrix_norm": float(np.linalg.norm(M.F)),
                         "quadrupole_error": float(np.linalg.norm(rt.tensor(M.source(q))-expected)),
                         "source": rt.tensor(M.source(q)).tolist()})
    progress("pointwise improved stress and independent scalar jets")
    stress = []
    for n in (1, 3):
        M = rt.ReciprocalModel(n)
        b, P, q, p = M.unpack(rt.primary_data(M))
        q = .13*rng.normal(size=M.multiplet.dimension)
        p = .07*rng.normal(size=M.multiplet.dimension)
        state = M.pack(b, P, q, p)
        for nr, na in ((8, 16), (12, 24)):
            d = rt.inherited_stress_diagnostics(M, state, nr, na)
            stress.append({"degree": n, **d})
    b, P, q, p = model.unpack(initial)
    points, _ = rt.sphere_quadrature(8, 16)
    qddot = model.unpack(model.rhs(0., initial))[3]
    jets, _ = rt.scalar_jets(model, q, p, qddot, points)
    wave_residual = (jets["dtt"]-jets["laplacian"]+jets["phi"]
                     + 2*np.einsum("ij,ptij->pt", rt.tensor(b), jets["hess"]))
    wave_error = float(np.max(np.abs(wave_residual)))
    progress("action variation, reciprocal force, and covariance")
    fd = finite_difference_checks(model, initial)
    expansion = action_expansion(model, initial)
    covariance = covariance_checks(model, initial)
    progress("coupled histories and one-way control")
    times = np.linspace(0, 4, 401)
    coarse = model.integrate(initial, times)
    fine = model.integrate(initial, times, rtol=1e-12, atol=1e-14)
    one_way = model.integrate(initial, times, rtol=1e-12, atol=1e-14, reciprocal=False)
    b, P, q, p = model.unpack(fine)
    beta = rt.tensor(b)
    energies = model.hamiltonian(fine)
    energy_drift = float(np.max(np.abs(energies-energies[0]))/abs(energies[0]))
    coarse_drift = float(np.max(np.abs(model.hamiltonian(coarse)-energies[0]))/abs(energies[0]))
    parts = model.energy_parts(fine)
    uniaxial = [nearest_uniaxial(B)["distance"] for B in beta]
    Q = np.einsum("i,tij,j->t", m, beta, m)
    q_one_way = model.unpack(one_way)[2]
    history = {
        "times": times.tolist(), "beta_components": b.tolist(),
        "scalar_coefficients": q.tolist(), "Q_m": Q.tolist(),
        "distance_to_uniaxial_cone": uniaxial,
        "energy_parts": {k: v.tolist() for k, v in parts.items()},
        "fine_relative_energy_drift": energy_drift,
        "coarse_relative_energy_drift": coarse_drift,
        "max_absolute_state_difference_between_tolerances": float(np.max(np.abs(coarse-fine))),
        "max_beta_frobenius": float(np.max(np.linalg.norm(beta, axis=(1, 2)))),
        "max_scalar_difference_from_one_way": float(np.max(np.linalg.norm(q-q_one_way, axis=1))),
        "max_tensor_difference_from_one_way": float(np.max(np.linalg.norm(b-one_way[:, :5], axis=1))),
        "one_way_relative_defect_in_reciprocal_H": float(np.max(np.abs(
            model.hamiltonian(one_way)-energies[0]))/abs(energies[0])),
        "scope": "initial-value solutions of the projected Hamiltonian; no future setting or triangle conditioning",
    }
    progress("free controls and omitted Einstein constraints")
    controls = {"n1_decoupled_error": free_control(rt.ReciprocalModel(1)),
                "scalar_zero_TT_error": free_control(model, scalar_zero=True)}
    constraints = []
    for idx in (0, 100, 400):
        rows = [rt.inherited_stress_diagnostics(model, fine[idx], nr, na)
                for nr, na in ((8, 16), (12, 24))]
        keys = ("mean_energy_density", "inhomogeneous_energy_rms", "momentum_density_rms")
        constraints.append({"time": float(times[idx]), "grids": rows,
                            "quadrature_disagreement": max(abs(rows[0][k]-rows[1][k]) for k in keys)})
    certificate = constraint_certificate(model)
    n = np.array([0., 0., 1.])
    identity_error = float(abs(Q[0]-.01*((m @ n)**2-1/3)))
    algebra_values = [v for row in algebra.values() for key, v in row.items()
                      if key not in ("dimension", "expected_dimension")]
    checks = {
        "complete harmonic multiplets and invariant derivative algebra":
            all(r["dimension"] == r["expected_dimension"] for r in algebra.values())
            and max(algebra_values) < 1e-10,
        "n1 null coupling and n3 field quadrupole":
            coherent[0]["interaction_matrix_norm"] < 1e-10
            and max(r["quadrupole_error"] for r in coherent) < 1e-10,
        "action source matches inherited improved stress on both grids":
            max(r["source_scaled_error"] for r in stress) < 1e-9,
        "pointwise scalar equation matches the modal equation": wave_error < 1e-9,
        "static matter-action remainder is quadratic":
            all(3.5 < r < 4.5 for r in expansion["halving_ratios"]),
        "Hamilton equations from independent energy differences":
            fd["hamilton_equations_scaled_error"] < 1e-7,
        "reciprocity holds and fails in the one-way control":
            fd["mixed_derivative_scaled_error"] < 1e-10
            and fd["one_way_mixed_derivative_defect"] > 1e-3,
        "ODE refinement and Hamiltonian conservation":
            history["max_absolute_state_difference_between_tolerances"] < 1e-8
            and max(energy_drift, coarse_drift) < 1e-8,
        "tensor remains within the frozen small-field range": history["max_beta_frobenius"] < .05,
        "free controls agree with independent harmonic solutions": max(controls.values()) < 1e-9,
        "action and source are covariant under common SO3 rotations": max(covariance.values()) < 1e-9,
        "omitted constraint sources converge and match an analytic certificate":
            max(r["quadrature_disagreement"] for r in constraints) < 1e-9
            and certificate["certificate_error"] < 1e-10,
        "Q_m has the uniaxial initial identity": identity_error < 1e-10,
    }
    checks = {k: bool(v) for k, v in checks.items()}
    return {"public_preregistration": rt.PUBLIC_PREREG, "seed": rt.SEED,
            "baseline": "22f77a373a67fcda91078ad0284be6a37c7ca20b",
            "scope": "Reciprocal variational ESU TT–scalar projection; full Einstein constraints, "
                     "localized source, and two-boundary triangle map remain open.",
            "model": {"degree": 3, "scalar_modes": 16, "tensor_components": 5,
                      "radius": 1., "kappa": 1., "C": model.C,
                      "omega_tensor2": model.omega_tensor2, "omega_scalar2": model.omega_scalar2},
            "algebra": algebra, "coherent_modes": coherent, "stress_checks": stress,
            "pointwise_wave_error": wave_error, "action_expansion": expansion,
            "hamiltonian_derivatives": fd, "covariance": covariance, "history": history,
            "free_controls": controls, "omitted_constraints": constraints,
            "constraint_certificate": certificate, "initial_Q_identity_error": identity_error,
            "checks": checks, "checks_passed": all(checks.values()), "verdict": rt.verdict(checks)}


def render(report):
    v = report["verdict"]
    lines = ["# Reciprocal scalar–TT projection", "", report["scope"], "",
             f"Public freeze: `{report['public_preregistration']}`. Seed: `{report['seed']}`.", "",
             "| Question | Verdict |", "|---|---|"]
    lines += [f"| {k} | {v[k]} |" for k in rt.VERDICT_FIELDS]
    if v["failed_checks"]:
        lines += ["", "Required checks failed; verdicts are UNRESOLVED:"]
        lines += [f"- {k}" for k in v["failed_checks"]]
    h = report["history"]
    lines += ["", "| History diagnostic | Value |", "|---|---:|"]
    for k in ("fine_relative_energy_drift", "max_absolute_state_difference_between_tolerances",
              "max_beta_frobenius", "max_scalar_difference_from_one_way",
              "max_tensor_difference_from_one_way", "one_way_relative_defect_in_reciprocal_H"):
        lines.append(f"| {k} | {h[k]:.12g} |")
    lines += ["", "| Time | Mean energy density | Inhomogeneous energy RMS | Momentum RMS |",
              "|---:|---:|---:|---:|"]
    for c in report["omitted_constraints"]:
        d = c["grids"][-1]
        lines.append(f"| {c['time']:g} | {d['mean_energy_density']:.12g} | "
                     f"{d['inhomogeneous_energy_rms']:.12g} | {d['momentum_density_rms']:.12g} |")
    lines += ["", "| Required check | Pass |", "|---|---|"]
    lines += [f"| {k} | {ok} |" for k, ok in report["checks"].items()]
    lines += ["", f"Passed {sum(report['checks'].values())}/{len(report['checks'])} numerical checks.", "",
              "The CLI failure/overwrite path is verified separately by an end-to-end test.", "",
              "No counting function, Born law, physical source-local readout, or retrocausal channel is inferred."]
    return "\n".join(lines)+"\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    report = run_probe(progress=lambda s: print(s, flush=True))
    summary = render(report)
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir/"probe.json").write_text(json.dumps(
            report, indent=2, allow_nan=False, default=lambda v: v.tolist())+"\n")
        (args.output_dir/"probe.md").write_text(summary)
    print(summary)
    return 0 if report["checks_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

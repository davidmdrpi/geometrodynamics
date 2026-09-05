"""Reproduce the conditional classical-equilibrium derivation (76ed50e).

Run with ``python -m experiments.closure_ledger.closure_equilibrium_probe``.
Use --output-dir to retain JSON and Markdown evidence. Failed numerical
criteria produce exit status 1; neither a narrative label nor proximity to
the quantum target is an acceptance criterion.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import erf

from geometrodynamics.bulk import closure_equilibrium as ce
from geometrodynamics.bulk.closure_current import minimal_rotation_lift
from geometrodynamics.bulk.mouth_spin_frame import qinv, qmul, spin_frame


GAMMAS = (math.pi / 4, 1.0, 3 * math.pi / 4)
BETAS = (16, 64, 256, 1024, 4096)


def _distance(a, b):
    return max(abs(a[k] - b[k]) for k in a)


def _table(p):
    return {f"{sa:+d},{sb:+d}": float(v) for (sa, sb), v in p.items()}


def _frame_residuals():
    rng = np.random.default_rng(281)
    errors = {"transport": 0.0, "lift_sign": 0.0,
              "spatial_rotation": 0.0, "initial_frame": 0.0}
    for _ in range(80):
        qs = rng.normal(size=(3, 4))
        qs /= np.linalg.norm(qs, axis=1)[:, None]
        q0, q1, qr = qs
        u, w = rng.normal(size=(2, 3))
        u, w = u / np.linalg.norm(u), w / np.linalg.norm(w)
        x = spin_frame(q0)[0]
        G = qmul(minimal_rotation_lift(w, x),
                 qmul(minimal_rotation_lift(u, w), minimal_rotation_lift(x, u)))
        initial = np.asarray(spin_frame(q0))
        returned = np.asarray(spin_frame(qmul(q0, qinv(G))))
        measured = float(np.sum((returned - initial) ** 2) / 16)
        errors["transport"] = max(errors["transport"], abs(measured - ce.rotor_potential(x, u, w)))
        negative = np.asarray(spin_frame(qmul(q0, qinv(-G))))
        errors["lift_sign"] = max(errors["lift_sign"], float(np.max(abs(returned - negative))))
        R = np.asarray(spin_frame(qr))
        errors["spatial_rotation"] = max(errors["spatial_rotation"],
                                         abs(measured - ce.rotor_potential(R @ x, R @ u, R @ w)))
        other = np.asarray(spin_frame(qmul(q1, qinv(G)))) - np.asarray(spin_frame(q1))
        errors["initial_frame"] = max(errors["initial_frame"], abs(float(np.sum(other ** 2) / 16) - measured))
    return {k: float(v) for k, v in errors.items()}


def run_probe():
    checks = []

    def check(name, value, bound):
        checks.append({"name": name, "error": float(value), "bound": bound,
                       "passed": bool(np.isfinite(value) and value < bound)})

    frame = _frame_residuals()
    check("P1 quaternion frame identity and invariances", max(frame.values()), 1e-11)
    thermal, refinements, quartic = [], [], []
    for gamma in GAMMAS:
        target = ce.joint_probabilities(gamma, None)
        for beta in BETAS:
            p = ce.joint_probabilities(gamma, beta)
            mass_errors = []
            for c in (-math.cos(gamma), math.cos(gamma)):
                scaled = math.sqrt(beta / (2 * math.pi)) * ce.canonical_partition(c, beta)
                mass_errors.append(scaled / (ce.limiting_mass(c) / (4 * math.pi)) - 1)
            thermal.append({"gamma": gamma, "beta": beta, "probabilities": _table(p),
                            "E": ce.correlation(p), "joint_error": _distance(p, target),
                            "scaled_mass_relative_errors": mass_errors})
        base = ce.joint_probabilities(gamma, 4096)
        refinements.append({"gamma": gamma,
                            "normal_error": _distance(base, ce.joint_probabilities(gamma, 4096, n_normal=256)),
                            "azimuth_error": _distance(base, ce.joint_probabilities(gamma, 4096, n_azimuth=2048))})
        p4 = ce.joint_probabilities(gamma, 4096, quartic=2)
        rel4 = [math.sqrt(4096 / (2 * math.pi)) * ce.canonical_partition(c, 4096, quartic=2)
                / (ce.limiting_mass(c) / (4 * math.pi)) - 1
                for c in (-math.cos(gamma), math.cos(gamma))]
        quartic.append({"gamma": gamma, "joint_error": _distance(p4, target),
                        "scaled_mass_relative_errors": rel4,
                        "finite_beta_16_change": _distance(ce.joint_probabilities(gamma, 16),
                                                          ce.joint_probabilities(gamma, 16, quartic=2))})
    final = [r for r in thermal if r["beta"] == 4096]
    check("P2 limiting joint law", max(r["joint_error"] for r in final), 2e-3)
    check("P2 scaled partition masses", max(abs(e) for r in final for e in r["scaled_mass_relative_errors"]), 0.01)
    check("P2 independent coordinate refinements",
          max(max(r["normal_error"], r["azimuth_error"]) for r in refinements), 1e-4)
    mc = ce.monte_carlo_joint(1.0, 16, n_samples=120000)
    reference = ce.joint_probabilities(1.0, 16)
    zscore = max(abs(mc["probabilities"][k] - reference[k]) / mc["standard_errors"][k] for k in reference)
    check("P2 independent full-sphere Monte Carlo (standard errors)", zscore, 6.0)
    check("P2 Monte Carlo correlation (standard errors)",
          abs(mc["correlation"] - ce.correlation(reference)) / mc["correlation_standard_error"], 6.0)
    mc["probabilities"], mc["standard_errors"] = _table(mc["probabilities"]), _table(mc["standard_errors"])

    gaussian_error, normal_E, covariance_error = 0.0, 0.0, 0.0
    for beta in (0, 1, *BETAS):
        normal_E = max(normal_E, abs(ce.correlation(ce.joint_probabilities(1.0, beta, model="normal"))))
        for c in (-0.8, 0.0, 0.8):
            a = beta * (1 - c * c) / 2
            exact = 1.0 if a == 0 else math.sqrt(math.pi) * erf(math.sqrt(a)) / (2 * math.sqrt(a))
            normal = ce.canonical_partition(c, beta, model="normal")
            gaussian_error = max(gaussian_error, abs(normal - exact))
            covariance_error = max(covariance_error, abs(normal - ce.canonical_partition(c, beta, model="normal_covariant")))
    check("P3 normal residual exact Gaussian", gaussian_error, 1e-12)
    check("P3 normal residual zero correlation", normal_E, 1e-12)
    check("P4 covariance of finite-temperature partition", covariance_error, 1e-12)
    rng = np.random.default_rng(282)
    x = rng.normal(size=(1000, 3))
    x /= np.linalg.norm(x, axis=1)[:, None]
    u, w = np.array([0., 0., 1.]), np.array([math.sin(1.), 0., math.cos(1.)])
    energy_error = float(np.max(abs(ce.rotor_potential(x, u, w, model="normal")
                                    - ce.rotor_potential(x, u, w, model="normal_covariant"))))
    check("P4 covariance of pointwise energy", energy_error, 1e-12)
    c, eps = math.cos(1.), 0.25
    independent_mass = quad(lambda t: 1 / ((1 + eps * math.sqrt(2 * (1 + c)) * math.cos(t))
                                          * math.sqrt(1 - c * c)), 0, 2 * math.pi, epsabs=1e-11)[0]
    check("P4 fixed-stiffness analytic mass", abs(independent_mass - ce.limiting_mass(c, model="normal_rescaled")), 1e-10)
    rescaled_limit = ce.joint_probabilities(1.0, None, model="normal_rescaled")
    rescaled = {"E_limit": ce.correlation(rescaled_limit),
                "E_beta4096": ce.correlation(ce.joint_probabilities(1.0, 4096, model="normal_rescaled"))}
    check("P4 fixed-stiffness limiting joint law",
          _distance(rescaled_limit, ce.joint_probabilities(1.0, 4096, model="normal_rescaled")), 2e-3)
    checks.append({"name": "P4 fixed stiffness changes the ensemble", "passed": abs(rescaled["E_limit"]) > 1e-3})
    check("P5 quartic limiting joint law", max(r["joint_error"] for r in quartic), 2e-3)
    check("P5 quartic scaled partition masses", max(abs(e) for r in quartic for e in r["scaled_mass_relative_errors"]), 0.01)
    checks.append({"name": "P5 quartic changes finite-temperature response",
                   "passed": min(r["finite_beta_16_change"] for r in quartic) > 1e-4})

    root = math.acos(-math.sqrt((1 + c) / 2))
    punctures = []
    for radius in (0.1, 0.05, 0.025):
        mass = sum(quad(lambda t: ce.limiting_density(t, c), center - radius, center + radius,
                        points=[center], epsabs=1e-12)[0] for center in (root, 2 * math.pi - root))
        punctures.append({"radius": radius, "mass": mass, "relative_error_to_2eps2": mass / (2 * radius ** 2) - 1})
    check("P2 puncture excision asymptotic", max(abs(r["relative_error_to_2eps2"]) for r in punctures), 1e-3)
    priors, marginal_error = [], 0.0
    for ratio in (0.5, 1.0, 2.0):
        p = ce.joint_probabilities(1.0, 64, prior_ratio=ratio)
        priors.append({"like_unlike_prior_ratio": ratio, "E": ce.correlation(p), "probabilities": _table(p)})
        marginal_error = max(marginal_error, *(abs(sum(v for s, v in p.items() if s[i] == 1) - 0.5) for i in (0, 1)))
    check("P6 paired half marginals", marginal_error, 1e-12)
    checks.append({"name": "P6 prior dependence survives", "passed": priors[-1]["E"] - priors[0]["E"] > 0.4})
    return {"local_preregistration_commit": "76ed50e",
            "published_preregistration_commit": "d83d46aab8cfe87b0a0adc7f9674401be455d74f",
            "original_base_commit": "92a915bfaaabd02564accd7467b81cedb1ee8c16",
            "integrated_main_commit": "e96d48abf57fc97f135bb323fe0da116793d1077",
            "assumptions": ["round rotor with identical sector inertia", "canonical equilibrium preparation",
                            "specified classical frame-restoring energy", "equal sector priors except P6",
                            "geodesic triangle with singlet partner sign"],
            "open": ["BAM origin of the coupling", "equilibration and local detector implementation",
                     "sector prior", "composition and event readout", "source-readout causality"],
            "quadrature": {"n_normal": 128, "n_azimuth": 1024, "domain": "whole sphere, no closure cutoff"},
            "checks": checks, "passed": sum(r["passed"] for r in checks), "total": len(checks),
            "frame_residuals": frame, "thermal": thermal, "refinements": refinements,
            "monte_carlo": mc, "normal_gaussian_error": float(gaussian_error),
            "rescaled": rescaled, "quartic": quartic, "punctures": punctures, "priors": priors,
            "standard_angle_CHSH_limit": ce.standard_chsh(),
            "standard_angle_CHSH_beta4096": ce.standard_chsh(4096)}


def render(report):
    lines = ["# Classical closure equilibrium probe", "",
             f"Local preregistration: `{report['local_preregistration_commit']}`. "
             f"Checks: **{report['passed']}/{report['total']}**.", "",
             f"GitHub publication after calculation: `{report['published_preregistration_commit']}`.", "",
             "Conditional model: " + "; ".join(report["assumptions"]) + ".", "",
             "| criterion | passed | error | bound |", "|---|---|---|---|"]
    for r in report["checks"]:
        lines.append(f"| {r['name']} | {r['passed']} | {r.get('error', '')} | {r.get('bound', '')} |")
    lines += ["", "Full-sphere partition integrals; beta = K/(k_B T).", "",
              "| gamma | beta | E | max joint error to limit | max scaled mass relative error |",
              "|---|---|---|---|---|"]
    for r in report["thermal"]:
        lines.append(f"| {r['gamma']:.8f} | {r['beta']} | {r['E']:.10f} | {r['joint_error']:.3e} | "
                     f"{max(abs(e) for e in r['scaled_mass_relative_errors']):.3e} |")
    mc = report["monte_carlo"]
    lines += ["", f"Independent Monte Carlo: E = {mc['correlation']:.9f} ± "
              f"{mc['correlation_standard_error']:.9f} (one standard error; "
              f"n = {mc['n_samples']}, seed = {mc['seed']}). This integrates an equilibrium law; "
              "it does not simulate equilibration.", "",
              f"Standard-angle CHSH: limit {report['standard_angle_CHSH_limit']:.10f}; "
              f"beta = 4096: {report['standard_angle_CHSH_beta4096']:.10f}. No global maximum asserted.", "",
              "Open: " + "; ".join(report["open"]) + ".", "",
              "The positive measure follows from the stated energy and preparation. "
              "The closure locus alone does not fix it. Detailed controls are in probe.json."]
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, help="retain probe.json and probe.md")
    args = parser.parse_args(argv)
    report = run_probe()
    output = render(report)
    print(output)
    if args.output_dir is not None:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "probe.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        (args.output_dir / "probe.md").write_text(output, encoding="utf-8")
    return 0 if all(r["passed"] for r in report["checks"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())

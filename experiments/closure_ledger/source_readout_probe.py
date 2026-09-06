"""Frozen source-readout, conditional grounding, and cross-round controls.

Run with python -m experiments.closure_ledger.source_readout_probe.
An informative synthetic record is not an operational BAM channel.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from scipy.spatial.transform import Rotation
from scipy.special import roots_legendre

from geometrodynamics.bulk.source_readout import (
    FIXED_ANALYZER, FUTURE_ANALYZERS, source_circle, source_sphere,
    record_statistics, orthogonal_noiseless_tail,
)
from geometrodynamics.bulk.closure_grounding import (
    tube_frame_energy, momentum_integrals, inertia_controls,
)
from geometrodynamics.bulk.closure_consistency import run_consistency


PUBLIC_PREREG = "d258bb14e73674fd6ecd7ff6f5d2a46edf20eeab"


def duffing_pointer_demo(pointer_momentum=0., strength=.4):
    """Append h(t) P q_osc^2 to the repository's actual conservative source.

    q_osc is the LOCAL SCALAR OSCILLATOR, not the quaternion frame variable.
    Shoot a periodic orbit, then preserve both its boundary values at P=0.
    Final pointer position is an output; imposing it would be another problem.
    """
    from experiments.closure_ledger.hamiltonian_source_eigenhistory_probe import (
        _rhs_red, _h_red, _W0,
    )
    wr, coupling, amplitude = 2.738858, .3203, .4
    rhs = _rhs_red(wr, coupling)
    matrix = np.array([[wr ** 2, coupling], [coupling, _W0 ** 2]])
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    seed = [amplitude * eigenvectors[1, 0] / eigenvectors[0, 0],
            math.pi / math.sqrt(eigenvalues[0])]
    options = dict(method="DOP853", rtol=2e-12, atol=2e-14)

    def shoot(parameters):
        q0, half_period = parameters
        sol = solve_ivp(rhs, (0, half_period), [amplitude, 0, q0, 0], **options)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.y[[1, 3], -1]

    shooting = root(shoot, seed, tol=1e-10)
    if not shooting.success and np.max(np.abs(shoot(shooting.x))) > 1e-10:
        raise RuntimeError(shooting.message)
    q0, half_period = shooting.x
    period = 2 * half_period
    initial = np.array([amplitude, 0., q0, 0.])
    times = np.linspace(0., period, 1001)
    baseline = solve_ivp(rhs, (0., period), initial, t_eval=times, **options)

    def pulse(t):
        return 2 * strength / period * np.sin(math.pi * t / period) ** 2

    def coupled(t, y):
        source = np.array(rhs(t, y[:4]))
        source[3] -= 2 * pulse(t) * y[5] * y[2]  # retain reciprocal force
        return np.r_[source, y[5] + pulse(t) * y[2] ** 2, 0.]

    measured = solve_ivp(coupled, (0., period), np.r_[initial, 0., pointer_momentum],
                         t_eval=times, **options)
    if not baseline.success or not measured.success:
        raise RuntimeError("periodic source integration failed")
    energy = np.array([_h_red(y[:4], wr, coupling) + .5 * y[5] ** 2
                       + pulse(t) * y[5] * y[2] ** 2
                       for t, y in zip(times, measured.y.T)])
    return {"period": float(period), "source_q0": float(q0),
            "pointer_initial_momentum": float(pointer_momentum),
            "baseline_closure_error": float(np.max(np.abs(baseline.y[:, -1] - initial))),
            "source_endpoint_change": float(np.max(np.abs(measured.y[:4, -1] - baseline.y[:, -1]))),
            "source_trajectory_change": float(np.max(np.abs(measured.y[:4] - baseline.y))),
            "pointer_shift": float(measured.y[4, -1]),
            "pointer_momentum_change": float(np.max(np.abs(measured.y[5] - pointer_momentum))),
            "energy_range": float(np.ptp(energy)),
            "observable": "q_osc^2, not a function of the spin-frame x"}


def run_probe():
    checks, records = {}, []

    def check(name, value):
        checks[name] = bool(value)

    for beta in (None, 64.):
        pair = []
        for choice, b in enumerate(FUTURE_ANALYZERS):
            x, w = (source_circle(FIXED_ANALYZER, b) if beta is None
                    else source_sphere(FIXED_ANALYZER, b, beta))
            sharp = record_statistics(x, w)
            noisy = record_statistics(x, w, noise=.15)
            row = {"beta": beta, "choice": choice, "noiseless": sharp, "noisy": noisy}
            records.append(row)
            pair.append(row)
            check(f"A1 mean beta={beta} choice={choice}", abs(sharp["mean"]) < 1e-12)
            check(f"A1 sign beta={beta} choice={choice}", abs(sharp["positive_probability"] - .5) < 1e-12)
            if beta is None:
                check(f"A1 exact variance choice={choice}", abs(sharp["variance"] - (.1, .4)[choice]) < 1e-12)
                check(f"A2 noiseless antiderivative choice={choice}",
                      abs(sharp["tail_probability"] - orthogonal_noiseless_tail()[choice]) < 1e-4)
            else:
                for coordinate, sizes in (("normal", (256, 512)), ("azimuth", (128, 1024))):
                    xr, wr = source_sphere(FIXED_ANALYZER, b, beta, *sizes)
                    refined = record_statistics(xr, wr)
                    check(f"A1 refine {coordinate} choice={choice}",
                          abs(sharp["variance"] - refined["variance"]) < 1e-4)
        check(f"A1 variance gap beta={beta}",
              pair[1]["noiseless"]["variance"] - pair[0]["noiseless"]["variance"] > .15)
        check(f"A2 noisy tail gap beta={beta}",
              pair[1]["noisy"]["tail_probability"] - pair[0]["noisy"]["tail_probability"] > .2)
    pointer, recoil = duffing_pointer_demo(), duffing_pointer_demo(.001)
    check("A3 periodic source orbit", pointer["baseline_closure_error"] < 1e-9)
    check("A3 zero-P source preserved", pointer["source_trajectory_change"] < 1e-9
          and pointer["source_endpoint_change"] < 1e-9)
    check("A3 nonzero pointer record", pointer["pointer_shift"] > 1e-6)
    check("A3 zero-P energy and momentum", pointer["energy_range"] < 1e-9
          and pointer["pointer_momentum_change"] < 1e-12)
    check("A3 finite-P recoil", recoil["source_endpoint_change"] > 1e-6)
    rng = np.random.default_rng(20260905)
    tubes = [tube_frame_energy(R) for R in Rotation.random(12, random_state=rng).as_matrix()]
    check("B1 finite-difference parent", max(t["fd_error"] for t in tubes) < 1e-11)
    check("B1 existing DtN static limit", max(t["dtn_error"] for t in tubes) < 1e-10)
    momenta = momentum_integrals()
    check("B3 Gaussian momentum integral", abs(momenta["gaussian"] - 2 * math.pi) < 1e-10)
    check("B3 quartic momentum integral", abs(momenta["quartic"] - math.pi ** 1.5) < 1e-10)
    metrics = [inertia_controls(x) for x in rng.normal(size=(12, 3))]
    check("B3 variable inertia volume", max(abs(m["variable_volume"] - m["expected_variable_volume"]) for m in metrics) < 1e-12)
    check("B3 anisotropic unit volume", max(abs(m["anisotropic_volume"] - 1) for m in metrics) < 1e-12)
    z, wz = roots_legendre(48)
    volumes = np.array([inertia_controls([math.sqrt(1 - zi * zi), 0, zi])["variable_volume"] for zi in z])
    inertia_mean = float((wz * z) @ volumes / (wz @ volumes))
    check("B3 nonuniform prior mean", abs(inertia_mean - 1 / 6) < 1e-12)
    consistency = run_consistency()
    for key, error in consistency["max_residuals"].items():
        tolerance = 1e-7 if key in ("positive_mass", "phase_gradient", "frame_hessian") else 1e-9
        check("C " + key, error < tolerance)
    return {"public_preregistration": PUBLIC_PREREG, "checks": checks,
            "passed": all(checks.values()), "records": records,
            "pointer": pointer, "finite_P_control": recoil, "tube_controls": tubes,
            "momentum_integrals": momenta, "variable_inertia_mean_z": inertia_mean,
            "consistency": consistency,
            "not_supplied": ["frame x to local BAM field family", "derived physical pointer/preparation",
                             "instrument-modified ensemble", "independently selectable future settings"]}


def render(report):
    lines = ["# Source-readout and grounding probe", "",
             f"Public preregistration: `{report['public_preregistration']}`.", "",
             "These are conditional laws and added classical couplings. The operational BAM gate remains open.", "",
             "| beta | choice | mean F | variance F | P(F>0) | P(abs(Y)>0.6), noise=0.15 |",
             "|---|---|---:|---:|---:|---:|"]
    for row in report["records"]:
        r = row["noiseless"]
        lines.append(f"| {row['beta']} | {row['choice']} | {r['mean']:.8g} | {r['variance']:.10f} | {r['positive_probability']:.10f} | {row['noisy']['tail_probability']:.10f} |")
    lines += ["", "The finite-pulse pointer reads q_osc^2, not the spin-frame variable.", "",
              "```json", json.dumps(report["pointer"], indent=2), "```", "",
              "| criterion | pass |", "|---|---|"]
    lines += [f"| {name} | {ok} |" for name, ok in report["checks"].items()]
    lines += ["", f"Passed {sum(report['checks'].values())}/{len(report['checks'])} criteria.", "",
              "The JSON includes all cross-round residuals, masses, recoil, and field-parent controls.", ""]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    report = run_probe()
    summary = render(report)
    if args.output_dir is not None:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "probe.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        (args.output_dir / "probe.md").write_text(summary)
    print(summary)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

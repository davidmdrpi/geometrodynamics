"""Reproduce the public normal-reversal experiment (freeze 533f925)."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from geometrodynamics.bulk import antipodal_normal_throat as ant

REQUIRED = ("geometry", "spectral_selection", "continuum_convergence",
            "angular_symmetry", "susceptibility", "seam", "momentum",
            "mouth_pair", "energy", "response_scaling", "independent_integrator")


def structural_report():
    rng = np.random.default_rng(ant.SEED)
    errors = dict(spectrum=0., selection=0., angular=0., canonical=0., response=0.)
    forbidden = 0.
    count = 0
    for size in (16, 24):
        for n in (1, 2):
            model = ant.Throat(size, n)
            k = model.frozen_matrix()
            for kind in ("antipodal", "reflective"):
                z = rng.normal(size=(4, model.dim))*0.05
                d = model.deck(np.eye(model.dim), kind).T
                w = model.symplectic()
                errors["canonical"] = max(errors["canonical"], float(abs(d.T@w@d-w).max()))
                sign = 1 if kind == "antipodal" else -1
                errors["angular"] = max(errors["angular"], float(abs(
                    model.angular_momentum(model.deck(z, kind))
                    - sign*model.angular_momentum(z)).max()))
                for eta in (1, -1):
                    if kind == "reflective":
                        forbidden = max(forbidden, float(abs(model.angular_momentum(
                            model.project(z, kind, eta))).max()))
            for eta in (1, -1):
                for m in range(size//2+1):
                    for sine in ([False, True] if 0 < m < size//2 else [False]):
                        epsilon = ant.normal_parity(m, n, eta)
                        g = model.mode(m, epsilon, sine).ravel()
                        errors["spectrum"] = max(errors["spectrum"], float(abs(
                            k@g-model.frequency_sq(m, epsilon)*g).max()))
                        x = np.zeros((1, 2, 2, size)); x[0, 0] = g.reshape(2, size)
                        z = model.pack(x, np.zeros_like(x), np.zeros((1, 2)), np.zeros((1, 2)))
                        errors["selection"] = max(errors["selection"], float(abs(
                            model.deck(z, eta=eta)-z).max()))
                        count += 1
                for m in (0, 1):
                    epsilon = ant.normal_parity(m, n, eta)
                    g = model.mode(m, epsilon).ravel()
                    e = np.array([1., epsilon])/np.sqrt(2)
                    for omega in (0.3, 1.4):
                        direct = g@np.linalg.solve(k-omega**2*np.eye(2*size), g)
                        predicted = e@model.susceptibility(m, omega)@e
                        errors["response"] = max(errors["response"], abs(float(direct-predicted)))
    convergence = []
    for m in (1, 2, 3):
        sizes = np.array([16, 24, 32, 48])
        error = [abs(ant.Throat(int(s)).frequency_sq(m, 1)-(m*m+1+.25+2)) for s in sizes]
        convergence.append({"m": m, "errors": error,
                            "order": float(-np.polyfit(np.log(sizes), np.log(error), 1)[0])})
    return {"errors": errors, "resolved_modes": count,
            "reflective_fixed_total_momentum": forbidden, "convergence": convergence,
            "surface_determinants": [float(np.linalg.det(ant.jacobian(k)[:2, :2]))
                                     for k in ("antipodal", "reflective")],
            "throat_determinants": [float(np.linalg.det(ant.jacobian(k)))
                                    for k in ("antipodal", "reflective")],
            "theta_pushforward": [ant.jacobian(k)[:, 1].tolist()
                                  for k in ("antipodal", "reflective")]}


def responding_report(size, n, duration=40.):
    model = ant.Throat(size, n)
    cases = [(m, a) for m in (0, 1) for a in ant.AMPLITUDES]
    starts = np.array([model.travelling(a, m) for m, a in cases])
    hist = ant.flow(model, starts, duration)
    rows = []
    for j, (m, a) in enumerate(cases):
        z = hist[:, j]
        x, _, q, _ = model.unpack(z)
        momentum = model.angular_momentum(z)
        per_face = model.angular_momentum(z, per_face=True)
        loads = model.loads(z)
        e = model.energy(z)
        rows.append({"m": m, "n": n, "amplitude": a,
                     "normal_parity": ant.normal_parity(m, n),
                     "initial_energy": float(e[0]), "initial_momentum": float(momentum[0]),
                     "max_q": float(abs(q).max()), "min_q": float(q.min()),
                     "seam_error": float(abs(model.deck(z)-z).max()/np.linalg.norm(z[0])),
                     "momentum_error": float(abs(momentum-momentum[0]).max()/abs(momentum[0])),
                     "energy_error": float(abs(e-e[0]).max()/e[0]),
                     "mouth_difference": float(abs(q[:, 0]-q[:, 1]).max()),
                     "relative_face_momentum_difference": float(abs(per_face[:, 0]-per_face[:, 1]).max()/abs(momentum[0])),
                     "signed_load_sum": float(abs(loads[:, 0]-loads[:, 1]).max()),
                     "minimum_face_momentum_ratio": float((per_face/per_face[0]).min()),
                     "max_local_density": float(np.sum(x*x, axis=1).max())})
    # Independent comparison to the same initial states with frozen q=0.
    short = ant.flow(model, starts, 4.)
    frozen = ant.flow(model, starts, 4., respond=False)
    scaling = []
    for mi, m in enumerate((0, 1)):
        sl = slice(mi*4, mi*4+4)
        mouth = abs(short[:, sl, -4:-2]).max(axis=(0, 2))
        correction = np.linalg.norm(short[:, sl, :model.field_size]
                                    - frozen[:, sl, :model.field_size], axis=-1).max(axis=0)
        scaling.append({"m": m, "mouth_peaks": mouth.tolist(), "field_corrections": correction.tolist(),
                        "mouth_slope": float(np.polyfit(np.log(ant.AMPLITUDES), np.log(mouth), 1)[0]),
                        "field_slope": float(np.polyfit(np.log(ant.AMPLITUDES), np.log(correction), 1)[0])})
    return {"sites": size, "n": n, "duration": duration, "rows": rows, "scaling": scaling}


def numerical_report(duration=40.):
    model = ant.Throat(16, 1)
    z0 = model.travelling(.4)
    coarse = ant.flow(model, z0, duration)[:, 0]
    fine = ant.flow(model, z0, duration, h=.005, sample=20)[:, 0]
    times = np.linspace(0., duration, len(coarse))
    ref = solve_ivp(lambda t, z: model.rhs(z)[0], (0., duration), z0,
                    method="DOP853", rtol=1e-10, atol=1e-12, t_eval=times)
    if not ref.success or len(ref.t) != len(times):
        raise RuntimeError(ref.message)
    scale = np.linalg.norm(ref.y[:, -1])
    ec = np.linalg.norm(coarse[-1]-ref.y[:, -1])/scale
    ef = np.linalg.norm(fine[-1]-ref.y[:, -1])/scale
    return {"duration": duration, "coarse_final_error": float(ec), "fine_final_error": float(ef),
            "improvement": float(ec/ef),
            "dop853_momentum_error": float(abs(model.angular_momentum(ref.y.T)
                                               - model.angular_momentum(z0)[0]).max()
                                              /abs(model.angular_momentum(z0)[0]))}


def evaluate(structure, runs, numerics):
    err = structure["errors"]
    rows = [r for run in runs for r in run["rows"]]
    slopes = [r for run in runs for r in run["scaling"]]
    checks = {
        "geometry": structure["surface_determinants"] == [1., -1.]
                    and structure["throat_determinants"] == [-1., -1.]
                    and structure["theta_pushforward"] == [[0., 1., 0.], [0., -1., 0.]],
        "spectral_selection": err["spectrum"] < 1e-10 and err["selection"] < 1e-12,
        "continuum_convergence": all(1.8 < r["order"] < 2.2 for r in structure["convergence"]),
        "angular_symmetry": max(err["angular"], err["canonical"], structure["reflective_fixed_total_momentum"]) < 1e-10,
        "susceptibility": err["response"] < 1e-10,
        "seam": bool(rows) and max(r["seam_error"] for r in rows) < 1e-9,
        "momentum": bool(rows) and max(r["momentum_error"] for r in rows) < 1e-9,
        "mouth_pair": bool(rows) and max(max(r["mouth_difference"], r["relative_face_momentum_difference"], r["signed_load_sum"]) for r in rows) < 1e-9,
        "energy": bool(rows) and max(r["energy_error"] for r in rows) < 2e-3,
        "response_scaling": bool(slopes) and all(abs(r["mouth_slope"]-2) < .15 and abs(r["field_slope"]-3) < .2 for r in slopes),
        "independent_integrator": numerics["fine_final_error"] < 1e-2
                                  and numerics["coarse_final_error"] < 1e-2
                                  and numerics["improvement"] > 1,
    }
    return {k: bool(v) for k, v in checks.items()}


def transfer_report():
    """Post-freeze exploration: does circulation mean net inter-mouth transfer?

    Add generic deck-fixed perturbations to excite both normal channels.
    Check a non-deck-fixed control to exclude a trivially zero flux routine.
    These measurements are not additional preregistered predictions.
    """
    model = ant.Throat(16, 1)
    rng = np.random.default_rng(ant.SEED+1)
    z = model.travelling(.2)
    z[:2*model.field_size] += .01*rng.normal(size=2*model.field_size)
    z = model.project(z)[0]
    hist = ant.flow(model, z, 4.)[:, 0]
    face = model.angular_momentum(hist, per_face=True)
    x = np.zeros((1, 2, 2, model.sites))
    x[0, 0, 0] = .1
    x[0, 1, 1] = .1
    control = model.pack(x, np.zeros_like(x), np.zeros((1, 2)), np.zeros((1, 2)))
    ch = ant.flow(model, control, 1.)[:, 0]
    return {"scope": "POST_FREEZE_EXPLORATION_NOT_A_FROZEN_GATE",
            "max_local_normal_flux": float(abs(model.normal_angular_flux(hist, local=True)).max()),
            "max_net_normal_flux": float(abs(model.normal_angular_flux(hist)).max()),
            "max_face_momentum_change": float(abs(face-face[0]).max()),
            "control_initial_net_flux": float(model.normal_angular_flux(control)[0]),
            "control_face_momentum_change": float(abs(model.angular_momentum(ch, per_face=True)
                                                       - model.angular_momentum(control, per_face=True)).max()),
            "interpretation": "CIRCULATION_SURVIVES_BUT_EXACT_DECK_SYMMETRY_FORBIDS_NET_NORMAL_ANGULAR_TRANSFER"}


def run_probe(duration=40., progress=lambda s: None):
    structure = structural_report()
    runs = []
    for size in (16, 24):
        for n in (1, 2):
            progress(f"N={size}, n={n}: responding mouths")
            runs.append(responding_report(size, n, duration))
    numerics = numerical_report(duration)
    checks = evaluate(structure, runs, numerics)
    return {"preregistration": ant.PREREG, "seed": ant.SEED, "structure": structure,
            "runs": runs, "numerics": numerics, "checks": checks,
            "transfer_exploration": transfer_report(),
            "verdict": ant.classify(checks, REQUIRED)}


def render(r):
    rows = [row for run in r["runs"] for row in run["rows"]]
    lines = ["# Normal reversal and circulating real waves", "",
             f"Public freeze: `{r['preregistration']}`.", "",
             f"Numerical validity: **{r['verdict']['numerical_validity']}**.", "",
             "| Check | Pass |", "|---|---|"]
    lines += [f"| {k} | {v} |" for k, v in r["checks"].items()]
    lines += ["", "| N | n | m | amplitude | normal parity | max opening | relative L error | energy excursion |",
              "|---|---|---|---|---|---|---|---|"]
    for run in r["runs"]:
        for v in run["rows"]:
            lines.append(f"| {run['sites']} | {v['n']} | {v['m']} | {v['amplitude']} | {v['normal_parity']} | {v['max_q']:.6g} | {v['momentum_error']:.3g} | {v['energy_error']:.3g} |")
    lines += ["", "Angular momentum is that of one real harmonic. Signed loads are waist",
              "loads, not a demonstrated bulk momentum balance or electric charge.", "",
              "Independent integrator:", "```json", json.dumps(r["numerics"], indent=2), "```", ""]
    lines += ["Post-freeze transfer exploration:", "```json",
              json.dumps(r["transfer_exploration"], indent=2), "```", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    r = run_probe(progress=lambda s: print(s, flush=True))
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir/"probe.json").write_text(json.dumps(r, indent=2, allow_nan=False)+"\n")
        (args.output_dir/"probe.md").write_text(render(r))
    print(render(r))
    return 0 if r["verdict"]["numerical_validity"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())

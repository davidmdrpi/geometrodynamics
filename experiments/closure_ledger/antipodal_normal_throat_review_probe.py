"""Amended four-cell experiment, publicly registered at 31cfc67.

The reviewer supplied the pattern before this amendment. These are
confirmation and stronger validation runs, not blind discoveries.
The original runner and first archive retain their original thresholds.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from geometrodynamics.bulk import antipodal_normal_throat as ant
from experiments.closure_ledger import antipodal_normal_throat_probe as original

AMENDMENT = "31cfc67e5ef92286c3893f2a8667845cf59b544d"
SEED = 2026111202
REQUIRED = (
    "spectral_selection", "angular_symmetry", "seam", "momentum", "mouth_pair",
    "response_scaling", "independent_integrator", "coordinate_geometry",
    "one_face_susceptibility", "phi_symbol_convergence", "original_runs_tightened_energy",
    "generic_coverage", "generic_multimode_and_moving_opening",
    "generic_seam_and_momentum", "generic_energy", "four_cell_circulation",
    "four_cell_transfer", "strong_response", "strong_independent_integrator",
    "validation_negative_controls",
)


def angle_difference(a, b):
    diff = np.asarray(a)-np.asarray(b)
    diff[..., :2] = (diff[..., :2]+np.pi) % (2*np.pi)-np.pi
    return diff


def geometry_report():
    rng = np.random.default_rng(SEED)
    rows = []
    for kind, (c, r) in ant.CELLS.items():
        points = rng.uniform(.2, 5.8, (12, 3))
        numeric = np.empty((len(points), 3, 3))
        for j in range(3):
            e = np.eye(3)[j]*1e-6
            numeric[:, :, j] = angle_difference(ant.identify(points+e, kind),
                                                ant.identify(points-e, kind))/(2e-6)
        pullback, energy, canonical, involution = 0., 0., 0., 0.
        for size, n in ((16, 1), (24, 2)):
            model = ant.Throat(size, n)
            z = rng.normal(size=(3, model.dim))*.1
            d = model.deck(np.eye(model.dim), kind).T
            w = model.symplectic()
            canonical = max(canonical, float(abs(d.T@w@d-w).max()))
            involution = max(involution, float(abs(d@d-np.eye(model.dim)).max()))
            energy = max(energy, float((abs(model.energy(z@d.T)-model.energy(z))/model.energy(z)).max()))
            for eta in (-1, 1):
                x, *_ = model.unpack(z)
                dx, *_ = model.unpack(model.deck(z, kind, eta))
                for face, u in enumerate((1., -1.)):
                    for i, phi in enumerate(model.phi):
                        theta = rng.uniform(0, 2*np.pi)
                        mapped = ant.identify([phi, theta, u], kind)
                        j = int(round(mapped[0]*size/(2*np.pi))) % size
                        target = 0 if mapped[2] > 0 else 1
                        actual = eta*(x[:, 0, target, j]*np.cos(n*mapped[1])
                                      + x[:, 1, target, j]*np.sin(n*mapped[1]))
                        represented = dx[:, 0, face, i]*np.cos(n*theta)+dx[:, 1, face, i]*np.sin(n*theta)
                        pullback = max(pullback, float(abs(actual-represented).max()))
        rows.append({"cell": kind, "theta_sign": c, "normal_sign": r,
                     "measured_determinant": float(np.linalg.det(numeric[0])),
                     "jacobian_error": float(abs(numeric-ant.jacobian(kind)).max()),
                     "coordinate_square_error": float(abs(angle_difference(
                         ant.identify(ant.identify(points, kind), kind), points)).max()),
                     "state_square_error": involution, "pullback_error": pullback,
                     "energy_symmetry_error": energy, "canonical_error": canonical})
    return rows


def one_face_response_report():
    worst = 0.
    count = 0
    for size in (16, 24):
        for n in (1, 2):
            model = ant.Throat(size, n)
            k = model.frozen_matrix()
            for m in (0, 1):
                wave = np.cos(m*model.phi); wave /= np.linalg.norm(wave)
                for face in (0, 1):
                    source = np.zeros((2, size)); source[face] = wave
                    for omega in (.3, 1.4):
                        solved = np.linalg.solve(k-omega**2*np.eye(2*size), source.ravel())
                        measured = solved.reshape(2, size)@wave
                        predicted = model.susceptibility(m, omega)[:, face]
                        worst = max(worst, float(abs(measured-predicted).max()))
                        count += 1
    return {"one_face_sources": count, "max_two_face_response_error": worst}


def phi_convergence_report():
    out = []
    sizes = np.array([16, 24, 32, 48])
    for n in (1, 2):
        for m in (1, 2, 3):
            errors = [abs(ant.Throat(int(size), n).frequency_sq(m, 1)
                          -(m*m+n*n+.25+2)) for size in sizes]
            out.append({"m": m, "n": n, "errors": errors,
                        "order": float(-np.polyfit(np.log(sizes), np.log(errors), 1)[0])})
    return out


def generic_start(model, kind, eta, rng):
    z = rng.normal(size=model.dim)
    z[-4:-2] *= .03
    z[-2:] *= .07
    z = model.project(z, kind, eta)[0]
    z[:2*model.field_size] *= .4/np.linalg.norm(z[:2*model.field_size])
    return z


def trajectory_metrics(model, hist, kind, eta):
    start = hist[0]
    scale = float(np.sum(start[:2*model.field_size]**2))
    total = model.angular_momentum(hist)
    face = model.angular_momentum(hist, per_face=True)
    energy = model.energy(hist)
    x, p, _, _ = model.unpack(start)
    powers = np.sum(abs(np.fft.fft(x, axis=-1))**2
                    + abs(np.fft.fft(p, axis=-1))**2, axis=(0, 1, 2))
    return {"initial_total": float(total[0]/scale),
            "max_total": float(abs(total).max()/scale),
            "momentum_drift": float(abs(total-total[0]).max()/scale),
            "face_change": float(abs(face-face[0]).max()/scale),
            "face_max": float(abs(face).max()/scale),
            "net_flux": float(abs(model.normal_angular_flux(hist)).max()/scale),
            "seam_error": float(abs(model.deck(hist, kind, eta)-hist).max()/np.linalg.norm(start)),
            "energy_excursion": float(abs(energy-energy[0]).max()/energy[0]),
            "outside_strongest_phi_bin": float(1-powers.max()/powers.sum()),
            "initial_mouth_momentum": float(np.linalg.norm(start[-2:])),
            "initial_mouth_position": float(np.linalg.norm(start[-4:-2])),
            "max_opening": float(abs(hist[:, -4:-2]).max())}


def generic_report(size, n, draws=4, duration=40.):
    model = ant.Throat(size, n)
    rng = np.random.default_rng(SEED+size+n)
    rows = []
    for kind, (c, r) in ant.CELLS.items():
        labels = [(eta, draw) for eta in (-1, 1) for draw in range(draws)]
        starts = np.array([generic_start(model, kind, eta, rng) for eta, _ in labels])
        hist = ant.flow(model, starts, duration, h=.005, sample=20)
        per_draw = [dict(eta=eta, draw=draw,
                        **trajectory_metrics(model, hist[:, j], kind, eta))
                    for j, (eta, draw) in enumerate(labels)]
        rows.append({"sites": size, "n": n, "cell": kind,
                     "theta_sign": c, "normal_sign": r, "orientation": c*r,
                     "duration": duration, "draws": per_draw})
    return rows


def generic_checks(rows):
    all_draws = [d for row in rows for d in row["draws"]]
    covered = {(r["sites"], r["n"], r["cell"]) for r in rows}
    expected = {(s, n, c) for s, n in ((16, 1), (24, 2)) for c in ant.CELLS}
    coverage = covered == expected and len(rows) == len(expected) and all(
        {(d["eta"], d["draw"]) for d in r["draws"]}
        == {(eta, draw) for eta in (-1, 1) for draw in range(4)}
        and len(r["draws"]) == 8 and r["duration"] == 40.
        and (r["theta_sign"], r["normal_sign"]) == ant.CELLS[r["cell"]]
        and r["orientation"] == r["theta_sign"]*r["normal_sign"] for r in rows)
    circulation, transfer = True, True
    for row in rows:
        ds = row["draws"]
        if not ds:
            circulation, transfer = False, False
            continue
        if row["theta_sign"] == -1:
            circulation &= max(d["max_total"] for d in ds) < 1e-9
        else:
            circulation &= max(abs(d["initial_total"]) for d in ds) > 1e-4
        if row["orientation"] == -1:
            transfer &= max(max(d["face_change"], d["net_flux"]) for d in ds) < 1e-9
        else:
            transfer &= min(d["face_change"] for d in ds) > 1e-5
    return {"generic_coverage": coverage,
            "generic_multimode_and_moving_opening": bool(all_draws) and all(
                d["outside_strongest_phi_bin"] > .1 and d["initial_mouth_momentum"] > 1e-6
                and d["initial_mouth_position"] > 1e-6 for d in all_draws),
            "generic_seam_and_momentum": bool(all_draws) and all(
                max(d["seam_error"], d["momentum_drift"]) < 1e-9 for d in all_draws),
            "generic_energy": bool(all_draws) and all(d["energy_excursion"] < 5e-4 for d in all_draws),
            "four_cell_circulation": bool(circulation),
            "four_cell_transfer": bool(transfer)}


def strong_report():
    model = ant.Throat(16, 1)
    start = model.travelling(2.5)
    coarse = ant.flow(model, start)[:, 0]
    fine = ant.flow(model, start, h=.005, sample=20)[:, 0]
    ref = solve_ivp(lambda t, z: model.rhs(z)[0], (0., 40.), start,
                    method="DOP853", rtol=1e-10, atol=1e-12,
                    t_eval=np.linspace(0, 40, len(coarse)))
    if not ref.success or ref.y.shape[1] != len(coarse):
        raise RuntimeError(ref.message)
    scale = np.linalg.norm(ref.y[:, -1])
    ec = float(np.linalg.norm(coarse[-1]-ref.y[:, -1])/scale)
    ef = float(np.linalg.norm(fine[-1]-ref.y[:, -1])/scale)
    momentum = model.angular_momentum(coarse)
    return {**trajectory_metrics(model, coarse, "antipodal", 1),
            "amplitude": 2.5, "duration": 40.,
            "relative_momentum_drift": float(abs(momentum-momentum[0]).max()/abs(momentum[0])),
            "coarse_final_error": ec, "fine_final_error": ef,
            "improvement": ec/ef, "independent_opening_coordinates": 1}


def negative_controls(rows):
    import copy
    model = ant.Throat(16, 1)
    rng = np.random.default_rng(SEED+99)
    # Break an actual compatible state, then evolve it without projection.
    start = generic_start(model, "antipodal", -1, rng)
    start[0] += .1
    history = ant.flow(model, start, duration=.2, h=.005, sample=20)[:, 0]
    broken_metrics = trajectory_metrics(model, history, "antipodal", -1)
    broken = copy.deepcopy(rows)
    target = next(r for r in broken if r["cell"] == "antipodal" and r["sites"] == 16)
    target["draws"][0].update(broken_metrics)
    # A repeated real state has no transfer even in the allowed orientable cell.
    start = generic_start(model, "shift_kept", -1, rng)
    frozen_metrics = trajectory_metrics(model, np.repeat(start[None], 2, axis=0),
                                        "shift_kept", -1)
    frozen = copy.deepcopy(rows)
    target = next(r for r in frozen if r["cell"] == "shift_kept" and r["sites"] == 16)
    target["draws"][0].update(frozen_metrics)
    return {"broken_seam_rejected": not generic_checks(broken)["generic_seam_and_momentum"],
            "frozen_transfer_rejected": not generic_checks(frozen)["four_cell_transfer"],
            "broken_history_metrics": broken_metrics, "frozen_history_metrics": frozen_metrics}


def run_probe(progress=lambda s: None):
    base = original.run_probe(progress=progress)
    geometry = geometry_report()
    response = one_face_response_report()
    convergence = phi_convergence_report()
    rows = []
    for size, n in ((16, 1), (24, 2)):
        progress(f"generic four-cell runs: N={size}, n={n}")
        rows += generic_report(size, n)
    strong = strong_report()
    negatives = negative_controls(rows)
    checks = {k: v for k, v in base["checks"].items()
              if k not in ("geometry", "continuum_convergence", "susceptibility", "energy")}
    checks.update({
        "coordinate_geometry": all(r["jacobian_error"] < 1e-8 and
            max(r[k] for k in ("coordinate_square_error", "state_square_error", "pullback_error",
                               "energy_symmetry_error", "canonical_error")) < 1e-10 for r in geometry),
        "one_face_susceptibility": response["one_face_sources"] == 32 and response["max_two_face_response_error"] < 1e-10,
        "phi_symbol_convergence": len(convergence) == 6 and all(1.8 < r["order"] < 2.2 for r in convergence),
        "original_runs_tightened_energy": max(d["energy_error"] for r in base["runs"] for d in r["rows"]) < 2e-5,
        **generic_checks(rows),
        "strong_response": strong["max_opening"] > .5 and strong["energy_excursion"] < 5e-5
            and max(strong["seam_error"], strong["relative_momentum_drift"]) < 1e-9,
        "strong_independent_integrator": max(strong["coarse_final_error"], strong["fine_final_error"]) < 1e-2
            and strong["improvement"] > 1,
        "validation_negative_controls": negatives["broken_seam_rejected"]
            and negatives["frozen_transfer_rejected"]})
    return {"amendment": AMENDMENT, "seed": SEED,
            "prior_knowledge": "FOUR_CELL_PATTERN_AND_STRONG_RESPONSE_REPORTED_BY_REVIEWER",
            "baseline": base, "geometry": geometry, "one_face_response": response,
            "phi_symbol_convergence": convergence, "generic_runs": rows,
            "strong_response": strong, "negative_controls": negatives,
            "checks": {k: bool(v) for k, v in checks.items()},
            "verdict": ant.classify(checks, REQUIRED)}


def render(report):
    lines = ["# Four-cell review: circulation and transfer are different constraints", "",
             f"Public amendment: `{AMENDMENT}`. Pattern already supplied by review.", "",
             f"Verdict: **{report['verdict']['numerical_validity']}**.", "",
             "Scale is the initial squared field-phase-space norm (0.16 for generic runs).",
             "Each row covers eight draws; face change is measured relative to its initial value.", "",
             "| N | n | cell | max total L / scale | min face change / scale | max face change / scale | max net flux / scale |",
             "|---|---|---|---:|---:|---:|---:|"]
    for r in report["generic_runs"]:
        ds = r["draws"]
        lines.append(f"| {r['sites']} | {r['n']} | {r['cell']} | {max(d['max_total'] for d in ds):.5g} | {min(d['face_change'] for d in ds):.5g} | {max(d['face_change'] for d in ds):.5g} | {max(d['net_flux'] for d in ds):.5g} |")
    lines += ["", "| Check | Pass |", "|---|---|"]
    lines += [f"| {k} | {v} |" for k, v in report["checks"].items()]
    lines += ["", "Strong response (one independent opening coordinate):", "```json",
              json.dumps(report["strong_response"], indent=2), "```", "",
              "The orientation/transfer correspondence holds for this four-cell Hamiltonian family.",
              "It is not a general non-orientability obstruction.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    report = run_probe(progress=lambda s: print(s, flush=True))
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir/"probe.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
        (args.output_dir/"probe.md").write_text(render(report))
    print(render(report))
    return 0 if report["verdict"]["numerical_validity"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())

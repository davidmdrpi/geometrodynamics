"""What crosses an Alice handle, booked at its two mouths. Freeze ``36ba9c2``."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from geometrodynamics.bulk import alice_handle_ports as ah

PREREG = "36ba9c2"
SCOPE = ("A real scalar on the unit sphere with one handle joining antipodal "
         "mouths on the z-axis, glued torus-type (twist c) or Klein-type. "
         "Fixed geometry, separable finite-element solver per azimuthal "
         "number, evolution by mode sum; #338's mesh as a cross-check. "
         "Classical; not GR, not S3.")
SIZES = (2000, 4000)
ANGLES = (0.3, 0.5)
MS = (1, 2, 3)


def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return jsonable(x.tolist())
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    return x


def measure_angle(n, a, progress):
    out = {"P1": [], "P2": [], "P3": [], "P4": []}
    progress(f"N={n} a={a} m=0")
    r = ah.p1_m0(a, n)
    out["P1"].append({k: v for k, v in r.items()})
    loop = ah.Loop(a, n)
    for m in MS:
        progress(f"N={n} a={a} m={m}")
        sec = {name: ah.Sector(loop, m, g) for name, g in (
            ("c0", ah.TORUS), ("c1", ah.Gluing("torus", 1.0)),
            ("cpi", ah.Gluing("torus", np.pi / m)), ("klein", ah.KLEIN))}
        for name in ("c0", "c1"):
            out["P2"].append(dict(a=a, m=m, twist=name, **ah.p2_torus(sec[name])))
        out["P3"].append(dict(a=a, m=m, **ah.p3_union(
            sec["klein"], sec["c0"], sec["cpi"])))
        out["P4"].append(dict(a=a, m=m, expected_split=ah.expected_split(a, m),
                              pairs=ah.p4_p5(sec["klein"])))
        del sec
    return out


def measure_packet(n, m, progress):
    progress(f"N={n} packet m={m}")
    return ah.p6_packet(0.5, n, m)


def measure(n, progress, cache=None):
    """All of P1-P6 at ``n`` nodes.  With ``cache``, each opening angle and
    each packet is its own stage, so a run stopped partway keeps what it
    finished and the stages can be computed in separate processes."""
    out = {"P1": [], "P2": [], "P3": [], "P4": [], "P6": []}
    for a in ANGLES:
        part = cached(cache, f"measure_{n}_a{a}", lambda: measure_angle(n, a, progress))
        for key, rows in part.items():
            out[key].extend(rows)
    for m in (1, 2):
        out["P6"].append(cached(cache, f"measure_{n}_packet{m}",
                                lambda: measure_packet(n, m, progress)))
    return out


def gates(meas):
    """Pass/fail of each frozen statement at one resolution."""
    g = {}
    g["P1"] = all(r["spectrum_relative"] <= 1e-10 and r["crossing_relative"] <= 1e-10
                  for r in meas["P1"])
    g["P2 degeneracy"] = all(r["degeneracy"] <= 1e-10 for r in meas["P2"])
    g["P2 conservation"] = all(r["L_drift"] <= 1e-9 for r in meas["P2"])
    g["P3"] = all(max(r["union"], r["cos_vs_c0"], r["sin_vs_cpi"]) <= 1e-10
                  for r in meas["P3"])
    runs = [p for r in meas["P4"] for p in r["pairs"] if "L0" in p]
    g["P4 split"] = all(p["split"] > 1e-6 for r in meas["P4"] for p in r["pairs"])
    g["P4 timing"] = all(p["L0"] != 0 and p["t_found"] is not None
                         and p["timing_relative"] <= p["timing_bound"] for p in runs)
    g["P5 identity"] = (all(p["identity"] <= 1e-8 for p in runs)
                        and all(r["identity"] <= 1e-8 for r in meas["P6"]))
    g["P5 non-vacuous"] = all(p["port_fraction"] >= 0.1 for p in runs)
    g["P6 sign"] = all(r["antisymmetry"] < r["symmetry"] for r in meas["P6"])
    g["P6 tolerance"] = all(r["antisymmetry"] <= 1e-3 and r["peak_fraction"] > 1e-4
                            for r in meas["P6"])
    return g


def worst(meas):
    """The quantity behind each gate, for the shrink-at-finer-resolution rule."""
    runs = [p for r in meas["P4"] for p in r["pairs"] if "L0" in p]
    return {
        "P1": max(max(r["spectrum_relative"], r["crossing_relative"]) for r in meas["P1"]),
        "P2 degeneracy": max(r["degeneracy"] for r in meas["P2"]),
        "P2 conservation": max(r["L_drift"] for r in meas["P2"]),
        "P3": max(max(r["union"], r["cos_vs_c0"], r["sin_vs_cpi"]) for r in meas["P3"]),
        "P4 split": min(p["split"] for r in meas["P4"] for p in r["pairs"]),
        "P4 timing": max(p["timing_relative"] / p["timing_bound"] for p in runs),
        "P5 identity": max([p["identity"] for p in runs]
                           + [r["identity"] for r in meas["P6"]]),
        "P5 non-vacuous": min(p["port_fraction"] for p in runs),
        "P6 sign": max(r["antisymmetry"] / r["symmetry"] for r in meas["P6"]),
        "P6 tolerance": max(r["antisymmetry"] for r in meas["P6"]),
    }


STRUCTURAL = ("P1", "P2 degeneracy", "P3", "P4 split", "P6 sign")
HIGHER_IS_BETTER = ("P4 split", "P5 non-vacuous")


def shrinks(key, coarse, fine):
    """Whether the quantity behind gate ``key`` moved toward passing."""
    return fine > coarse if key in HIGHER_IS_BETTER else fine < coarse


def verdict(per_size, p7, registered=SIZES):
    """The freeze's verdict rules, applied only to the evidence they name.

    Returns (verdict, reasons, notes). ``COUNTEREXAMPLE`` needs a structural
    statement to fail at both registered resolutions; ``UNRESOLVED_NUMERICALLY``
    needs every miss to shrink at the finer resolution (with at most one
    further doubling of N, compared against the resolution before it). A run
    that lacks the registered resolutions or the mesh check, or that misses a
    tolerance whose miss does not shrink, gets ``NO_VERDICT``: the freeze
    names no verdict for it, and none is invented.
    """
    sizes = sorted(int(n) for n in per_size)
    ps = {int(n): v for n, v in per_size.items()}
    notes = []
    absent = [n for n in registered if n not in ps]
    if absent:
        return "NO_VERDICT", [f"registered resolutions {list(registered)} not all "
                              f"run (missing {absent})"], notes
    extra = [n for n in sizes if n not in registered]
    if extra and (len(extra) > 1 or extra[0] != 2 * max(registered)):
        return "NO_VERDICT", [f"resolutions {extra} beyond the one allowed "
                              f"doubling to {2 * max(registered)}"], notes
    if p7.get("geometry") is None:
        return "NO_VERDICT", ["P7 mesh check not run"], notes
    if not p7["geometry"]:
        return "INCOMPLETE", ["P7 geometry"], notes
    keys = list(ps[registered[0]]["gates"])
    structural = [k for k in keys if k in STRUCTURAL
                  and all(not ps[n]["gates"][k] for n in registered)]
    missed = [k for k in keys if not all(ps[n]["gates"][k] for n in sizes)]
    last, before = sizes[-1], sizes[-2]
    not_shrinking = [k for k in missed
                     if not shrinks(k, ps[before]["worst"][k], ps[last]["worst"][k])]
    if not p7["spectrum"]:
        missed.append("P7 spectrum")
    if structural:
        notes = [f"also missed: {k}" + (" (not shrinking)" if k in not_shrinking else "")
                 for k in missed if k not in structural]
        return "COUNTEREXAMPLE", structural, notes
    if not missed:
        return "CONFIRMED", [], notes
    if not_shrinking:
        return "NO_VERDICT", missed, [f"{k}: miss does not shrink from N = {before} "
                                      f"to N = {last}" for k in not_shrinking]
    return "UNRESOLVED_NUMERICALLY", missed, notes


def cached(cache, name, compute):
    """Reuse a finished stage from ``cache`` (a directory) so an interrupted
    run resumes where it stopped; the numbers are those of the stage's run."""
    if cache is None:
        return compute()
    path = Path(cache) / f"{name}.json"
    if path.exists():
        return json.loads(path.read_text())
    out = jsonable(compute())
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out))
    return out


def run_probe(progress=lambda s: None, sizes=SIZES, mesh=True, cache=None):
    started = time.time()
    per_size = {}
    for n in sizes:
        meas = cached(cache, f"measure_{n}", lambda: measure(n, progress, cache))
        per_size[n] = dict(measured=meas, gates=gates(meas), worst=worst(meas))
    p7 = {"geometry": None, "spectrum": None}
    if mesh:
        progress("P7 mesh")
        res = cached(cache, "p7_mesh", ah.p7_mesh)
        tor, kle = res["torus"], res["klein"]
        geometry = (tor["geometry"]["difference_spread"] < 1e-9
                    and kle["geometry"]["sum_spread"] < 1e-9)
        spectrum = all(r["runs"][0]["max_relative"] <= 0.02
                       and r["runs"][1]["max_relative"] < r["runs"][0]["max_relative"]
                       and r["largest_m"] < 8 for r in (tor, kle))
        p7 = dict(geometry=geometry, spectrum=spectrum, detail=res)
    v, why, notes = verdict(per_size, p7)
    return jsonable(dict(prereg=PREREG, scope=SCOPE, verdict=v, reasons=why, notes=notes,
                         sizes=list(sizes), per_size=per_size, P7=p7,
                         seconds=time.time() - started))


def markdown(rep):
    lines = [f"# Alice handle ports: probe against freeze `{rep['prereg']}`", "",
             rep["scope"], "", f"**Verdict: {rep['verdict']}**"
             + (f" ({', '.join(rep['reasons'])})" if rep["reasons"] else ""), ""]
    lines += [f"* {x}" for x in rep.get("notes", [])] + ([""] if rep.get("notes") else [])
    lines += [
             "| statement | " + " | ".join(f"N = {n}" for n in rep["sizes"]) + " |",
             "|---|" + "---|" * len(rep["sizes"])]
    first = rep["per_size"][str(rep["sizes"][0])] if str(rep["sizes"][0]) in rep["per_size"] \
        else rep["per_size"][rep["sizes"][0]]
    for key in first["gates"]:
        cells = []
        for n in rep["sizes"]:
            ps = rep["per_size"].get(str(n), rep["per_size"].get(n))
            cells.append(f"{'pass' if ps['gates'][key] else 'FAIL'} ({ps['worst'][key]:.3g})")
        lines.append(f"| {key} | " + " | ".join(cells) + " |")
    p7 = rep["P7"]
    lines += ["", f"P7 geometry: {p7['geometry']}; spectrum: {p7['spectrum']}"]
    if p7.get("detail"):
        for kind, r in p7["detail"].items():
            runs = ", ".join(f"{x['n_sphere']}/{x['n_ring']}: {x['max_relative']:.3g}"
                             for x in r["runs"])
            lines.append(f"* {kind}: largest m among the lowest 20 = {r['largest_m']}; "
                         f"max relative difference {runs}")
    lines += ["", "Lowest Klein splittings Δ/ω̄ (first resolution):", ""]
    for r in first["measured"]["P4"]:
        lines.append(f"* a = {r['a']}, m = {r['m']}: "
                     + ", ".join(f"{p['split']:.3g}" for p in r["pairs"])
                     + f"  (expectation {r['expected_split']:.2g})")
    lines += ["", f"Run time {rep['seconds']:.0f} s.", ""]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path(
        "experiments/closure_ledger/runs/20261118_alice_handle_ports"))
    parser.add_argument("--sizes", type=int, nargs="+", default=list(SIZES))
    parser.add_argument("--no-mesh", action="store_true")
    parser.add_argument("--cache", type=Path, default=None,
                        help="directory for finished stages, to resume a stopped run")
    args = parser.parse_args(argv)
    rep = run_probe(lambda s: print(s, flush=True), tuple(args.sizes),
                    mesh=not args.no_mesh, cache=args.cache)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "probe.json").write_text(json.dumps(rep, indent=1) + "\n")
    (args.out / "probe.md").write_text(markdown(rep))
    print(markdown(rep))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

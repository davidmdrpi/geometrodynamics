"""Two routes to the sheet pairing, on equal terms. Freeze ``fd78693``.

Round 14's model (``klein_handle_electric.Electric``), doubly reversed seam.
Two kinds of charge-conjugation-odd initial data drive a sheet charge:

    F   flux trapped around the handle loop, with a real field
    M   a matter current, with no trapped flux

They are compared at matched total energy. Only gauge-invariant observables
are used: sheet charge, handle-loop flux, link currents, loop holonomy.

See ``docs/klein_handle_routes_prereg.md`` and
``docs/klein_handle_routes.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from geometrodynamics.bulk import klein_handle_electric as ke
from geometrodynamics.bulk import klein_handle_response as kh

SEED = 2026111101
DRAWS = 12
ENERGIES = (2.5e-4, 1e-3, 4e-3, 1.6e-2)
REFERENCE_ENERGY = 1e-3
FLUX_FRACTION = 0.5
FLUX_FRACTIONS = (0.1, 1.0 / 3.0, 0.5, 0.9)
COUPLINGS = (0.5, 1.0, 2.0)
FLUXES = (1e-3, 2e-3, 4e-3)
STEP = kh.STEP
SHORT_STEPS = 200                    # [0, 4]
LONG_STEPS = 10_000                  # [0, 200]
EARLY_STEPS = 5                      # t = 0.1
PRIMARY, REPLICATE = kh.PRIMARY, kh.REPLICATE

INCOMPLETE = kh.INCOMPLETE
UNRESOLVED_NUMERICALLY = kh.UNRESOLVED_NUMERICALLY
COUNTEREXAMPLE = kh.COUNTEREXAMPLE


# --------------------------------------------------------------------------
# Preparations
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Draw:
    """Unit-norm components with the parities of the doubly reversed seam."""

    x: np.ndarray          # swap even
    y: np.ndarray          # swap odd
    px: np.ndarray         # swap odd,  orthogonal to y
    py: np.ndarray         # swap even, orthogonal to x


def draw(model: ke.Electric, rng: np.random.Generator) -> Draw:
    p = model.net.swap

    def part(sign):
        r = rng.standard_normal(model.n_sites)
        return 0.5 * (r + sign * p @ r)

    def unit(v):
        return v / np.linalg.norm(v)

    x, y = unit(part(+1)), unit(part(-1))
    px, py = part(-1), part(+1)
    # Total charge <x, p_y> - <y, p_x> is made zero term by term.
    px = unit(px - (px @ y) * y)
    py = unit(py - (py @ x) * x)
    return Draw(x, y, px, py)


def _state(model, x=None, y=None, px=None, py=None, flux=0.0, gauss=True):
    s = model.slices
    z = np.zeros(model.dim)
    for name, value in (("x", x), ("y", y), ("px", px), ("py", py)):
        if value is not None:
            z[s[name]] = value
    if gauss:
        rho = model.charge_density(z)[0]
        z[s["E"]] = np.linalg.lstsq(model.divergence, model.coupling * rho,
                                    rcond=None)[0]
    z[s["E"]] += flux * model.handle_cycle
    return z


def _scale_to_energy(model, build, energy: float) -> np.ndarray:
    """``build(lam)`` with ``lam`` chosen so the total energy is ``energy``."""
    def gap(lam):
        return float(model.energy(build(lam))[0]) - energy

    hi = 1.0
    while gap(hi) < 0:
        hi *= 2.0
    return build(brentq(gap, 0.0, hi, xtol=1e-15, rtol=1e-14))


def loop_edges(model: ke.Electric) -> int:
    return int(round(model.handle_cycle.sum()))


def flux_for(model: ke.Electric, energy: float, fraction: float) -> float:
    return float(np.sqrt(2.0 * fraction * energy / loop_edges(model)))


def matter(model, d: Draw, energy: float) -> np.ndarray:
    """M: complex admissible field, zero charge, loop-free Gauss field."""
    return _scale_to_energy(model, lambda lam: _state(
        model, lam * d.x, lam * d.y, lam * d.px, lam * d.py), energy)


def matter_at_rest(model, d: Draw, energy: float) -> np.ndarray:
    """M0: the review's subclass, ``Pi = a = E = 0``."""
    return _scale_to_energy(model, lambda lam: _state(
        model, lam * d.x, lam * d.y), energy)


def flux_route(model, d: Draw, energy: float, fraction: float = FLUX_FRACTION,
               at_rest: bool = False, flux: float | None = None) -> np.ndarray:
    """F (or F0): real field with energy ``(1 - f) energy`` plus loop flux."""
    px = None if at_rest else d.px
    field = _scale_to_energy(
        model, lambda lam: _state(model, lam * d.x,
                                  px=None if px is None else lam * px),
        (1.0 - fraction) * energy)
    phi = flux_for(model, energy, fraction) if flux is None else flux
    field[model.slices["E"]] += phi * model.handle_cycle
    return field


def both(model, d: Draw, energy: float, fraction: float = FLUX_FRACTION
         ) -> np.ndarray:
    """B: the matter-current field with trapped flux added."""
    z = matter(model, d, (1.0 - fraction) * energy)
    z[model.slices["E"]] += flux_for(model, energy, fraction) \
        * model.handle_cycle
    return z


def draws(model, rng, count: int = DRAWS) -> list:
    return [draw(model, rng) for _ in range(count)]


# --------------------------------------------------------------------------
# Observables
# --------------------------------------------------------------------------

def loop_current(model, z: np.ndarray) -> np.ndarray:
    return model.link_current(z) @ model.handle_cycle


def swap_even_field(model, z: np.ndarray) -> np.ndarray:
    e = model.part(z, "E")
    return 0.5 * (e + e @ model.edge_swap.T)


def peak_sheet_charge(model, starts: np.ndarray, steps: int = SHORT_STEPS
                      ) -> np.ndarray:
    """``max_t |Q_1|`` for each row of ``starts``."""
    hist = ke.flow(model, starts, steps)
    t, b = hist.shape[:2]
    return np.abs(model.sheet_charge(hist.reshape(t * b, -1))
                  ).reshape(t, b).max(axis=0)


def sheet_history(model, start: np.ndarray, steps: int = SHORT_STEPS
                  ) -> np.ndarray:
    return model.sheet_charge(ke.flow(model, start, steps)[:, 0])


def _slope(xs, ys) -> float:
    return float(np.polyfit(np.log(xs), np.log(ys), 1)[0])


def _log2_slopes(values) -> list:
    return [float(np.log2(values[i + 1] / values[i]))
            for i in range(len(values) - 1)]


# --------------------------------------------------------------------------
# P1  A flux with the field at rest
# --------------------------------------------------------------------------

def rest_flux_report(model, ds: list) -> dict:
    d0 = model.deck(-1, -1)
    s = model.slices
    # C o T_0: a -> -a, E -> E, y -> -y, p_x -> -p_x.
    sign = np.ones(model.dim)
    for name in ("y", "px", "a"):
        sign[s[name]] = -1.0
    worst, moving, flux_motion, fixed, sym = 0.0, np.inf, np.inf, 0.0, 0.0
    for energy in ENERGIES:
        rest = np.array([flux_route(model, d, energy, at_rest=True)
                         for d in ds])
        move = np.array([flux_route(model, d, energy) for d in ds])
        fixed = max(fixed, float(np.abs(rest @ d0.T - rest).max()))
        sym = max(sym, float(np.abs(rest * sign - rest).max()))
        worst = max(worst, float(peak_sheet_charge(model, rest).max()))
        moving = min(moving, float(peak_sheet_charge(model, move).min()))
        hist = ke.flow(model, rest, SHORT_STEPS)
        t, b = hist.shape[:2]
        flux = model.handle_flux(hist.reshape(t * b, -1)).reshape(t, b)
        flux_motion = min(flux_motion,
                          float(np.abs(flux - flux[0]).max(axis=0).min()))
    return {"max_sheet_charge_at_rest": worst,
            "min_sheet_charge_moving": moving,
            "min_flux_change_at_rest": flux_motion,
            "deck_fixed_error": fixed, "ct0_fixed_error": sym}


# --------------------------------------------------------------------------
# P2 - P3  Energy scaling and ranking
# --------------------------------------------------------------------------

def scaling_report(model, ds: list) -> dict:
    peaks = {"M": [], "M0": [], "F": []}
    for energy in ENERGIES:
        peaks["M"].append(peak_sheet_charge(model, np.array(
            [matter(model, d, energy) for d in ds])))
        peaks["M0"].append(peak_sheet_charge(model, np.array(
            [matter_at_rest(model, d, energy) for d in ds])))
        peaks["F"].append(peak_sheet_charge(model, np.array(
            [flux_route(model, d, energy) for d in ds])))
    peaks = {k: np.array(v) for k, v in peaks.items()}       # (energy, draw)
    out = {"energies": list(ENERGIES)}
    for name, values in peaks.items():
        med = np.median(values, axis=1)
        per_draw = [_slope(ENERGIES, values[:, j])
                    for j in range(values.shape[1])]
        out[name] = {"median_peak": med.tolist(),
                     "slope_of_median": _slope(ENERGIES, med),
                     "per_draw_slopes": per_draw}
    ratio = peaks["F"] / peaks["M"]
    med_ratio = np.median(ratio, axis=1)
    fit = np.polyfit(np.log(ENERGIES), np.log(med_ratio), 1)
    out["ratio"] = {
        "median": med_ratio.tolist(),
        "min": ratio.min(axis=1).tolist(), "max": ratio.max(axis=1).tolist(),
        "slope_of_median": float(fit[0]),
        "draws_above_one": [int(np.sum(r > 1)) for r in ratio],
        "crossover_energy": float(np.exp(-fit[1] / fit[0])),
    }
    return out


# --------------------------------------------------------------------------
# P4  Flux fraction
# --------------------------------------------------------------------------

def flux_fraction_report(model, ds: list) -> dict:
    peaks = np.array([peak_sheet_charge(model, np.array(
        [flux_route(model, d, REFERENCE_ENERGY, f) for d in ds]))
        for f in FLUX_FRACTIONS])
    reference = FLUX_FRACTIONS.index(1.0 / 3.0)
    ratios = peaks / peaks[reference]

    def predicted(f):
        return np.sqrt(3 * f) * (1 - f) * 1.5

    return {"fractions": list(FLUX_FRACTIONS),
            "median_ratio": np.median(ratios, axis=1).tolist(),
            "predicted_ratio": [float(predicted(f)) for f in FLUX_FRACTIONS],
            "best_fraction_by_draw": [float(FLUX_FRACTIONS[i]) for i in
                                      np.argmax(peaks, axis=0)]}


# --------------------------------------------------------------------------
# P5  Coupling dependence
# --------------------------------------------------------------------------

def coupling_report(net: kh.Network, ds: list) -> dict:
    base = ke.Electric(net, coupling=1.0)
    m0 = np.array([matter_at_rest(base, d, REFERENCE_ENERGY) for d in ds])
    fl = np.array([flux_route(base, d, REFERENCE_ENERGY) for d in ds])
    out = {"couplings": list(COUPLINGS), "M0": [], "F": [],
           "gauss_identity_error": 0.0}
    for g in COUPLINGS:
        model = ke.Electric(net, coupling=g)
        out["M0"].append(float(np.median(peak_sheet_charge(model, m0))))
        out["F"].append(float(np.median(peak_sheet_charge(model, fl))))
        hist = ke.flow(model, fl[:1], SHORT_STEPS)[:, 0]
        out["gauss_identity_error"] = max(
            out["gauss_identity_error"], float(np.abs(
                g * model.sheet_charge(hist) - model.seam_flux(hist)).max()))
    out["M0_slopes"] = _log2_slopes(out["M0"])
    out["F_slopes"] = _log2_slopes(out["F"])
    return out


# --------------------------------------------------------------------------
# P6  Electric sectors
# --------------------------------------------------------------------------

def sector_report(model, ds: list) -> dict:
    out = {}
    for name, build in (("M0", matter_at_rest), ("M", matter),
                        ("F", flux_route)):
        starts = np.array([build(model, d, REFERENCE_ENERGY) for d in ds])
        hist = ke.flow(model, starts, SHORT_STEPS)
        flat = hist.reshape(-1, model.dim)
        e = model.part(flat, "E")
        even = swap_even_field(model, flat)
        out[name] = {
            "handle_flux": float(np.abs(model.handle_flux(flat)).max()),
            "loop_holonomy": float(np.abs(model.loop_holonomy(flat)).max()),
            "loop_current": float(np.abs(loop_current(model, flat)).max()),
            "swap_even_field": float(np.abs(even).max()),
            "swap_odd_field": float(np.abs(e - even).max()),
            "sheet_charge": float(np.abs(model.sheet_charge(flat)).max()),
        }
    return out


# --------------------------------------------------------------------------
# P7  Interference
# --------------------------------------------------------------------------

def interference_report(model, ds: list) -> dict:
    s = model.slices
    cycle = model.handle_cycle
    defects, differences, symmetry = [], [], 0.0
    scale = 0.0
    for phi in FLUXES:
        defect, diff = [], []
        for d in ds:
            m = matter(model, d, REFERENCE_ENERGY)

            def with_flux(z, sign=1.0):
                out = z.copy()
                out[s["E"]] += sign * phi * cycle
                return out

            real = np.zeros(model.dim)
            imag = np.zeros(model.dim)
            for name in ("x", "px"):
                real[s[name]] = m[s[name]]
            for name in ("y", "py"):
                imag[s[name]] = m[s[name]]
            starts = np.array([m, with_flux(m), with_flux(m, -1.0),
                               with_flux(real), with_flux(imag),
                               model.conjugation @ with_flux(m)])
            hist = ke.flow(model, starts, SHORT_STEPS)
            t, b = hist.shape[:2]
            q = model.sheet_charge(hist.reshape(t * b, -1)).reshape(t, b)
            q_m, q_b, q_bm, q_fx, q_fy, q_conj = q.T
            defect.append(float(np.abs(q_b - q_m - q_fx - q_fy).max()))
            diff.append(float(np.abs(q_b - q_bm).max()))
            symmetry = max(symmetry, float(np.abs(q_conj + q_b).max()))
            scale = max(scale, float(np.abs(q_m).max()))
        defects.append(float(np.median(defect)))
        differences.append(float(np.median(diff)))
    return {"fluxes": list(FLUXES), "median_defect": defects,
            "defect_slopes": _log2_slopes(defects),
            "median_flux_reversal_difference": differences,
            "reversal_slopes": _log2_slopes(differences),
            "joint_reversal_error": symmetry, "matter_scale": scale}


# --------------------------------------------------------------------------
# P8  Screening rule
# --------------------------------------------------------------------------

def screening_report(model, ds: list) -> dict:
    s, g = model.slices, model.coupling
    cycle = model.handle_cycle
    edges = np.flatnonzero(cycle > 0.5)
    n_loop = loop_edges(model)
    t1 = EARLY_STEPS * STEP
    rows = []
    for d in ds:
        z = flux_route(model, d, REFERENCE_ENERGY)
        x = z[s["x"]]
        bond = float(np.sum(x[model.tail[edges]] * x[model.head[edges]]))
        phi = float(model.handle_flux(z)[0])
        hist = ke.flow(model, z, SHORT_STEPS)[:, 0]
        flux = model.handle_flux(hist)
        early = float(flux[EARLY_STEPS] - phi)
        predicted = -0.5 * g ** 2 * phi * bond / n_loop * t1 ** 2
        late = flux - phi
        peak = float(late[np.argmax(np.abs(late))])
        rows.append({
            "bond_order": bond / float(x @ x),
            "measured": early, "predicted": predicted,
            "ratio": early / predicted,
            "sign_agrees": bool(early * predicted > 0),
            "resolved": bool(abs(bond) > 0.02 * float(x @ x)),
            "peak_change_sign_matches_early": bool(peak * early > 0),
        })
    resolved = [r for r in rows if r["resolved"]]
    return {
        "time": t1, "rows": rows, "resolved_draws": len(resolved),
        "worst_ratio_error": float(max(abs(r["ratio"] - 1) for r in resolved))
        if resolved else float("nan"),
        "sign_agreement": int(sum(r["sign_agrees"] for r in resolved)),
        "positive_bond_orders": int(sum(r["bond_order"] > 0 for r in rows)),
        "peak_sign_matches_early": int(sum(
            r["peak_change_sign_matches_early"] for r in rows)),
    }


def screening_convergence_report(model, ds: list,
                                 step_counts=(2, 5, 10)) -> dict:
    """Post-freeze. Is the frozen formula the leading order?

    The frozen 5 % tolerance at ``t = 0.1`` was missed. The first correction
    I blamed, the rate of the bond order, vanishes identically: ``x`` is swap
    even, ``p_x`` swap odd and the handle loop swap invariant. The next one
    is second order in ``t``. If the formula is the leading order, its
    relative error must fall as ``t^2``.
    """
    s, g = model.slices, model.coupling
    edges = np.flatnonzero(model.handle_cycle > 0.5)
    n_loop = loop_edges(model)
    errors, rates, signs = [], [], []
    for d in ds:
        z = flux_route(model, d, REFERENCE_ENERGY)
        x, px = z[s["x"]], z[s["px"]]
        tail, head = model.tail[edges], model.head[edges]
        bond = float(np.sum(x[tail] * x[head]))
        if abs(bond) <= 0.02 * float(x @ x):
            continue
        rates.append(abs(float(np.sum(px[tail] * x[head]
                                      + x[tail] * px[head]))) / abs(bond))
        phi = float(model.handle_flux(z)[0])
        flux = model.handle_flux(ke.flow(model, z, max(step_counts))[:, 0])
        row = []
        for n in step_counts:
            t = n * STEP
            predicted = -0.5 * g ** 2 * phi * bond / n_loop * t ** 2
            row.append((flux[n] - phi) / predicted - 1.0)
        errors.append(row)
        signs.append(all(abs(e) < 1 for e in row))
    errors = np.abs(np.array(errors))
    median = np.median(errors, axis=0)
    times = [n * STEP for n in step_counts]
    return {"times": times, "resolved_draws": int(errors.shape[0]),
            "median_error": median.tolist(),
            "max_error": errors.max(axis=0).tolist(),
            "error_slope": _slope(times, median),
            "bond_rate_over_bond": float(max(rates)),
            "sign_agreement_all_times": int(sum(signs))}


# --------------------------------------------------------------------------
# P9  Which reversal each route needs
# --------------------------------------------------------------------------

def reversal_report(model, ds: list, rng: np.random.Generator) -> dict:
    alice = model.deck(1, -1)            # time kept, fibre reversed
    klein = model.deck(-1, -1)
    starts = np.array([matter_at_rest(model, d, REFERENCE_ENERGY)
                       for d in ds])
    hist = ke.flow(model, starts, SHORT_STEPS)
    flat = hist.reshape(-1, model.dim)
    general = np.array([matter(model, d, REFERENCE_ENERGY) for d in ds])

    # Ordinary handle: everything swap even, zero charge.
    p = model.net.swap
    ordinary = []
    for flux in (0.0, FLUXES[-1]):
        for _ in range(len(ds)):
            parts = [0.5 * (r + p @ r) for r in
                     rng.standard_normal((4, model.n_sites))]
            x, y, px, py = parts
            px = px - (px @ y) / (y @ y) * y
            py = py - (py @ x) / (x @ x) * x
            z = _state(model, x, y, px, py, flux=flux)
            z[:4 * model.n_sites] *= 0.1 / np.linalg.norm(
                z[:4 * model.n_sites])
            z[model.slices["E"]] = np.linalg.lstsq(
                model.divergence, model.coupling * model.charge_density(z)[0],
                rcond=None)[0] + flux * model.handle_cycle
            ordinary.append(z)
    ordinary = np.array(ordinary)
    plain = model.deck(1, 1)
    ohist = ke.flow(model, ordinary, SHORT_STEPS).reshape(-1, model.dim)
    return {
        "rest_data_fixed_by_alice_deck": float(np.abs(
            starts @ alice.T - starts).max()),
        "rest_data_fixed_by_klein_deck": float(np.abs(
            starts @ klein.T - starts).max()),
        "alice_fixed_set_invariance": float(np.abs(
            flat @ alice.T - flat).max()),
        "rest_sheet_charge": float(np.abs(model.sheet_charge(flat)).max()),
        "general_matter_not_alice_fixed": float(np.abs(
            general @ alice.T - general).max()),
        "ordinary_fixed_error": float(np.abs(
            ordinary @ plain.T - ordinary).max()),
        "ordinary_gauss": float(np.abs(model.gauss(ordinary)).max()),
        "ordinary_sheet_charge": float(np.abs(
            model.sheet_charge(ohist)).max()),
        "ordinary_local_charge": float(np.abs(
            model.charge_density(ohist)).max()),
    }


# --------------------------------------------------------------------------
# P10 and the declared measurements: the long window
# --------------------------------------------------------------------------

def long_observables(model, starts: np.ndarray, steps: int = LONG_STEPS,
                     h: float = STEP) -> dict:
    """Sheet charge, handle flux and link-current overlap, step by step.

    Histories are not stored; the long window would not fit in memory.
    """
    z = np.atleast_2d(np.asarray(starts, dtype=float))
    j0 = model.link_current(z)
    norm = np.sum(j0 ** 2, axis=1)
    norm = np.where(norm > 0, norm, 1.0)
    q = np.empty((steps + 1, z.shape[0]))
    flux = np.empty_like(q)
    overlap = np.empty_like(q)
    for n in range(steps + 1):
        q[n] = model.sheet_charge(z)
        flux[n] = model.handle_flux(z)
        overlap[n] = np.sum(model.link_current(z) * j0, axis=1) / norm
        if n < steps:
            z = kh.midpoint_step(model, z, h)
    return {"sheet": q, "flux": flux, "overlap": overlap,
            "final_gauss": float(np.abs(model.gauss(z)).max())}


def _stats(q: np.ndarray) -> dict:
    rms = np.sqrt(np.mean(q ** 2, axis=0))
    mean = np.mean(q, axis=0)
    return {"rms": rms, "dc": np.abs(mean) / rms}


def long_window_report(model, ds: list, steps: int = LONG_STEPS) -> dict:
    out = {"window": steps * STEP, "energies": list(ENERGIES), "rows": []}
    short = min(SHORT_STEPS, steps)
    for energy in ENERGIES:
        builders = {"M": matter, "M0": matter_at_rest, "F": flux_route,
                    "B": both}
        obs = {name: long_observables(model, np.array(
            [fn(model, d, energy) for d in ds]), steps)
            for name, fn in builders.items()}
        st = {name: _stats(o["sheet"]) for name, o in obs.items()}
        short_rms = {name: _stats(o["sheet"][:short + 1])["rms"]
                     for name, o in obs.items()}
        peak = {name: np.abs(o["sheet"][:short + 1]).max(axis=0)
                for name, o in obs.items()}
        flux = obs["F"]["flux"]
        retention = flux / flux[0]
        overlap = obs["M0"]["overlap"]
        below = overlap < 0.5
        half = np.where(below.any(axis=0), below.argmax(axis=0), steps) * STEP
        out["rows"].append({
            "energy": energy,
            "median_rms": {k: float(np.median(v["rms"]))
                           for k, v in st.items()},
            "median_dc_fraction": {k: float(np.median(v["dc"]))
                                   for k, v in st.items()},
            "ratio_long_rms": float(np.median(st["F"]["rms"]
                                              / st["M"]["rms"])),
            "ratio_short_peak": float(np.median(peak["F"] / peak["M"])),
            "ratio_short_rms": float(np.median(short_rms["F"]
                                               / short_rms["M"])),
            "draws_flux_route_larger_long": int(np.sum(
                st["F"]["rms"] > st["M"]["rms"])),
            "flux_retention_min": float(np.median(retention.min(axis=0))),
            "flux_retention_final": float(np.median(retention[-1])),
            "draws_flux_changes_sign": int(np.sum(retention.min(axis=0) < 0)),
            "median_current_half_life": float(np.median(half)),
            "final_gauss": max(o["final_gauss"] for o in obs.values()),
        })
    return out


def numerical_report(model, d: Draw, steps: int = LONG_STEPS) -> dict:
    """Halved step and an independent integrator, one draw, long window."""
    from scipy.integrate import solve_ivp

    out = {}
    for name, fn in (("M", matter), ("F", flux_route)):
        z0 = fn(model, d, REFERENCE_ENERGY)
        coarse = long_observables(model, z0, steps)["sheet"][:, 0]
        fine = long_observables(model, z0, 2 * steps, STEP / 2)["sheet"][:, 0]
        ref = solve_ivp(lambda t, y: model.vector_field(y)[0],
                        (0.0, steps * STEP), z0, method="DOP853",
                        rtol=1e-10, atol=1e-13,
                        t_eval=np.linspace(0.0, steps * STEP, steps + 1))
        ref_q = model.sheet_charge(ref.y.T)
        rms = lambda v: float(np.sqrt(np.mean(v ** 2)))
        out[name] = {
            "rms": rms(coarse), "rms_halved_step": rms(fine),
            "rms_dop853": rms(ref_q),
            "rms_change_halved_step": abs(rms(fine) - rms(coarse)) / rms(fine),
            "rms_change_dop853": abs(rms(ref_q) - rms(fine)) / rms(ref_q),
            "pointwise_dop853_vs_halved": float(
                np.abs(ref_q - fine[::2]).max() / np.abs(ref_q).max()),
            "dop853_gauss": float(np.abs(model.gauss(ref.y.T)).max()),
        }
    return out


# --------------------------------------------------------------------------
# Verdict
# --------------------------------------------------------------------------

VERDICT_FIELDS = ("numerical_validity", "flux_route_requirement",
                  "energy_scaling", "short_time_ranking", "flux_fraction",
                  "coupling_dependence", "electric_sectors", "interference",
                  "screening_rule", "screening_leading_order",
                  "reversal_requirements", "long_time_ranking",
                  "long_window_numerics", "support_dynamics",
                  "orienting_datum", "positive_response", "event_law")
CONFIRMED = {
    "flux_route_requirement": "FLUX_ROUTE_NEEDS_FIELD_MOMENTUM",
    "energy_scaling": "MATTER_LINEAR_FLUX_THREE_HALVES",
    "short_time_ranking": "MATTER_CURRENT_DOMINATES_BY_INVERSE_SQRT_ENERGY",
    "flux_fraction": "OPTIMUM_AT_ONE_THIRD_IN_THE_FLUX",
    "coupling_dependence": "MATTER_ROUTE_IS_TRANSPORT_FLUX_ROUTE_IS_ELECTRIC",
    "electric_sectors": "ROUTES_OCCUPY_DIFFERENT_ELECTRIC_SECTORS",
    "interference": "NO_FIRST_ORDER_INTERFERENCE_RELATIVE_SIGN_PHYSICAL",
    "screening_rule": "EARLY_FLUX_CHANGE_SET_BY_LOOP_BOND_ORDER",
    # Added after the frozen tolerance was missed; not a frozen prediction.
    "screening_leading_order":
        "BOND_ORDER_FIXES_THE_SIGN_AND_LEADING_ORDER_OF_EARLY_FLUX_CHANGE",
    "long_window_numerics": "STATISTICS_STABLE_UNDER_STEP_AND_INTEGRATOR",
    "reversal_requirements": "MATTER_ROUTE_NEEDS_ONLY_FIBRE_REVERSAL",
    "long_time_ranking": "FLUX_ROUTE_GAINS_AT_LONG_TIMES",
}
FIXED = {
    "support_dynamics": "EFFECTIVE_SUPPORT_ASSUMED",
    "orienting_datum": "NOT_DERIVED",
    "positive_response": "NOT_DERIVED",
    "event_law": "NOT_DERIVED",
}


def verdict(checks: dict, groups: dict, counterexamples: dict | None = None
            ) -> dict:
    """Status per field.

    A check that misses its tolerance makes its field
    ``UNRESOLVED_NUMERICALLY``. A frozen prediction contradicted by a
    resolved measurement is registered in ``counterexamples`` and makes its
    field ``COUNTEREXAMPLE``; that is a statement about the prediction, not
    about the numerics, so it does not affect ``numerical_validity``.
    """
    counterexamples = counterexamples or {}
    out, failed, missing_all, contradicted = {}, [], [], []
    for name in CONFIRMED:
        names = groups.get(name, [])
        missing = [c for c in names if c not in checks]
        bad = [c for c in names if c in checks and not checks[c]]
        if name in counterexamples:
            out[name] = COUNTEREXAMPLE
            contradicted += bad
        elif missing or not names:
            out[name] = INCOMPLETE
            missing_all += missing or [f"<no checks for {name}>"]
        elif bad:
            out[name] = UNRESOLVED_NUMERICALLY
            failed += bad
        else:
            out[name] = CONFIRMED[name]
    out.update(FIXED)
    out["numerical_validity"] = (
        INCOMPLETE if missing_all
        else UNRESOLVED_NUMERICALLY if failed else "PASS")
    out["failed_checks"] = failed
    out["missing_checks"] = missing_all
    out["contradicted_checks"] = contradicted
    out["counterexamples"] = dict(counterexamples)
    return out

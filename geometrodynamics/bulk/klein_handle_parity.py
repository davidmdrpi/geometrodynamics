"""Temporal symmetry of exchange across the seam. Freeze ``95d7a80``.

For a symmetry generator ``X`` with sheet Noether value ``N_1``:

    seam reverses time, X anticommutes with the deck:  N_1(t) + N_1(-t) = N_tot   (odd)
    seam reverses time, X commutes:                    N_tot = 0, N_1 even
    seam keeps time,    X commutes:                    N_1 = N_tot / 2            (constant)
    seam keeps time,    X anticommutes:                N_tot = 0, N_1 unconstrained

Tested on round 15's ungauged doublet (``klein_handle_matching.Multiplet``).
See ``docs/klein_handle_parity_prereg.md`` and ``docs/klein_handle_parity.md``.
"""

from __future__ import annotations

import numpy as np

from geometrodynamics.bulk import klein_handle_matching as km
from geometrodynamics.bulk import klein_handle_response as kh

SEED = 2026111102
DRAWS = 6
STEPS, STEP = kh.STEPS, kh.STEP
GENERATORS = ("i", "sigma_x", "sigma_y", "sigma_z")
SEAMS = (("J", -1), ("identity", 1))
PRIMARY, REPLICATE = kh.PRIMARY, kh.REPLICATE

# Predicted class for (internal map, s_t, generator kind).
PREDICTED = {
    ("J", -1, "charge"): "odd", ("J", -1, "spin"): "even",
    ("J", 1, "charge"): "unconstrained", ("J", 1, "spin"): "constant",
    ("identity", -1, "charge"): "even", ("identity", -1, "spin"): "even",
    ("identity", 1, "charge"): "constant", ("identity", 1, "spin"): "constant",
}


def sheet_noether(model: km.Multiplet, name: str, z: np.ndarray) -> np.ndarray:
    """Sheet-1 part of ``N_X = 1/2 Omega(z, X z) = -Re sum conj(Pi) X Psi``."""
    psi, mom, *_ = model.unpack(z)
    n = model.net.n
    acted = model.internal(name)(psi)
    return -np.sum(np.real(np.conj(mom[:, :, :n]) * acted[:, :, :n]),
                   axis=(1, 2))


def admissible(model: km.Multiplet, deck: np.ndarray, s_t: int,
               rng: np.random.Generator, amplitude: float = 0.125
               ) -> np.ndarray:
    f = model.field_dim
    basis = kh.null_space(deck[:f, :f] - np.eye(f))
    z = np.zeros(model.dim)
    z[:f] = basis @ rng.standard_normal(basis.shape[1])
    z[:f] *= amplitude / np.linalg.norm(z[:f])
    z[f:f + 2] = 0.02
    z[f + 2], z[f + 3] = 0.01, s_t * 0.01
    return z


def classify(n_f: np.ndarray, n_b: np.ndarray, total: float, scale: float
             ) -> dict:
    """Even and odd parts in time of the sheet value about ``total / 2``."""
    dev_f, dev_b = n_f - total / 2, n_b - total / 2
    even, odd = 0.5 * (dev_f + dev_b), 0.5 * (dev_f - dev_b)
    return {"even_part": float(np.abs(even).max() / scale),
            "odd_part": float(np.abs(odd).max() / scale),
            "variation": float((n_f.max() - n_f.min()) / scale),
            "total": float(total / scale)}


def parity_report(net: kh.Network, rng: np.random.Generator,
                  draws: int = DRAWS) -> list:
    rows = []
    for internal, loop in SEAMS:
        model = km.Multiplet(net, components=2, loop_sign=loop)
        for s_t in (-1, 1):
            deck = model.deck(internal, s_t)
            per_gen = {g: [] for g in GENERATORS}
            defect = 0.0
            for _ in range(draws):
                z0 = admissible(model, deck, s_t, rng)
                fwd = km.flow(model, z0, STEPS)[:, 0]
                bwd = km.flow(model, z0, STEPS, -STEP)[:, 0]
                scale = float(np.sum(z0[:model.field_dim] ** 2))
                if s_t == -1:
                    defect = max(defect, float(np.abs(fwd - bwd @ deck.T).max()
                                               / np.linalg.norm(z0)))
                else:
                    defect = max(defect, float(np.abs(fwd @ deck.T - fwd).max()
                                               / np.linalg.norm(z0)))
                for g in GENERATORS:
                    total = float(model.noether(g, z0[None])[0])
                    drift = float(np.abs(model.noether(g, fwd) - total).max()
                                  / scale)
                    c = classify(sheet_noether(model, g, fwd),
                                 sheet_noether(model, g, bwd), total, scale)
                    c["drift"] = drift
                    per_gen[g].append(c)
            for g in GENERATORS:
                kind = "charge" if g == "i" else "spin"
                cs = per_gen[g]
                rows.append({
                    "internal": internal, "loop_sign": loop, "s_t": s_t,
                    "generator": g, "kind": kind,
                    "predicted": PREDICTED[(internal, s_t, kind)],
                    "max_even_part": max(c["even_part"] for c in cs),
                    "min_even_part": min(c["even_part"] for c in cs),
                    "max_odd_part": max(c["odd_part"] for c in cs),
                    "min_odd_part": min(c["odd_part"] for c in cs),
                    "max_variation": max(c["variation"] for c in cs),
                    "min_variation": min(c["variation"] for c in cs),
                    "max_abs_total": max(abs(c["total"]) for c in cs),
                    "max_drift": max(c["drift"] for c in cs),
                    "symmetry_defect": defect,
                })
    return rows


def satisfies(row: dict, zero: float = 1e-9, present: float = 1e-6) -> dict:
    """S1 (the relation) and S2 (not vacuous) for one row."""
    p = row["predicted"]
    if p == "odd":
        rule = row["max_even_part"] < zero
        live = (row["min_odd_part"] > present
                and row["min_variation"] > present)
    elif p == "even":
        rule = row["max_odd_part"] < zero and row["max_abs_total"] < zero
        # A nonzero even part is not exchange: a constant sheet value has
        # one. Exchange means the value varies (corrected after review).
        live = (row["min_even_part"] > present
                and row["min_variation"] > present)
    elif p == "constant":
        rule = row["max_variation"] < zero
        live = True
    else:                                   # unconstrained
        rule = row["max_abs_total"] < zero
        live = (row["min_even_part"] > 1e-7 and row["min_odd_part"] > 1e-7
                and row["min_variation"] > present)
    return {"rule": bool(rule), "live": bool(live)}


VERDICT_FIELDS = ("numerical_validity", "parity_rule", "non_vacuous",
                  "charge_against_spin", "selection", "bam_cell",
                  "positive_response", "event_law")
CONFIRMED = {
    "parity_rule": "EXCHANGE_PARITY_FIXED_BY_TIME_REVERSAL_AND_COMMUTATION",
    "non_vacuous": "ALLOWED_EXCHANGE_OCCURS_EXCLUDED_PARITY_ABSENT",
    "charge_against_spin":
        "OPPOSITE_PARITY_ON_THE_TRANSPORT_SEAM_WITH_TIME_REVERSAL",
}
FIXED = {
    "selection": "BY_FIELD_TYPE_AND_TIME_EXTENSION_NOT_BY_GEOMETRY_ALONE",
    "bam_cell": "DEPENDS_ON_FIELD_TYPE_AND_TIME_EXTENSION",
    "positive_response": "NOT_DERIVED",
    "event_law": "NOT_DERIVED",
}


def verdict(checks: dict, groups: dict) -> dict:
    out, failed, missing_all = {}, [], []
    for name in CONFIRMED:
        names = groups.get(name, [])
        missing = [c for c in names if c not in checks]
        bad = [c for c in names if c in checks and not checks[c]]
        if missing or not names:
            out[name] = kh.INCOMPLETE
            missing_all += missing or [f"<no checks for {name}>"]
        elif bad:
            out[name] = kh.UNRESOLVED_NUMERICALLY
            failed += bad
        else:
            out[name] = CONFIRMED[name]
    out.update(FIXED)
    out["numerical_validity"] = (
        kh.INCOMPLETE if missing_all
        else kh.UNRESOLVED_NUMERICALLY if failed else "PASS")
    out["failed_checks"], out["missing_checks"] = failed, missing_all
    return out

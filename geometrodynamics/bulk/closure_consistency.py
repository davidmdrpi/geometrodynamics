"""Independent overlaps of rounds 5--8 and the classical equilibrium model.

All comparisons use the chosen geodesic triangle family, not full BAM
history space. Half-phase, singlet partner sign, and priors are explicit.
"""

import math

import numpy as np

from geometrodynamics.bulk import closure_equilibrium as equilibrium
from geometrodynamics.bulk import closure_measurement as measurement
from geometrodynamics.bulk import closure_current as current
from geometrodynamics.bulk import history_action as action
from geometrodynamics.bulk.conditioning_variable import n_window_is_sector_blind
from geometrodynamics.bulk.mouth_spin_frame import qmul
from geometrodynamics.bulk.source_readout import source_circle


def run_consistency(angles=(.3, .7, 1., 1.4, 2., 2.6), n=32768):
    errors = dict.fromkeys(("quaternion_phase", "holonomy_trace", "positive_mass",
                           "oriented_mass", "positive_joint", "oriented_joint",
                           "phase_gradient", "frame_hessian", "source_density",
                           "source_parity", "source_normalization", "sector_prior",
                           "normal_model_E", "normal_window_E"), 0.)
    rows = []

    def record(key, error):
        errors[key] = max(errors[key], float(np.max(np.abs(error))))

    rng = np.random.default_rng(20260905)
    for gamma in angles:
        a = np.array([0., 0., 1.])
        b = np.array([math.sin(gamma), 0., math.cos(gamma)])
        x, weights = source_circle(a, b, n)
        record("source_density", weights - measurement.source_density_on_circle(a, b, x) / n)
        record("source_parity", weights - np.roll(weights, n // 2))
        record("source_normalization", weights.sum() - 1)
        positives, oriented = {}, {}
        for sa in (-1, 1):
            for sb in (-1, 1):
                u, w = sa * a, -sb * b
                c, cross = float(u @ w), float(np.linalg.norm(np.cross(u, w)))
                masses = action.morse_bott_component_masses(u, w, n=n)
                pos, ori = masses["M_0"] + masses["M_pi"], masses["M_0"] - masses["M_pi"]
                positives[sa, sb], oriented[sa, sb] = pos, ori
                record("positive_mass", pos - equilibrium.limiting_mass(c))
                record("oriented_mass", ori - 2 * math.pi * (1 + c) / cross)
                rows.append({"gamma": gamma, "s_A": sa, "s_B": sb, "c": c,
                             "M_0": masses["M_0"], "M_pi": masses["M_pi"],
                             "equilibrium_mass": equilibrium.limiting_mass(c)})
                # Direct quaternion multiplication, independently of N/D.
                for point in rng.normal(size=(4, 3)):
                    point /= np.linalg.norm(point)
                    lift = current.minimal_rotation_lift
                    G = qmul(qmul(lift(w, point), lift(u, w)), lift(point, u))
                    N, D = equilibrium.triangle_data(point, u, w)
                    expected = np.r_[D, N * point] / math.hypot(N, D)
                    record("quaternion_phase", G - expected)
                    record("holonomy_trace", action.holonomy_trace(point, u, w) + G[0])
                normal = np.cross(u, w) / cross
                for point in x[::max(1, n // 12)]:
                    _, D = equilibrium.triangle_data(point, u, w)
                    if abs(D) < .15:
                        continue  # singular punctures are outside this regular test
                    h = 2e-6
                    plus = math.cos(h) * point + math.sin(h) * normal
                    minus = math.cos(h) * point - math.sin(h) * normal
                    phases = measurement.solid_angle(np.array([plus, minus]), u, w) / 2
                    delta = (phases[0] - phases[1] + math.pi) % (2 * math.pi) - math.pi
                    record("phase_gradient", delta / (2 * h) - cross / D)
                    hessian = (equilibrium.rotor_potential(plus, u, w)
                               + equilibrium.rotor_potential(minus, u, w)) / h ** 2
                    record("frame_hessian", hessian - (cross / D) ** 2)
        eq_joint = equilibrium.joint_probabilities(gamma, None)
        old_joint = measurement.closed_form_weights(math.pi - gamma)  # partner sign
        old_oriented = current.singlet_loop_law(gamma, n=n)["P"]
        for sector in positives:
            record("positive_joint", positives[sector] / sum(positives.values()) - eq_joint[sector])
            record("positive_joint", old_joint[sector] - eq_joint[sector])
            record("oriented_joint", oriented[sector] / sum(oriented.values()) - old_oriented[sector])
        for row in current.sector_prior_control(math.pi - gamma)["rows"]:
            E = equilibrium.correlation(equilibrium.joint_probabilities(
                gamma, None, prior_ratio=row["ratio"]))
            record("sector_prior", E - row["E"])
        record("normal_model_E", equilibrium.correlation(equilibrium.joint_probabilities(
            gamma, 64., model="normal")))
        record("normal_window_E", n_window_is_sector_blind(gamma, eps=.05, n=10000)["E"])
    return {"max_residuals": errors, "mass_rows": rows, "n_circle": n,
            "angle_count": len(angles), "sector_count": len(rows)}

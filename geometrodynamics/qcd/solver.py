"""
Störmer–Verlet solver for the hadronic network.

Implements the wave equation on branched 1D domains with:
- SAT (Simultaneous Approximation Term) junction boundary treatment
- Coupled nonlinear junction solve via scipy.optimize.root
- Loop-current evolution for hybrid excitations
- Bridge field nucleation and pair creation
"""

from __future__ import annotations

from math import pi, sqrt, isfinite
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import root as scipy_root

# NumPy 2.0 renamed np.trapz → np.trapezoid; support both.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz

from geometrodynamics.qcd.constants import (
    SAT_PENALTY,
    SAT_DAMPING,
    KAPPA_J,
    HBAR_C,
    SIGMA_FM,
    A_COULOMB,
)
from geometrodynamics.qcd.network import HadronicNetwork
from geometrodynamics.qcd.bridge import cornell_equilibrium_amplitude


class HadronicNetworkSolver:
    """Coupled wave-equation solver on a hadronic flux-tube network."""

    def __init__(self, net: HadronicNetwork, antipodal_coupling: float = 0.15):
        self.net = net
        self.dt = net.dt
        self.N = net.N_grid
        self.alpha_ap = antipodal_coupling
        self._init_bufs()
        self.time = 0.0
        self.E_hist: List[float] = []
        self.t_hist: List[float] = []
        self.nucleation_events: List[Dict] = []
        self._junction_u: Dict[int, float] = {}
        self._pair_creation_penalty: float = 0.0
        self._junction_state: Dict[int, Dict[int, Dict[str, float]]] = {}
        self._junction_flux: Dict[int, Tuple[Optional[float], Optional[float]]] = {}
        self._junction_sat: Dict[int, Tuple[Optional[float], Optional[float]]] = {}

    def _init_bufs(self) -> None:
        self._buf: Dict[int, np.ndarray] = {}
        self._ptr: Dict[int, int] = {}
        for bid, b in self.net.branches.items():
            n = max(2, int(round(b.tau_delay / self.dt)))
            self._buf[bid] = np.zeros((n, self.N))
            self._ptr[bid] = 0

    def _antipodal_delayed(self, bid: int, psi: np.ndarray) -> np.ndarray:
        b = self.net.branches[bid]
        buf, ptr = self._buf[bid], self._ptr[bid]
        old = buf[ptr].copy()
        if b.attached_ring or b.periodic:
            buf[ptr] = np.roll(psi, self.N // 2)
        elif b.mobius:
            buf[ptr] = -psi[::-1]
        else:
            buf[ptr] = psi[::-1]
        self._ptr[bid] = (ptr + 1) % len(buf)
        return old

    def _is_ep(self, nid: int) -> bool:
        return self.net.nodes[nid].endpoint

    def _branch_acc(self, bid: int, psi: np.ndarray) -> np.ndarray:
        b = self.net.branches[bid]
        N = self.N
        ds = b.length / (N - 1)
        v2 = b.v_speed ** 2
        s = np.linspace(0, b.length, N)
        U = np.array([b.U_potential(si, b.length) for si in s])
        mu2 = self.net.confinement_mass_sq(bid)
        lap = np.empty(N)
        lap[1:-1] = (psi[:-2] - 2 * psi[1:-1] + psi[2:]) / ds ** 2

        # Left end
        if b.periodic:
            lap[0] = (psi[-1] - 2 * psi[0] + psi[1]) / ds ** 2
        elif b.attached_ring and b.attach_end == "left":
            lap[0] = (-2 * psi[0] + psi[1]) / ds ** 2
        elif self._is_ep(b.node_left):
            lap[0] = 0.0
        else:
            lap[0] = (-2 * psi[0] + 2 * psi[1]) / ds ** 2

        # Right end
        if b.periodic:
            lap[-1] = (psi[-2] - 2 * psi[-1] + psi[0]) / ds ** 2
        elif b.attached_ring and b.attach_end == "left":
            lap[-1] = (psi[-2] - 2 * psi[-1] + psi[0]) / ds ** 2
        elif b.mobius:
            lap[-1] = (psi[-2] - 2 * psi[-1] + psi[-2]) / ds ** 2
        elif self._is_ep(b.node_right):
            lap[-1] = 0.0
        else:
            lap[-1] = (2 * psi[-2] - 2 * psi[-1]) / ds ** 2

        anti = self._antipodal_delayed(bid, psi)
        acc = (
            v2 * lap
            - U * psi
            - b.lambda_nl * psi ** 3
            + self.alpha_ap * float(b.orientation) * anti
            - mu2 * psi
        )

        # SAT penalty + velocity damping
        pen_left, pen_right = self._junction_sat.get(bid, (None, None))
        if pen_left is not None or pen_right is not None:
            psi_old = b.psi_old if b.psi_old is not None else psi
            dt = self.dt
            gamma_eff = SAT_DAMPING * b.sat_scale
            if pen_left is not None:
                vel_l = (psi[0] - psi_old[0]) / (dt + 1e-30)
                acc[0] += pen_left - gamma_eff * vel_l
            if pen_right is not None:
                vel_r = (psi[-1] - psi_old[-1]) / (dt + 1e-30)
                acc[-1] += pen_right - gamma_eff * vel_r

        return acc

    def _junction_role_coeffs(self, b) -> Tuple[float, float]:
        k_role = 1.0
        source_role = 0.0
        if b.role == "junction_bridge":
            k_role *= 1.25
        elif b.role == "ring":
            k_role *= 0.85
        elif b.role == "mobius_arm":
            k_role *= 1.10
        return k_role, source_role

    def _mobius_kirchhoff_term(self, b, u_x, g_x, ds) -> float:
        if not b.mobius or b.kappa_mobius <= 0.0:
            return 0.0
        return float(b.orientation) * b.kappa_mobius * (u_x + g_x)

    def _solve_junction_monolithic(self, junc, active, new) -> bool:
        """Monolithic implicit solve for one junction."""
        m = len(active)
        loop_slots = [
            (i, a)
            for i, a in enumerate(active)
            if a["b"].attached_ring and a["b"].kappa_loop > 0.0
        ]
        n_loop = len(loop_slots)

        def unpack(z):
            u = z[0]
            g = z[1 : 1 + m]
            loops = {}
            idx = 1 + m
            for i, a in loop_slots:
                loops[i] = (z[idx], z[idx + 1])
                idx += 2
            return u, g, loops

        def F(z):
            u, g, loops = unpack(z)
            res = np.zeros(1 + m + 2 * n_loop)
            for i, a in enumerate(active):
                res[i] = a["k"] * (g[i] - u) + a["mreg"] * (g[i] - a["inner_val"])
            kir = 0.0
            row = 1 + m
            for i, a in enumerate(active):
                b = a["b"]
                kir += a["flux_sign"] * a["k"] * (u - g[i]) + a["loop_src"]
                if a["mobius_k"] > 0.0:
                    kir += float(b.orientation) * a["mobius_k"] * (u + g[i])
                if i in loops:
                    I1, P1 = loops[i]
                    I0, P0 = b.I_loop, b.Phi_loop
                    L_ind = b.L_loop_inductance
                    I_half = 0.5 * (I0 + I1)
                    P_half = 0.5 * (P0 + P1)
                    drive = b.loop_couple * (u - g[i])
                    gamma_LC = 0.45
                    kir += b.loop_couple * I_half / (L_ind + 1e-30)
                    res[row] = I1 - I0 - self.dt / (L_ind + 1e-30) * (
                        drive - b.kappa_loop * P_half - gamma_LC * I_half
                    )
                    res[row + 1] = P1 - P0 - self.dt * I_half
                    row += 2
            kir -= junc.J_nonlinear * u ** 3
            res[m] = kir
            return res

        u0 = float(np.mean([a["face_val"] for a in active]))
        z0 = [u0] + [a["inner_val"] for a in active]
        for i, a in loop_slots:
            z0 += [a["b"].I_loop, a["b"].Phi_loop]
        z0 = np.array(z0, dtype=float)

        try:
            sol = scipy_root(F, z0, method="hybr")
        except Exception:
            return False
        if (not sol.success) or float(np.max(np.abs(F(sol.x)))) > 1.0:
            return False

        u, g, loops = unpack(sol.x)
        u = float(np.clip(u, -10.0, 10.0))
        g = np.clip(g, -30.0, 30.0)

        for i, a in enumerate(active):
            p = new[a["bid"]][0]
            ds_a = a["ds"]
            psi_face = float(p[a["face_idx"]])
            _sat_scale = a["b"].sat_scale
            pen = (_sat_scale * SAT_PENALTY / self.dt ** 2) * (u - psi_face)
            prev_sat = self._junction_sat.get(a["bid"], (None, None))
            if a["face_idx"] == 0:
                self._junction_sat[a["bid"]] = (pen, prev_sat[1])
            else:
                self._junction_sat[a["bid"]] = (prev_sat[0], pen)
            # Flux diagnostic
            g_val = float(g[i])
            inner2_cur = float(p[a["inner2_idx"]]) if self.N > 3 else g_val
            if a["face_idx"] == 0:
                J = (-3 * psi_face + 4 * g_val - inner2_cur) / (2 * ds_a + 1e-30)
            else:
                J = (3 * psi_face - 4 * g_val + inner2_cur) / (2 * ds_a + 1e-30)
            prev_flux = self._junction_flux.get(a["bid"], (None, None))
            if a["face_idx"] == 0:
                self._junction_flux[a["bid"]] = (J, prev_flux[1])
            else:
                self._junction_flux[a["bid"]] = (prev_flux[0], J)

        for i, (I1, P1) in loops.items():
            b = active[i]["b"]
            b.I_loop = float(I1)
            b.Phi_loop = float(P1)

        self._junction_u[junc.junction_id] = u
        self._junction_state[junc.junction_id] = {
            a["bid"]: {
                "u": float(u),
                "g": float(g[i]),
                "a_inner": float(a["inner_val"]),
                "k": float(a["k"]),
                "mreg": float(a["mreg"]),
                "mobius_k": float(a["mobius_k"]),
            }
            for i, a in enumerate(active)
        }
        return True

    def _apply_junctions(self, new: Dict) -> None:
        """Full nonlinear coupled junction solve + loop-current evolution."""
        net = self.net
        for junc in net.junctions.values():
            active = []
            for bid, sgn in zip(junc.branch_ids, junc.side_signs):
                b = net.branches.get(bid)
                if b is None or not b.active:
                    continue
                p = new[bid][0]
                ds = b.length / (self.N - 1)
                v2 = b.v_speed ** 2
                if sgn == +1:
                    face_idx = self.N - 1
                    inner_idx = self.N - 2
                    inner2_idx = max(self.N - 3, 0)
                else:
                    face_idx = 0
                    inner_idx = 1
                    inner2_idx = min(2, self.N - 1)
                face_val = float(p[face_idx])
                inner_val = float(p[inner_idx])

                k_role, _ = self._junction_role_coeffs(b)
                if b.attached_ring:
                    k_base = (2 * pi * b.v_speed / b.length) / (ds + 1e-30)
                else:
                    k_base = v2 / (ds ** 2)
                k = k_role * k_base
                mreg = 0.50 * k

                loop_src = 0.0
                mobius_k = 0.0
                if b.mobius and b.kappa_mobius > 0.0:
                    mobius_k = b.kappa_mobius

                active.append(dict(
                    bid=bid, b=b, sgn=sgn, ds=ds,
                    face_idx=face_idx, inner_idx=inner_idx,
                    inner2_idx=inner2_idx,
                    face_val=face_val, inner_val=inner_val,
                    k=k, mreg=mreg, flux_sign=sgn,
                    loop_src=loop_src, mobius_k=mobius_k,
                ))
            if not active:
                continue
            self._solve_junction_monolithic(junc, active, new)

    def step(self) -> None:
        """Single Störmer–Verlet time step."""
        net = self.net
        dt = self.dt
        new = {}
        for bid, b in net.branches.items():
            if not b.active or b.psi is None:
                continue
            psi = b.psi
            psi_old = b.psi_old if b.psi_old is not None else psi.copy()
            acc = self._branch_acc(bid, psi)
            psi_new = np.clip(2 * psi - psi_old + dt ** 2 * acc, -50, 50)
            # Endpoint Dirichlet
            if not (b.periodic or b.attached_ring):
                if self._is_ep(b.node_left):
                    psi_new[0] = 0.0
                if b.mobius or self._is_ep(b.node_right):
                    psi_new[-1] = 0.0
            new[bid] = (psi_new, psi.copy())

        self._apply_junctions(new)

        for bid, (psi_new, psi_prev) in new.items():
            b = net.branches[bid]
            if not b.active:
                continue
            # Bridge field step
            if b.bridge is not None and not b.bridge.broken:
                mid = self.N // 2
                psi2 = float(psi_new[mid] ** 2)
                b.bridge.step(psi2, dt)
                if b.bridge.above_threshold:
                    b.bridge.broken = True
            b.psi_old = psi_prev
            b.psi = psi_new

        # Nucleation check
        for bid in list(net.branches.keys()):
            b = net.branches[bid]
            if not b.active or b.bridge is None or not b.bridge.broken:
                continue
            E_before_components = self.total_energy(return_components=True)
            B_val = (pi - A_COULOMB) * HBAR_C
            m_q = sqrt(B_val * SIGMA_FM)
            penalty = 2.0 * m_q
            E_parent_wave_before = 0.0
            if b.psi is not None:
                ds_p = b.length / (self.N - 1)
                vel_p = (b.psi - b.psi_old) / (dt + 1e-30) if b.psi_old is not None else np.zeros(self.N)
                grad_p = np.gradient(b.psi, ds_p)
                E_parent_wave_before = float(_trapezoid(
                    0.5 * vel_p ** 2 + 0.5 * b.v_speed ** 2 * grad_p ** 2, dx=ds_p
                ))
            E_parent_bridge_before = b.bridge.energy() if b.bridge else 0.0
            if E_parent_wave_before < penalty:
                b.bridge.broken = False
                continue
            # Drain energy locally
            if b.psi is not None:
                mid = self.N // 2
                sig = self.N // 6
                mask = np.exp(-0.5 * ((np.arange(self.N) - mid) / max(sig, 1)) ** 2)
                scale = 1.0 - (penalty / (E_parent_wave_before + 1e-30)) * mask
                b.psi *= np.clip(scale, 0.1, 1.0)
                if b.psi_old is not None:
                    b.psi_old *= np.clip(scale, 0.1, 1.0)

            ok, msg = net.nucleate_pair_on_branch(bid)
            if not ok:
                b.bridge.broken = False
                continue
            self._pair_creation_penalty += penalty
            daughter_bids = [d for d in sorted(net.branches.keys()) if d not in self._buf]

            self.nucleation_events.append({
                "bid": bid, "time": self.time,
                "E_before": E_before_components["total"],
                "E_parent_wave_before": E_parent_wave_before,
                "E_parent_bridge_before": E_parent_bridge_before,
                "penalty": penalty, "m_q": m_q,
            })
            for d_bid in daughter_bids:
                db = net.branches[d_bid]
                if db.psi is not None and len(db.psi) > 4:
                    N_d = len(db.psi)
                    s_d = np.linspace(0, db.length, N_d)
                    mode1 = np.sin(pi * s_d / db.length)
                    mode1 /= np.linalg.norm(mode1) + 1e-30
                    proj = float(np.dot(db.psi, mode1))
                    db.psi = proj * mode1
                    db.psi_old = proj * mode1
                n_delay = max(2, int(round(db.tau_delay / self.dt)))
                self._buf[d_bid] = np.zeros((n_delay, self.N))
                self._ptr[d_bid] = 0
                if db.psi is not None:
                    for k_buf in range(n_delay):
                        self._buf[d_bid][k_buf] = db.psi[::-1]

        self.time += self.dt

    def run(self, n_steps: int, record_every: int = 10) -> Dict:
        for k in range(n_steps):
            self.step()
            if k % record_every == 0:
                self.t_hist.append(self.time)
                self.E_hist.append(self.total_energy())
        return {"time": np.array(self.t_hist), "energy": np.array(self.E_hist)}

    def total_energy(self, return_components: bool = False):
        """Component-reporting energy diagnostic."""
        E_wave = E_bridge = E_junction = E_loop = 0.0
        dt = self.dt
        for bid, b in self.net.branches.items():
            if not b.active:
                continue
            N = self.N
            ds = b.length / (N - 1)
            v2 = b.v_speed ** 2
            s = np.linspace(0, b.length, N)
            U = np.array([b.U_potential(si, b.length) for si in s])
            mu2 = self.net.confinement_mass_sq(bid)
            vel = (b.psi - b.psi_old) / (dt + 1e-30)
            kin = 0.5 * vel ** 2
            grad = np.gradient(b.psi, ds)
            E_wave += float(_trapezoid(
                kin + 0.5 * v2 * grad ** 2 + 0.5 * U * b.psi ** 2
                + 0.25 * b.lambda_nl * b.psi ** 4 + 0.5 * mu2 * b.psi ** 2,
                dx=ds,
            ))
            if b.bridge is not None:
                E_bridge += b.bridge.energy()
            if b.attached_ring:
                L_ind = b.L_loop_inductance
                if L_ind > 0:
                    E_loop += 0.5 * L_ind * b.I_loop ** 2
                if b.kappa_loop > 0.0:
                    E_loop += 0.5 * b.kappa_loop * b.Phi_loop ** 2

        for junc in self.net.junctions.values():
            jid = junc.junction_id
            if jid in self._junction_state:
                for bid, st in self._junction_state[jid].items():
                    e_loc = (
                        0.5 * abs(st["k"]) * (st["u"] - st["g"]) ** 2
                        + 0.5 * st["mreg"] * (st["g"] - st["a_inner"]) ** 2
                    )
                    if st.get("mobius_k", 0.0) > 0.0:
                        e_loc += 0.5 * st["mobius_k"] * (st["u"] + st["g"]) ** 2
                    E_junction += e_loc
            else:
                vals = []
                for bid, sgn in zip(junc.branch_ids, junc.side_signs):
                    b = self.net.branches.get(bid)
                    if b and b.active and b.psi is not None:
                        vals.append(float(b.psi[-1] if sgn == +1 else b.psi[0]))
                u_fb = float(np.mean(vals)) if vals else 0.0
                for bid, sgn in zip(junc.branch_ids, junc.side_signs):
                    b = self.net.branches.get(bid)
                    if b and b.active and b.psi is not None:
                        face = float(b.psi[-1] if sgn == +1 else b.psi[0])
                        E_junction += 0.5 * KAPPA_J * (face - u_fb) ** 2

        E_pair = self._pair_creation_penalty
        E_total = E_wave + E_bridge + E_junction + E_loop + E_pair
        if return_components:
            return {
                "wave": E_wave, "bridge": E_bridge,
                "junction": E_junction, "loop": E_loop,
                "pair": E_pair, "total": E_total,
            }
        return E_total

"""
SU(2) spinor transport — spin-½ from Hopf holonomy.

Explicit demonstration that the state acquires a sign flip (e^{iπ} = −1)
after 2π rotation and returns to itself only after 4π.  This is the
geometric origin of spin-½: the holonomy of the Hopf connection at the
poles of S³.
"""

import numpy as np
from geometrodynamics.hopf.connection import hopf_holonomy


def compute_spinor_monodromy(n_pts: int = 401) -> dict:
    """SU(2) spinor transport around 2π and 4π rotations.

    Returns
    -------
    dict with keys:
        angle        : ndarray of rotation angles [0, 4π]
        psi_path     : (n_pts, 2) complex spinor trajectory
        overlap_2pi  : ⟨ψ₀|U(2π)|ψ₀⟩  (should be −1)
        overlap_4pi  : ⟨ψ₀|U(4π)|ψ₀⟩  (should be +1)
        signflip_err : |overlap_2pi + 1|
        return_err   : |overlap_4pi − 1|
        u1_phase_2pi : e^{i·holonomy(χ=0)} = e^{iπ} = −1
        u1_phase_eq  : e^{i·holonomy(χ=π/2)} = 1
    """
    ang = np.linspace(0.0, 4.0 * np.pi, n_pts)
    psi0 = np.array([1.0 + 0.0j, 0.0 + 0.0j])
    psi_path = np.zeros((n_pts, 2), dtype=complex)

    for i, a in enumerate(ang):
        U = np.array(
            [[np.exp(0.5j * a), 0.0j], [0.0j, np.exp(-0.5j * a)]],
            dtype=complex,
        )
        psi_path[i] = U @ psi0

    i2 = int(np.argmin(np.abs(ang - 2.0 * np.pi)))
    i4 = int(np.argmin(np.abs(ang - 4.0 * np.pi)))

    overlap_2pi = np.vdot(psi0, psi_path[i2])
    overlap_4pi = np.vdot(psi0, psi_path[i4])

    u1_phase_2pi = complex(np.exp(1j * hopf_holonomy(0.0)))
    u1_phase_eq = complex(np.exp(1j * hopf_holonomy(np.pi / 2.0)))

    return dict(
        angle=ang,
        psi_path=psi_path,
        overlap_2pi=overlap_2pi,
        overlap_4pi=overlap_4pi,
        signflip_err=abs(overlap_2pi + 1.0),
        return_err=abs(overlap_4pi - 1.0),
        u1_phase_2pi=u1_phase_2pi,
        u1_phase_eq=u1_phase_eq,
    )

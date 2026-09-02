"""
Throat transport: the map ``sigma`` and what it actually is.

CORRECTION (finite-mouth topology round, ``docs/finite_mouth_bundle_transport.md``)
──────────────────────────────────────────────────────────────────────────────
An earlier version of this docstring said that

    sigma(z1, z2) = (conj(z2), -conj(z1))

is "the UNIQUE orientation-reversing spinor map on S^3 that preserves the
Hopf bundle structure" and "a theorem of the Hopf fibration". That is false.
Its real 4x4 matrix has determinant +1: it is left multiplication by the unit
quaternion -j, an isoclinic rotation of S^3 by pi/2, orientation-PRESERVING
on S^3. What it does do, exactly:

    h(sigma z) = -h(z)                antipode on the Hopf base S^2
    sigma(e^{i phi} z) = e^{-i phi} sigma(z)     reversal of the Hopf fibre S^1

two orientation reversals whose product preserves the orientation of the
total space. No orientation-reversing isometry of S^3 preserves the Hopf
fibration at all (fibre-preserving isometries are U(2) and K.U(2), both of
real determinant +1).

In the spin-1/2 representation the same map is

    T = i sigma_y K      for the Hopf complex structure (fibre = scalar phase)
    T = i sigma_y        for the half-spinor complex structure of Spin(4)

i.e. the antilinear K is a statement about which complex structure the state
space carries, not about the geometry. ``T^2 = -I``, ``T^dagger T = I``,
``det T = 1`` hold; they are properties of the quaternionic structure of C^2.
The finite-mouth round (``geometrodynamics/bulk/mouth_topology.py``,
``docs/finite_mouth_bundle_transport.md`` section 5b) found that the deck-
generator holonomy of the RP^2 mouth on the restricted Pin+ spinor bundle
is left multiplication by a unit quaternion of the (i, j) plane -- the same
fibre-reversing component of Pin(2) that sigma = L_{-j} belongs to, with
the same square -1 -- determined up to a Spin(2) = U(1) direction and a
sign. So T^2 = -1 is geometric (the spin holonomy of a 2 pi rotation on the
round neck). The direction -j within that component is gauge (a Spin(2)
conjugation), and the Hopf fibre IS the mouth's Spin(2)
(``bulk/mouth_spin_frame.py``: the bulk mouth S^3 is the spin-frame bundle
of the brane mouth S^2, fibre angle phi = frame angle 2 phi), so the
antilinear K is canonical. The sign is the Pin^- sector of the mouth, and
Pin^- bordism requires the two mouths of a created pair to carry opposite
sectors. Still chosen: the antipodal quotient construction itself.

Verification checks retained:
    T^2 = -I, T^dagger T = I, det T = 1, T sigma_z T^dagger = -sigma_z.
"""

from __future__ import annotations

import numpy as np


# Pauli matrices
_SX = np.array([[0, 1], [1, 0]], dtype=complex)
_SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
_SZ = np.array([[1, 0], [0, -1]], dtype=complex)
_ID = np.eye(2, dtype=complex)


def derive_throat_transport() -> np.ndarray:
    """Return T = iσ_y, the matrix of ``sigma(z1, z2) = (conj z2, -conj z1)``
    on real coefficients.

    ``sigma`` is left multiplication by ``-j`` on ``S^3`` (orientation-
    preserving, ``det = +1``); it is NOT derived from an orientation reversal.
    The name is kept for the callers that depend on it.

    Returns the 2×2 transport matrix.
    """
    return 1j * _SY


def verify_transport_properties(T: np.ndarray) -> dict:
    """Verify all required properties of the throat transport.

    Returns a dict of property names → (value, passed) pairs.
    """
    results = {}

    # 1. T² = −I  (double cover / 4π periodicity)
    T2 = T @ T
    t2_check = np.allclose(T2, -_ID)
    results["T²=−I"] = (float(np.max(np.abs(T2 + _ID))), t2_check)

    # 2. Unitarity: T†T = I
    TdT = T.conj().T @ T
    unit_check = np.allclose(TdT, _ID)
    results["T†T=I"] = (float(np.max(np.abs(TdT - _ID))), unit_check)

    # 3. det(T) = 1  (SU(2), not U(2))
    det_T = complex(np.linalg.det(T))
    det_check = abs(det_T - 1.0) < 1e-12
    results["det=1"] = (abs(det_T - 1.0), det_check)

    # 4. Orientation reversal: T σ_z T† = −σ_z
    flipped = T @ _SZ @ T.conj().T
    flip_check = np.allclose(flipped, -_SZ)
    results["TσzT†=−σz"] = (float(np.max(np.abs(flipped + _SZ))), flip_check)

    # 5. Preserves σ_x and σ_y under conjugation (same as iσ_y action)
    sx_conj = T @ _SX @ T.conj().T
    sx_check = np.allclose(sx_conj, -_SX)
    results["TσxT†=−σx"] = (float(np.max(np.abs(sx_conj + _SX))), sx_check)

    # 6. Action on basis states
    up = np.array([1, 0], dtype=complex)
    dn = np.array([0, 1], dtype=complex)
    t_up = T @ up
    t_dn = T @ dn
    up_check = np.allclose(t_up, -dn)  # T|↑⟩ = −|↓⟩
    dn_check = np.allclose(t_dn, up)   # T|↓⟩ = +|↑⟩
    results["T|↑⟩=−|↓⟩"] = (float(np.max(np.abs(t_up + dn))), up_check)
    results["T|↓⟩=+|↑⟩"] = (float(np.max(np.abs(t_dn - up))), dn_check)

    return results


def orientation_reversal_on_s3(z1: complex, z2: complex) -> tuple[complex, complex]:
    """``sigma(z1, z2) = (conj z2, -conj z1)``: left multiplication by ``-j``.

    Despite the name (kept for callers), this map PRESERVES the orientation of
    ``S^3`` (real determinant +1). It reverses the orientation of the Hopf
    base and of the fibre. See the module docstring.
    """
    return (z2.conjugate(), -z1.conjugate())


def verify_hopf_preservation(n_samples: int = 1000, rng=None) -> float:
    """Verify σ preserves the Hopf fibration: π(σ(p)) is well-defined.

    The Hopf map π: S³ → S² sends (z₁,z₂) to the Bloch sphere point
    (2Re(z̄₁z₂), 2Im(z̄₁z₂), |z₁|²−|z₂|²).

    Under σ the base point on S² goes to its antipode and the fibre
    phase is reversed.  This check verifies that σ maps fibers to fibers
    (the base point of σ(z) is the antipode of the base point of z, for
    every point on the fibre) and that |σ(z)| = 1.

    Returns the maximum fiber-mapping error across random samples.

    If ``rng`` is None, a fixed-seed generator (seed 42) is used so the
    check is deterministic; pass an unseeded Generator to sample fresh
    points.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    max_err = 0.0
    for _ in range(n_samples):
        # Random point on S³
        v = rng.standard_normal(4)
        v = v / np.linalg.norm(v)
        z1 = complex(v[0], v[1])
        z2 = complex(v[2], v[3])

        # Apply orientation reversal
        w1, w2 = orientation_reversal_on_s3(z1, z2)

        # Check |w|² = 1
        norm_err = abs(abs(w1) ** 2 + abs(w2) ** 2 - 1.0)
        max_err = max(max_err, norm_err)

        # Fibre-to-fibre: the Hopf base point of sigma(e^{i phi} z) must be
        # the antipode of the base point of z, for every phase phi.
        def base(a: complex, b: complex) -> np.ndarray:
            return np.array([2.0 * (a.conjugate() * b).real,
                             2.0 * (a.conjugate() * b).imag,
                             abs(a) ** 2 - abs(b) ** 2])

        phase = np.exp(1j * rng.uniform(0.0, 2.0 * np.pi))
        u1, u2 = orientation_reversal_on_s3(phase * z1, phase * z2)
        fibre_err = float(np.max(np.abs(base(u1, u2) + base(z1, z2))))
        max_err = max(max_err, fibre_err)

    return max_err


def derive_singlet_from_transport() -> np.ndarray:
    """Construct the singlet state from the derived transport.

    |Ψ⟩ = N × Σ_s |s⟩ ⊗ T|s⟩

    where T is the Hopf-derived throat transport and the sum is
    over the spin basis.  Normalised to unit length.

    This is the SAME construction used in bell/analyzers.py.  The chain

        Hopf fibration → orientation reversal → T = iσ_y → singlet

    that this function's docstring once advertised is NOT a derivation:
    the first link is false (σ preserves the orientation of S³), and the
    singlet is this formula's definition.  See the module docstring.
    """
    T = derive_throat_transport()
    up = np.array([1, 0], dtype=complex)
    dn = np.array([0, 1], dtype=complex)
    psi = np.kron(up, T @ up) + np.kron(dn, T @ dn)
    return psi / np.linalg.norm(psi)

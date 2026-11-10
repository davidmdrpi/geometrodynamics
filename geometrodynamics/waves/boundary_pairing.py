"""Baseline qualifications for the prospective conjugate-mouth experiment.

The endpoint-indexed Maxwell DtN map, spin-frame epsilon pairing and detector
Gram form are distinct objects. Nothing here evolves a responding mouth or
derives probabilities. See docs/conjugate_boundary_pairing_protocol.md.
"""
from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

from geometrodynamics.embedding.transport import derive_throat_transport
from .maxwell_throat_response import Bridge, field_basis, solve_port


def _matrix(value, shape=None):
    value = np.asarray(value, dtype=complex)
    if value.ndim != 2 or (shape is not None and value.shape != shape):
        raise ValueError("Wrong matrix shape")
    if not np.isfinite(value).all():
        raise ValueError("Finite matrix required")
    return value


def transfer_to_dtn(transfer, pole_tolerance=1e-12):
    """Map (u,j)_L to (u,j)_R into outward flux (-j_L,j_R)=D(u_L,u_R).

    Retain the numerically computed determinant: imposing det=1 here would
    hide a reciprocity error. A coordinate pole is not a scattering failure.
    """
    transfer = _matrix(transfer, (2, 2))
    if not np.isfinite(pole_tolerance) or pole_tolerance <= 0:
        raise ValueError("Positive finite pole tolerance required")
    a, b, c, d = transfer.ravel()
    if abs(b) <= pole_tolerance * max(1., np.linalg.norm(transfer)):
        raise ValueError("Unresolved Dirichlet coordinate pole")
    return np.array([[a / b, -1 / b], [c - d * a / b, d / b]])


def interval_dtn(bridge: Bridge, omega, left, right, ell=1, rtol=2e-11):
    """Integrate a segment of the inherited bridge without mirroring its ends."""
    if not np.isfinite([omega, left, right, rtol]).all() or omega <= 0 or rtol <= 0:
        raise ValueError("Finite positive frequency/tolerance required")
    if not -bridge.limit <= left < right <= bridge.limit:
        raise ValueError("Ordered interval must lie inside the bridge")
    if not np.isfinite(ell) or int(ell) != ell or ell < 1:
        raise ValueError("Radiative Maxwell ell must be an integer >=1")

    def rhs(x, flat):
        p = bridge.profile(x)
        n, r = p["lapse"], p["radius"]
        generator = bridge.stretch * p["jacobian"] * np.array(
            [[0., 1 / n], [n * ell * (ell + 1) / r**2 - omega**2 / n, 0.]])
        return (generator @ flat.reshape(2, 2)).ravel()

    sol = solve_ivp(rhs, (left, right), np.eye(2).ravel(), method="DOP853",
                    rtol=rtol, atol=rtol * .01)
    if not sol.success:
        raise RuntimeError(sol.message)
    return transfer_to_dtn(sol.y[:, -1].reshape(2, 2))


def join_dtn(segments, condition_limit=1e12):
    """Join adjacent outward-flux maps by stationary interface elimination.

    No conjugation or sign is assigned to the common field value. Outward
    fluxes sum to zero. This is static coherent joining, not tensor composition.
    """
    segments = [_matrix(d, (2, 2)) for d in segments]
    if not segments or not np.isfinite(condition_limit) or condition_limit <= 1:
        raise ValueError("Segments and a finite condition limit >1 required")
    assembled = np.zeros((len(segments) + 1, len(segments) + 1), complex)
    for i, d in enumerate(segments):
        assembled[i:i + 2, i:i + 2] += d
    boundary = [0, len(segments)]
    if len(segments) == 1:
        return assembled
    interior = list(range(1, len(segments)))
    hbb = assembled[np.ix_(boundary, boundary)]
    hbi = assembled[np.ix_(boundary, interior)]
    hii = assembled[np.ix_(interior, interior)]
    hib = assembled[np.ix_(interior, boundary)]
    # A scalar block can have condition number 1 while being nearly zero.
    singular = np.linalg.svd(hii, compute_uv=False)
    scale = max(1., np.linalg.norm(assembled, 2))
    if singular[-1] <= scale / condition_limit:
        raise ValueError("Unresolved singular or ill-conditioned interface")
    return hbb - hbi @ np.linalg.solve(hii, hib)


def scattering_from_dtn(dtn, omega):
    dtn = _matrix(dtn, (2, 2))
    if not np.isfinite(omega) or omega <= 0:
        raise ValueError("Positive finite frequency required")
    z = 1j * omega * np.eye(2)
    return -np.linalg.solve(dtn - z, dtn + z)


def direct_waist_load(bridge, omega, incoming, ell=1, order=192):
    """Integrate the coherent internal field square, without a load matrix.

    Incoming amplitudes use the inherited incident-power normalization.
    This is mechanical load, not a count rate or a probability.
    """
    incoming = np.asarray(incoming, complex)
    if incoming.shape != (2,) or not np.isfinite(incoming).all():
        raise ValueError("Two finite incident amplitudes required")
    if not isinstance(order, (int, np.integer)) or order < 2:
        raise ValueError("Quadrature order must be an integer >=2")
    port = solve_port(bridge, omega, ell, rtol=2e-12)
    x, w = leggauss(order)
    x, w = x * bridge.limit, w * bridge.limit
    p = bridge.profile(x)
    u, _ = field_basis(port, x)
    weight = (w * bridge.stretch * p["jacobian"] * p["lapse"]
              * ell * (ell + 1) * p["shape"] / (omega**2 * p["radius"]**2))
    return float(np.sum(weight * abs(u @ incoming)**2))


def polarization_matrix(quadratic):
    """Recover a Hermitian two-input response from six real experiments."""
    a, b = np.eye(2, dtype=complex)
    real = (quadratic(a + b) - quadratic(a - b)) / 4
    imag = (quadratic(a - 1j * b) - quadratic(a + 1j * b)) / 4
    return np.array([[quadratic(a), real + 1j * imag],
                     [real - 1j * imag, quadratic(b)]])


def output_gram(scattering, inputs, port=1):
    """Positive outgoing-power Gram for ONE fixed-background experiment."""
    scattering = _matrix(scattering, (2, 2))
    inputs = _matrix(inputs)
    if inputs.shape[1] != 2 or port not in (0, 1):
        raise ValueError("Two input columns and port 0 or 1 required")
    outputs = inputs @ scattering[port, :]
    return np.outer(outputs.conj(), outputs)


def third_order_response(gram):
    """Inclusion-exclusion for three alternatives in one fixed Gram form."""
    gram = _matrix(gram, (3, 3))
    def response(indices):
        return np.sum(gram[np.ix_(indices, indices)]).real
    return float(response([0, 1, 2]) - response([0, 1]) - response([0, 2])
                 - response([1, 2]) + response([0]) + response([1]) + response([2]))


def context_gram_counterexample():
    """Four individually PSD outcome Gram matrices with CHSH=4.

    These are constructed controls, not BAM predictions or a joint history
    functional. Per-context PSD plus unbiased marginals does not bound CHSH.
    """
    rows = []
    correlations = []
    for x, y in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        target = (-1)**(x * y)
        p = np.array([(1 + s * t * target) / 4
                      for s in [1, -1] for t in [1, -1]]).reshape(2, 2)
        gram = np.diag(p.ravel())
        correlations.append(float(np.sum(p * np.array([[1, -1], [-1, 1]]))))
        rows.append(dict(settings=[x, y], probabilities=p.tolist(),
                         min_eigenvalue=float(np.linalg.eigvalsh(gram).min()),
                         normalization_error=float(abs(p.sum() - 1)),
                         marginal_error=float(np.max(abs(p.sum(axis=0) - .5)))))
    return dict(contexts=rows, chsh=abs(sum(correlations[:3]) - correlations[3]),
                common_joint_functional="NOT_SUPPLIED")


def conditional_spinor_checks(samples=256, seed=7281):
    """Known epsilon identity and controls; no amplitude is derived here."""
    if not isinstance(samples, (int, np.integer)) or samples < 1:
        raise ValueError("Positive integer sample count required")
    rng = np.random.default_rng(seed)
    epsilon = derive_throat_transport()
    sigma = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
    worst = dict(identity=0., normalization=0., marginal=0., rotation=0., phase=0.,
                 global_pin_sign=0., transport_square=0.)
    for _ in range(samples):
        z, w, u = rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2))
        z, w, u = [v / np.linalg.norm(v) for v in [z, w, u]]
        rotate = np.array([[u[0], u[1]], [-u[1].conjugate(), u[0].conjugate()]])
        a = np.einsum('i,kij,j->k', z.conj(), sigma, z).real
        b = np.einsum('i,kij,j->k', w.conj(), sigma, w).real
        za, wb = [z, epsilon @ z.conj()], [w, epsilon @ w.conj()]
        p = np.array([[abs(v @ epsilon @ t)**2 / 2 for t in wb] for v in za])
        target = (1 - np.outer([1, -1], [1, -1]) * (a @ b)) / 4
        worst['identity'] = max(worst['identity'], float(np.max(abs(p - target))))
        worst['normalization'] = max(worst['normalization'], float(abs(p.sum() - 1)))
        worst['marginal'] = max(worst['marginal'], float(np.max(abs(p.sum(axis=0) - .5))))
        value = abs(z @ epsilon @ w)**2
        worst['rotation'] = max(worst['rotation'], float(abs(abs((rotate @ z) @ epsilon @ (rotate @ w))**2 - value)))
        phases = np.exp(1j * rng.uniform(-np.pi, np.pi, 2))
        worst['phase'] = max(worst['phase'], float(abs(abs((phases[0]*z) @ epsilon @ (phases[1]*w))**2 - value)))
        worst['global_pin_sign'] = max(worst['global_pin_sign'], float(abs(abs(z @ (-epsilon) @ w)**2 - value)))
        worst['transport_square'] = max(worst['transport_square'], float(np.max(abs(epsilon @ (epsilon @ z.conj()).conj() + z))))
    # Chosen unequal efficiencies preserve unbiased marginals but change E.
    c = np.cos(1.)
    base = (1 - np.outer([1, -1], [1, -1]) * c) / 4
    weighted = base * np.array([[2., 1.], [1., 2.]])
    weighted /= weighted.sum()
    correlation = float(np.sum(weighted * np.array([[1, -1], [-1, 1]])))
    return dict(residuals=worst, unequal_efficiency_correlation=correlation,
                unequal_efficiency_gap=abs(correlation + c),
                scope="CONDITIONAL_SPINOR_IDENTITY")

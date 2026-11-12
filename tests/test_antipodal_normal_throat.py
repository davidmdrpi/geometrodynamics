"""Hamiltonian, quotient and response controls for freeze 533f925."""

import numpy as np
import pytest

from geometrodynamics.bulk import antipodal_normal_throat as ant
from experiments.closure_ledger import antipodal_normal_throat_probe as probe


@pytest.fixture
def model():
    return ant.Throat(8, 1)


def test_vector_field_comes_from_full_hamiltonian(model):
    z = np.random.default_rng(12).normal(size=model.dim)*.13
    eps = 1e-6
    e = np.eye(model.dim)
    gradient = (model.energy(z+eps*e)-model.energy(z-eps*e))/(2*eps)
    assert np.max(abs(model.symplectic()@model.rhs(z)[0]+gradient)) < 1e-7
    # The normal force must include the field variation, not only support.
    assert np.max(abs(model.loads(z))) > 1e-3


@pytest.mark.parametrize("kind", ant.CELLS)
def test_coordinate_identification_and_real_field_pullback(model, kind):
    rng = np.random.default_rng(19)
    points = rng.uniform(.1, 6., (10, 3))
    assert np.allclose(ant.identify(ant.identify(points, kind), kind), points)
    x = rng.normal(size=(1, 2, 2, model.sites))
    z = model.pack(x, np.zeros_like(x), np.zeros((1, 2)), np.zeros((1, 2)))
    xx, *_ = model.unpack(model.deck(z, kind))
    theta = .31
    for face, u in enumerate((1., -1.)):
        for i, phi in enumerate(model.phi):
            target = ant.identify([phi, theta, u], kind)
            j = int(round(target[0]*model.sites/(2*np.pi))) % model.sites
            other = 0 if target[2] > 0 else 1
            pulled = x[0, 0, other, j]*np.cos(model.n*target[1])+x[0, 1, other, j]*np.sin(model.n*target[1])
            represented = xx[0, 0, face, i]*np.cos(model.n*theta)+xx[0, 1, face, i]*np.sin(model.n*theta)
            assert pulled == pytest.approx(represented, abs=1e-12)


def test_fixed_space_remains_symplectic_when_time_is_kept(model):
    d = model.deck(np.eye(model.dim)).T
    values, vectors = np.linalg.eigh(d)
    fixed = vectors[:, values > .5]
    restricted = fixed.T@model.symplectic()@fixed
    assert np.linalg.matrix_rank(restricted) == fixed.shape[1]


@pytest.mark.parametrize("kind", ant.CELLS)
@pytest.mark.parametrize("eta", [-1, 1])
def test_decks_are_canonical_symmetries_at_nonzero_mouth_displacement(model, kind, eta):
    z = np.random.default_rng(3).normal(size=(3, model.dim))*.2
    d = model.deck(np.eye(model.dim), kind, eta).T
    assert np.max(abs(d@d-np.eye(model.dim))) < 1e-12
    assert np.max(abs(d.T@model.symplectic()@d-model.symplectic())) < 1e-12
    assert np.allclose(model.energy(z@d.T), model.energy(z), atol=1e-12)
    assert np.max(abs(model.rhs(z@d.T)-model.rhs(z)@d.T)) < 1e-12


@pytest.mark.parametrize("n", [1, 2])
@pytest.mark.parametrize("eta", [-1, 1])
def test_selection_has_a_negative_control(n, eta):
    model = ant.Throat(8, n)
    for m in range(5):
        epsilon = ant.normal_parity(m, n, eta)
        for chosen, expected in [(epsilon, 1), (-epsilon, -1)]:
            x = np.zeros((1, 2, 2, 8)); x[0, 0] = model.mode(m, chosen)
            z = model.pack(x, np.zeros_like(x), np.zeros((1, 2)), np.zeros((1, 2)))
            assert np.max(abs(model.deck(z, eta=eta)-expected*z)) < 1e-12
            g = x[0, 0].ravel()
            assert np.max(abs(model.frozen_matrix()@g-model.frequency_sq(m, chosen)*g)) < 1e-12


def test_circulation_allowed_by_antipodal_but_not_reflective(model):
    z = model.travelling(.3)
    assert np.max(abs(model.deck(z)-z)) < 1e-12
    assert model.angular_momentum(z)[0] > .1
    assert np.max(abs(model.deck(z, "reflective")-z)) > .1
    rng = np.random.default_rng(7)
    for eta in (-1, 1):
        states = model.project(rng.normal(size=(8, model.dim)), "reflective", eta)
        assert np.max(abs(model.angular_momentum(states))) < 1e-12
        # Zero total does not imply zero local current or empty fields.
        x, p, _, _ = model.unpack(states)
        assert np.max(abs(x[:, 0]*p[:, 1]-x[:, 1]*p[:, 0])) > .01


def test_initial_angular_momentum_is_continuous_not_quantized(model):
    low, high = model.travelling(.1), model.travelling(.2)
    assert model.angular_momentum(high)[0] == pytest.approx(4*model.angular_momentum(low)[0])


def test_normal_flux_is_the_face_noether_derivative(model):
    z = np.random.default_rng(41).normal(size=model.dim)*.2
    dz = model.rhs(z)[0]
    eps = 1e-6
    rate = (model.angular_momentum(z+eps*dz, per_face=True)
            - model.angular_momentum(z-eps*dz, per_face=True))[0]/(2*eps)
    flux = model.normal_angular_flux(z)[0]
    assert abs(flux) > 1e-4
    assert np.allclose(rate, [flux, -flux], atol=1e-10)
    fixed = model.project(z)
    assert abs(model.normal_angular_flux(fixed)[0]) < 1e-12
    assert np.max(abs(model.normal_angular_flux(fixed, local=True))) > 1e-4


def test_transfer_control_distinguishes_circulation_from_exchange():
    r = probe.transfer_report()
    assert r["max_local_normal_flux"] > 1e-5
    assert r["max_net_normal_flux"] < 1e-12
    assert r["max_face_momentum_change"] < 1e-12
    assert r["control_initial_net_flux"] > .01
    assert r["control_face_momentum_change"] > .001


def test_signed_normal_load_is_not_opposite_angular_momentum(model):
    z = model.travelling(.3)
    loads = model.loads(z)[0]
    momenta = model.angular_momentum(z, per_face=True)[0]
    assert loads[0] == pytest.approx(loads[1]) and loads[0] > 0
    assert momenta[0] == pytest.approx(momenta[1]) and momenta[0] > 0
    signed = loads*np.array([1., -1.])
    assert signed.sum() == pytest.approx(0.)


def test_kappa_shifts_only_odd_normal_channel():
    a, b = ant.Throat(8, 1, kappa=0.), ant.Throat(8, 1, kappa=1.2)
    assert b.frequency_sq(1, 1) == pytest.approx(a.frequency_sq(1, 1))
    assert b.frequency_sq(1, -1)-a.frequency_sq(1, -1) == pytest.approx(2.4)


def test_susceptibility_agrees_with_independent_grid_source(model):
    omega = .7
    for m in (0, 1, 2):
        wave = np.cos(m*model.phi); wave /= np.linalg.norm(wave)
        source = np.r_[wave, np.zeros(model.sites)]
        solution = np.linalg.solve(model.frozen_matrix()-omega**2*np.eye(2*model.sites), source)
        measured = solution.reshape(2, -1)@wave
        assert np.allclose(measured, model.susceptibility(m, omega)@np.array([1., 0.]), atol=1e-12)


def test_verlet_is_reversible_and_preserves_angular_momentum(model):
    z = np.random.default_rng(11).normal(size=model.dim)*.1
    forward = ant.step(model, z, .01)
    assert np.max(abs(ant.step(model, forward, -.01)[0]-z)) < 1e-14
    assert model.angular_momentum(forward)[0] == pytest.approx(model.angular_momentum(z)[0], abs=1e-14)


def test_responding_flow_moves_mouths_without_losing_the_seam(model):
    z = model.travelling(.4)
    hist = ant.flow(model, z, duration=4.)[:, 0]
    assert abs(hist[:, -4:-2]).max() > .01
    assert np.max(abs(model.deck(hist)-hist)) < 1e-12
    assert np.max(abs(model.angular_momentum(hist)-model.angular_momentum(z)[0])) < 1e-12
    frozen = ant.flow(model, z, duration=4., respond=False)[:, 0]
    assert np.max(abs(frozen[:, -4:])) == 0
    assert np.max(abs(hist[:, :model.field_size]-frozen[:, :model.field_size])) > 1e-5


def test_structural_probe_uses_all_resolved_modes_and_both_sizes():
    r = probe.structural_report()
    assert r["resolved_modes"] == 160
    assert max(r["errors"].values()) < 1e-10
    assert all(1.8 < c["order"] < 2.2 for c in r["convergence"])


def test_response_scaling_and_independent_integrator():
    runs = [probe.responding_report(16, 1, duration=4.)]
    num = probe.numerical_report(duration=4.)
    checks = probe.evaluate(probe.structural_report(), runs, num)
    assert all(checks.values()), checks
    assert num["improvement"] > 3.5


def test_missing_or_failed_checks_cannot_pass():
    complete = {k: True for k in probe.REQUIRED}
    assert ant.classify(complete, probe.REQUIRED)["numerical_validity"] == "PASS"
    complete["momentum"] = False
    assert ant.classify(complete, probe.REQUIRED)["numerical_validity"] == "UNRESOLVED_NUMERICALLY"
    del complete["geometry"]
    assert ant.classify(complete, probe.REQUIRED)["numerical_validity"] == "INCOMPLETE"


def test_invalid_grids_and_recording_intervals():
    with pytest.raises(ValueError):
        ant.Throat(7)
    with pytest.raises(ValueError):
        ant.Throat(8, 0)
    model = ant.Throat(8)
    with pytest.raises(ValueError):
        ant.flow(model, model.travelling(.1), duration=.031)

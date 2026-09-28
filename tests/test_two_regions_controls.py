"""Review controls must distinguish off-shell comparisons from Einstein data."""
import hashlib
import json
from pathlib import Path

import numpy as np

from geometrodynamics.waves import two_regions as tr
from experiments.closure_ledger import two_regions_controls as controls

ROOT = Path(__file__).resolve().parents[1]


def test_active_density_from_anisotropic_sigma_stress():
    phi = np.array([.3, .2, -.4, .1])
    pi = np.array([.5, -.2, .7, .4])
    gradient = np.array([[.2, .4, -.1, .3], [-.3, .5, .2, .1], [.6, .1, -.2, .4]])
    f = 1-phi @ phi/6
    G = np.eye(4)/f+np.outer(phi, phi)/(6*f*f)
    kinetic = pi @ G @ pi
    spatial_tensor = gradient @ G @ gradient.T
    spatial, potential = np.trace(spatial_tensor), 1.5/f**2
    rho = .5*(kinetic+spatial)+potential
    stress = spatial_tensor+np.eye(3)*(.5*(kinetic-spatial)-potential)
    np.testing.assert_allclose(rho+np.trace(stress), 2*kinetic-2*potential, atol=2e-15)
    np.testing.assert_allclose(controls.active_density(kinetic, spatial, potential), rho+np.trace(stress), atol=2e-15)
    np.testing.assert_allclose(controls.active_density(0., spatial, potential), -2*potential, atol=2e-15)
    assert controls.background_active_density(0.) < 0 < controls.background_active_density(np.pi/4)


def test_review_archive_provenance_preserves_original_failed_gate():
    original = json.loads((ROOT/controls.ARCHIVE).read_text())
    review = json.loads((ROOT/'experiments/closure_ledger/runs/20260927_two_regions/review_controls.json').read_text())
    for record in (original, review):
        for path, digest in record['source_sha256'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
    assert review['original_verdict'] == original['verdict'] == 'REGISTERED_GATE_FAILURE'
    assert original['gates']['all_successive_field_differences'] is False


def test_fixed_metric_control_replays_and_is_not_constraint_complete():
    original = json.loads((ROOT/controls.ARCHIVE).read_text())
    review = json.loads((ROOT/'experiments/closure_ledger/runs/20260927_two_regions/review_controls.json').read_text())
    data = controls.load_control(original, 'pair', 'round')
    assert data.amplitudes == (.02, .03)
    np.testing.assert_array_equal(data.coefficients, original['cases']['round'][-1]['coefficients'])
    d = tr.diagnostics(data)
    expected = review['grids'][0]['controls'][1]
    assert expected['constraint_solved'] is False
    assert d['H_normalized_max'] > .001
    np.testing.assert_allclose([r['direct_rate_01'] for r in d['ledger']], expected['rates_01'], rtol=1e-10, atol=1e-13)
    assert d['ledger'][0]['direct_rate_01'] > 1.5e-4
    response = review['grids'][-1]['energy_response']
    assert response['solved'][0][1] < -.02 < 0 < response['fixed_round'][0][1] < .002
    np.testing.assert_allclose(np.array(response['solved'])-response['fixed_round'], response['solved_minus_fixed_round'], atol=1e-14)


def test_seed_overlap_supremum_is_geometric_not_grid_dependent():
    separation, radius = 1.2, .45
    x = np.array([np.cos(radius), np.sin(radius), 0., 0.])
    boundary_value = np.cosh(8*(x @ tr.CENTERS[1]))/np.cosh(8)
    expected = np.cosh(8*np.cos(separation-radius))/np.cosh(8)
    np.testing.assert_allclose(boundary_value, expected, atol=1e-15)
    assert .11 < expected < .12

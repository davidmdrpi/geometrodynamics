import numpy as np
from geometrodynamics.waves import r3_family as rf
from geometrodynamics.waves import esu_floquet as fl


def test_section_round_trip_at_anisotropic_point():
    rng = np.random.default_rng(17)
    z = np.r_[1.002, .003, .1*rng.normal(size=10)]
    assert np.abs(rf.to_section(rf.to_state(z))-z).max() < 1e-13


def test_exact_jet_linearisation_at_the_esu():
    Z0 = np.r_[1., np.zeros(11)]
    p, J = rf.DP(Z0, 12)
    ev = np.linalg.eigvals(J)
    assert np.abs(p-Z0).max() < 1e-12
    assert abs(abs(ev).max()-85.0196952) < 1e-5 and abs(abs(ev).min()-.0117620) < 1e-6
    tr = np.trace(fl.monodromy('T', 2))
    assert max(abs(J[2+2*k, 2+2*k]+J[3+2*k, 3+2*k]-tr) for k in range(5)) < 1e-8


def test_exact_jet_linearisation_matches_finite_differences_off_family():
    zt = np.r_[1.0004, -.002, .05, .03, -.04, .1]
    p, Jj = rf.DP(zt, 6)
    Jf = rf._fd(lambda u: rf.P(np.r_[u, np.zeros(6)])[:6], zt, 1e-6)
    assert np.abs(Jj-Jf).max()/np.abs(Jj).max() < 1e-7


def test_loop_action_of_an_ellipse():
    t = np.linspace(0, 2*np.pi, 200, endpoint=False)
    Z = np.zeros((200, 6))
    Z[:, 2], Z[:, 3] = .3*np.cos(t), -.7*np.sin(t)      # area pi*.21 -> action .105
    assert abs(rf.loop_action(Z)-.105) < 1e-8


def test_pair_classification():
    out = rf.classify_pairs([85., 1/85., np.exp(.8j), np.exp(-.8j)])
    assert sorted(p['kind'] for p in out) == ['ELLIPTIC', 'HYPERBOLIC']

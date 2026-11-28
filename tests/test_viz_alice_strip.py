"""The four-panel charge-without-charge strip on six Alice handles: geometry of the schematic."""
import numpy as np

from geometrodynamics.viz import alice_strip as st
from geometrodynamics.viz.alice_circulation import RIM, flux_through_tube_allowed


def test_six_handles_twelve_distinct_mouths_on_one_surface():
    mouths = st.mouths_on_single_surface()
    assert len(mouths) == 6
    flat = [m for pair in mouths for m in pair]
    for a, b in mouths:
        assert np.allclose(a, -b)
    for i in range(12):
        for j in range(i + 1, 12):
            assert np.linalg.norm(flat[i] - flat[j]) > 2 * RIM
    assert all(np.linalg.norm(m) + RIM < st.SHEET_RADIUS for m in flat)


def test_three_sourceless_flux_loops():
    q = st.mouth_charges()
    assert len(q) == 12
    for sheet in ("top", "bottom"):
        assert sum(q[(sheet, i)] for i in range(6)) == 0
    for i in range(6):
        assert q[("top", i)] == -q[("bottom", i)]
    for down, up in st.CIRCUITS:
        assert q[("top", down)] == -1 and q[("top", up)] == +1


def test_every_seam_reverses_time_and_charge_and_admits_the_flux():
    for s_t, s_f in st.seam_signs().values():
        assert s_t == -1 and flux_through_tube_allowed(s_t, s_f)
        assert st.winding_after_seam(+1, s_f) == -1


def test_each_loop_is_one_closed_worldline_with_two_waist_flips():
    s = np.linspace(0, 1, 6001)
    for k in range(3):
        pts = [st.worldline(k, x) for x in s]
        P = np.array([p.xyz for p in pts])
        assert np.allclose(P[0], P[-1], atol=1e-9)
        assert np.linalg.norm(np.diff(P, axis=0), axis=1).max() < 0.05
        flips = [pts[i].xyz for i in range(1, len(pts) - 1) if pts[i].orientation != pts[i - 1].orientation]
        assert len(flips) == 2 and all(abs(x[2]) < 0.01 for x in flips)
        assert {p.orientation for p in pts if p.sheet == "top"} == {+1}
        assert {p.orientation for p in pts if p.sheet == "bottom"} == {-1}


def test_pulse_crosses_the_equator_and_refocuses_at_the_source_image():
    c = st.COLUMNS[st.WAVE_SOURCE]
    xy, top = st.wavefront(1e-3)
    assert top.all() and np.abs(xy - c).max() < 0.01
    xy, top = st.wavefront(np.pi / 2 + 0.2)
    assert top.any() and (~top).any()
    xy, top = st.wavefront(np.pi - 1e-3)
    assert (~top).all() and np.abs(xy - st.identify(c)).max() < 0.01   # bottom −c, which is the source


def test_render_smoke(tmp_path):
    st.render(tmp_path / "strip.png", dpi=30)
    assert (tmp_path / "strip.png").stat().st_size > 10_000

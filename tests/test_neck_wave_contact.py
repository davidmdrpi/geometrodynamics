"""Check the smooth field and its diagnostics, not a solved fluid/GR model."""
import json
from pathlib import Path
import shutil
import subprocess
import numpy as np
import pytest
from scripts.neck_wave_contact import ContactNeck, Parameters, H, X


@pytest.mark.parametrize('amplitude', [0., .08])
def test_null_controls_remain_open(amplitude):
    m = ContactNeck(Parameters(amplitude=amplitude))
    assert not m.band_events and not m.field_events
    assert max(m.axis_peak(t) for t in np.linspace(0, 9, 37)) < .5
    if amplitude == 0:
        np.testing.assert_array_equal(m.radius(X, 4.25), m.rest(X))


def test_threshold_includes_grazing_and_reopens():
    m = ContactNeck()
    assert m.band_events == pytest.approx([3.53009254206, 4.96990745794], abs=1e-9)
    assert m.field_events == pytest.approx([3.508636586, 4.991363414], abs=1e-9)
    assert m.field_events[0] < m.band_events[0] < m.band_events[1] < m.field_events[1]
    for t in m.field_events:
        assert m.axis_peak(t) == pytest.approx(.5, abs=1e-10)
    assert ContactNeck(Parameters(amplitude=.135)).band_events == [4.25]
    assert not ContactNeck(Parameters(amplitude=.134)).band_events
    assert len(ContactNeck(Parameters(amplitude=.136)).band_events) == 2
    assert m.gap(4.25) == pytest.approx(2*(.16-.18-.025))
    assert m.axis_peak(5.3) < .5
    for t in np.linspace(0,9,19):
        np.testing.assert_allclose(m.radius([-H,H],t), [1,1], atol=1e-14)


def test_density_is_even_bounded_and_fuses_without_switch():
    m = ContactNeck(); xs = np.linspace(-H,H,301)[None,:]; rho = np.linspace(-1.1,1.1,331)[:,None]
    for t in [0, m.field_events[0], m.tc, 4.25, m.field_events[-1], 9]:
        p = m.density(xs,rho,t)
        assert np.isfinite(p).all() and p.min() >= 0 and p.max() <= 1+1e-15
        np.testing.assert_allclose(p,m.density(xs,-rho,t),atol=1e-14)
    # The signed parameter crosses zero; the density formula remains well-defined.
    assert m.radius(0,4.25) < 0
    assert m.density(0,0,4.25) == pytest.approx(.7765475476329033)
    assert m.density(0,0,9) < .5
    # An event cache cannot alter the field: the visual is not a switched surface.
    before = m.density(xs,rho,4.)
    m.tc = None; m.field_events = []; m.band_events = []
    np.testing.assert_array_equal(before,m.density(xs,rho,4.))


@pytest.mark.parametrize('x,rho,t', [(-.129,.03,3.53),(0.,0.,4.25),(.2,.08,4.)])
def test_cartesian_hessian_matches_finite_differences(x,rho,t):
    m=ContactNeck(); step=2e-6
    def f(v):return float(m.density(v[0],np.hypot(v[1],v[2]),t))
    point=np.array([x,rho,0.]); eye=np.eye(3)*step; f0=f(point)
    diagonal=[(f(point+d)-2*f0+f(point-d))/step**2 for d in eye]
    mixed=(f(point+eye[0]+eye[1])-f(point+eye[0]-eye[1])-f(point-eye[0]+eye[1])+f(point-eye[0]-eye[1]))/(4*step**2)
    expected=m.hessian_components(x,rho,t)
    np.testing.assert_allclose([*diagonal,mixed],expected,rtol=2e-5,atol=.05)
    assert np.isfinite(expected).all()


def test_field_and_hessian_continuity_through_both_events():
    m=ContactNeck(); xs=np.linspace(-.23,.23,81)[None,:]; rho=np.linspace(0,.12,61)[:,None]
    for event in [*m.field_events,*m.band_events]:
        changes=[]
        for dt in [1e-3,1e-4,1e-5]:
            a=[m.density(xs,rho,event-dt),*m.hessian_components(xs,rho,event-dt)]
            b=[m.density(xs,rho,event+dt),*m.hessian_components(xs,rho,event+dt)]
            changes.append(max(float(np.max(abs(u-v))) for u,v in zip(a,b)))
        assert changes[2] < .11*changes[1] < .012*changes[0]


def test_metadata_and_html_controls():
    root=Path(__file__).resolve().parents[1]
    record=json.loads((root/'docs/images/finite_neck/contact.json').read_text())
    m=ContactNeck()
    assert record['regularization']=='C_SMOOTH_FIELD_DERIVED_SURFACE'
    assert record['field_contour_crossings']==pytest.approx(m.field_events,abs=1e-10)
    for s in record['samples']:
        assert s['psi']==pytest.approx(float(m.density(s['x'],s['rho'],s['t'])),abs=1e-14)
    node=shutil.which('node')
    if node is None:pytest.skip('Node is needed for the HTML math/control checks')
    subprocess.run([node,str(Path(__file__).with_name('finite_neck_ui.cjs'))],check=True)

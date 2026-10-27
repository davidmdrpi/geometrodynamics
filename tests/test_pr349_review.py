import numpy as np
import pytest
from geometrodynamics.waves.directional_flash_response import coasting_focus_audit
from geometrodynamics.waves.deformable_mode_conversion import DeformableDefect,Parameters


@pytest.mark.parametrize('k',[1,2])
def test_finite_aperture_capture_does_not_select_only_resonance(k):
    a=.3;n=2*k+1;upper=(np.pi-2*a)/(n*np.pi);lower=(np.pi-3*a)/(n*np.pi)
    for speed in [lower+.2*(upper-lower),lower+.8*(upper-lower)]:
        route=coasting_focus_audit(a,speed)
        hits=[e for e in route['events'] if e['target'] is not None]
        assert len(hits)==2 and all(e['target']==1-e['source'] for e in hits)
        assert all(e['return_number']==n and e['time']<route['contact_time'] for e in hits)
        t=np.linspace(0,route['contact_time'],37)
        np.testing.assert_allclose(np.exp(1j*(a+speed*t+np.pi)),
            np.exp(1j*(-a-speed*(route['contact_time']-t))),atol=2e-15)
    # The endpoint is contact and is excluded by the current first-flight rule.
    assert not any(e['target'] is not None for e in coasting_focus_audit(a,upper)['events'])


def test_moving_defects_do_not_source_the_zero_psi_channel():
    m=DeformableDefect(Parameters(length=80.,dx=.2));u,v=m.initial(pair=True,speed=.3)
    for _ in range(150):u,v=m.step(u,v,.04)
    np.testing.assert_array_equal(u[1:],0.)
    np.testing.assert_array_equal(v[1:],0.)

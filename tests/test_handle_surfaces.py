import json
import numpy as np
import pytest
from geometrodynamics.waves import handle_surfaces as s
from geometrodynamics.waves import handle_evolution as e
from experiments.closure_ledger import handle_surface_audit as a


def state(n=64):
    y=np.zeros((8,n));y[0]=1;y[1]=2
    return y


def test_flat_round_sphere_and_collapsing_control():
    y=state();dy=np.zeros_like(y);dy[1]=1
    np.testing.assert_allclose(s.expansions(y,dy)['Einstein'],np.tile([[1.],[-1.]],(1,64)))
    assert s.classify(s.expansions(y,dy)['Einstein'])['untrapped'].all()
    dy[1]=0;y[3]=.3
    assert s.classify(s.expansions(y,dy)['Einstein'])['future_trapped'].all()
    y[3]=-.3
    assert s.classify(s.expansions(y,dy)['Einstein'])['past_trapped'].all()
    y[3]=0
    labels=s.classify(s.expansions(y,dy)['Einstein'])
    assert labels['near_marginal'].all() and not labels['future_trapped'].any()


def test_conformal_time_derivative_can_change_trapping_class():
    y=state();dy=np.zeros_like(y)
    y[3]=-.01;y[4]=1;y[6]=-.3
    pair=s.expansions(y,dy)
    assert s.classify(pair['Einstein'])['past_trapped'].all()
    assert s.classify(pair['Jordan'])['future_trapped'].all()


def test_time_reversal_exchanges_and_negates_future_expansions():
    y=state();y[3]=.3;y[4]=.4;y[6]=.2
    dy=np.zeros_like(y);dy[1]=.07;dy[4]=.11
    reverse=y.copy();reverse[[2,3,6,7]]*=-1
    original=s.expansions(y,dy);back=s.expansions(reverse,dy)
    for frame in original:
        np.testing.assert_allclose(back[frame],-original[frame][::-1],atol=1e-15)
        np.testing.assert_allclose(s.expansions(y,-dy)[frame],original[frame][::-1],atol=1e-15)


def test_area_rate_route_against_analytic_inhomogeneous_conformal_geometry():
    n=128;x=2*np.pi*np.arange(n)/n
    y=state(n);y[0]=1+.1*np.cos(x);y[1]=2+.2*np.sin(x)
    y[3]=.1*np.cos(x);y[4]=.3*np.sin(x);y[5]=.2*np.cos(x)
    y[6]=.07*np.cos(x);y[7]=.05*np.sin(x)
    dy=np.zeros_like(y);dy[1]=.2*np.cos(x);dy[4]=.3*np.cos(x);dy[5]=-.2*np.sin(x)
    pair=s.expansions(y,dy);area=s.area_expansions(y,2*np.pi/n)
    for frame in pair:
        np.testing.assert_allclose(pair[frame],area[frame],atol=3e-10,rtol=0)


def test_cut_cancellation_is_orientation_not_independent_recoil():
    # Exact compact vacuum baseline C=1/2, psi=1: Q_X=4 pi.
    y=state();y[:2]=1;y[2]=1;y[3]=-.5
    assert s.section_flux(y,0)==pytest.approx(4*np.pi)
    assert s.section_flux(y,0,1)+s.section_flux(y,0,-1)==0
    # Independent sections need not cancel merely because their normals oppose.
    y[3,1]=-.25
    assert s.section_flux(y,0,1)+s.section_flux(y,1,-1)!=0


@pytest.mark.parametrize('damage',['nan','shape','f','metric'])
def test_invalid_surface_data_rejected(damage):
    y=state();dy=np.zeros_like(y)
    if damage=='nan':y[0,0]=np.nan
    elif damage=='shape':dy=dy[:7]
    elif damage=='f':y[4]=3
    else:y[0]=0
    with pytest.raises(ValueError):s.expansions(y,dy)


def test_archived_audit_without_any_evolution(monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError('audit must not rerun the evolution')
    monkeypatch.setattr(e,'run',forbidden)
    result=a.audit()
    archived=json.loads((a.DEFAULT/'audit.json').read_text())
    assert a.prior.agreement(result,archived)
    assert result['status']['FUTURE_TRAPPED_NECK_AT_FINAL_TIME']
    assert result['status']['EXTERIOR_BULK_RECIPROCAL_MOMENTUM_TRANSFER']=='NOT_TESTED'
    for row in result['rows']:
        for f in ('Einstein','Jordan'):
            initial=row['initial']['frames'][f]
            assert not initial['seam_future_trapped']
            assert initial['seam_past_trapped'] == (row['eta']==.3)
            final=row['final']['frames'][f]
            assert final['seam_future_trapped']
            assert final['area_route_max_difference']<3e-7
        assert not row['final']['cut_faces']['independent_surfaces']
    for row in result['convergence']:
        coarse,fine=row['seam_differences']
        assert 14<coarse/fine<18
        assert row['fine_negative_margin']>10
        assert row['margin_over_medium_fine']>1e9
    # Both frames find trapping at the neck, but differ at the central bulk.
    fine=result['rows'][-1]['final']['frames']
    assert min(fine['Einstein']['central_bulk_expansions'])>0
    assert max(fine['Jordan']['central_bulk_expansions'])<0


def test_tampered_archive_rejected(tmp_path):
    (tmp_path/'evolution.json').write_text('{"groups":[]}')
    with pytest.raises(ValueError,match='hash mismatch'):a.audit(tmp_path)


def test_incomplete_output_withdrawn(tmp_path,monkeypatch):
    def fail(*args,**kwargs):raise ValueError('damaged input')
    monkeypatch.setattr(a,'audit',fail)
    monkeypatch.setattr('sys.argv',['audit','--output-dir',str(tmp_path)])
    target=tmp_path/'audit.json';target.write_text('{"status":{"old_pass":true}}')
    with pytest.raises(ValueError,match='damaged input'):a.main()
    assert json.loads(target.read_text())=={'status':{},'error':'damaged input'}

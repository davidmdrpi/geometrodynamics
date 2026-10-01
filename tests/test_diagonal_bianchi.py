import numpy as np
from scipy.integrate import solve_ivp
from geometrodynamics.waves import diagonal_bianchi as d


def test_phase_map_jacobian_and_full_system():
    z=d.STAR+np.array([.001,.002,.003,.004,-.002,.001])
    image,jac=d.phase_map([z],True)
    h=1e-6
    f,_=d.phase_map(np.concatenate([z+h*np.eye(6),z-h*np.eye(6)]))
    np.testing.assert_allclose(jac[0],(f[:6]-f[6:]).T/(2*h),atol=1e-6,rtol=1e-6)
    hist=d.full_return(z)
    np.testing.assert_allclose(d.canonical(hist['returned']),image[0],atol=1e-10,rtol=0)
    for y in np.asarray(hist['samples'])[::32]:
        assert np.max(abs(d.full.constraints(y)['residual']))<1e-10


def test_diagonal_phase_vector_field_matches_full_system():
    z=d.STAR+np.array([.001,.002,.02,.03,-.01,.04])
    initial=d.full_initial(z)
    sol=solve_ivp(lambda t,y:d.full.conformal_rhs(y),(0,.2),initial,rtol=1e-12,atol=1e-14)
    y=sol.y[:,-1];A,u,q,qp,M,L=d.full.unpack(y)
    shape=d.E@np.log(np.diag(M))/2;velocity=d.E@np.diag(L)
    full=d.full.conformal_rhs(y)
    phase=np.arctan2(-qp[0]/2,q[0])
    clock=(qp[0]**2-q[0]*full[7])/(2*(q[0]**2+qp[0]**2/4))
    state=np.array([A,u,shape[0],velocity[0],shape[1],velocity[1],y[2]])
    accel=d.E@np.diag(d.full.unpack(full)[5])
    np.testing.assert_allclose(clock*d.phase_rhs(phase,state),[u,full[1],velocity[0],accel[0],velocity[1],accel[1],1.],atol=1e-11,rtol=1e-11)


def toy(points,jacobian=False):
    points=np.asarray(points);out=points.copy();J=np.zeros((len(points),6,6))
    out[:,0]=1+2*(points[:,0]-1);out[:,1]=points[:,1]/2
    J[:,0,0]=2;J[:,1,1]=.5
    for i,z in enumerate(points):
        v=z[2:];I=np.dot(v,v)/2;angle=2.9+.4*I
        R=np.array([[np.cos(angle),np.sin(angle)],[-np.sin(angle),np.cos(angle)]])
        R4=np.kron(np.eye(2),R);dv=np.kron(np.eye(2),np.array([[0.,1.],[-1.,0.]]))@R4@v
        out[i,2:]=R4@v
        J[i,2:,2:]=R4+.4*np.outer(dv,v)
    return (out,J) if jacobian else (out,np.ones(len(points)))


def test_toy_circular_family_has_known_action_and_rotation():
    t=2*np.pi*np.arange(63)/63;a=.05
    K=np.tile(d.STAR,(63,1));K[:,2:]=a*np.array([np.cos(t),-np.sin(t),-np.sin(t),-np.cos(t)]).T
    rec=d.circle(a,K,2.9,mapper=toy)
    assert 'error' not in rec
    assert abs(abs(d.action(K))-a*a)<1e-14
    assert abs(rec['omega']-2.9-.4*abs(d.action(rec['K'])))<1e-11
    assert d.circular_area(rec['K'])<0


def test_multiple_shooting_finds_known_nonzero_two_cycle(monkeypatch):
    def resonant(points,jacobian=False):
        points=np.asarray(points);out=points.copy();J=np.zeros((len(points),6,6))
        for i,z in enumerate(points):
            v=z[2:];angle=np.pi+.4*(v@v/2-.01)
            R=np.kron(np.eye(2),np.array([[np.cos(angle),np.sin(angle)],[-np.sin(angle),np.cos(angle)]]))
            dv=np.kron(np.eye(2),np.array([[0.,1.],[-1.,0.]]))@R@v
            out[i,:2]=[1+2*(z[0]-1),z[1]/2];out[i,2:]=R@v
            J[i,0,0]=2;J[i,1,1]=.5;J[i,2:,2:]=R+.4*np.outer(dv,v)
        return (out,J) if jacobian else (out,np.ones(len(points)))
    monkeypatch.setattr(d,'phase_map',resonant)
    z=d.STAR+np.array([0.,0.,.09,0.,0.,-.09])
    rec=d.shoot(np.array([z,2*d.STAR-z]))
    assert 'error' not in rec
    assert np.max(abs(rec['image']-rec['nodes'][::-1]))<1e-9
    assert abs(np.dot(rec['nodes'][0,2:],rec['nodes'][0,2:])/2-.01)<1e-8


import copy
import pytest
from experiments.closure_ledger import diagonal_bianchi_probe as probe


@pytest.fixture(scope='module')
def measured():
    return probe.load_run(probe.RUN)


def test_archive_replay_rebuilds_decisions():
    from experiments.closure_ledger.diagonal_bianchi_replay import replay
    assert replay()['replay']=='VERIFIED'


@pytest.mark.parametrize('damage',['image','missing_curve','nonfinite','missing_circle'])
def test_altered_circle_evidence_cannot_pass(measured,damage):
    raw=copy.deepcopy(measured)
    if damage=='image':
        raw['circles'][0]['image'][0][0]+=.001
        assert not probe.circle_metrics(raw['circles'][0])['accepted']
    else:
        if damage=='missing_curve':del raw['circles'][0]['K']
        elif damage=='nonfinite':raw['circles'][0]['K'][0][0]=float('nan')
        else:raw['circles'].pop(0)
        with pytest.raises((ValueError,KeyError)):
            probe.analyze(raw)


def test_changed_history_is_rejected(measured):
    rec=copy.deepcopy(measured['circles'][0])
    rec['full_histories'][0]['returned'][0]+=.01
    with pytest.raises(ValueError,match='endpoints'):
        probe.circle_metrics(rec)


def test_background_is_not_a_nontrivial_two_cycle():
    rec=dict(nodes=np.tile(d.STAR,(2,1)),image=np.tile(d.STAR,(2,1)))
    m=probe.closure_metrics(rec,.01)
    assert not m['candidate'] and not m['verified']


def test_changed_published_report_is_rejected(tmp_path):
    import shutil
    from experiments.closure_ledger.diagonal_bianchi_replay import replay
    for path in probe.RUN.iterdir():shutil.copy(path,tmp_path/path.name)
    with (tmp_path/'report.json').open('a') as stream:stream.write(' ')
    with pytest.raises(ValueError,match='fingerprint'):
        replay(tmp_path)

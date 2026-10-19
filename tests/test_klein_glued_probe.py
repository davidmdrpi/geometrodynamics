"""Independent energy assembly and a small actual glued-field run."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from klein_glued_probe import run,element_energy
from geometrodynamics.viz.klein_handles import build_surface,HandleSpec


def test_element_partition_includes_shared_rim_mass():
    s=build_surface(handles=[HandleSpec((1,0,0),(0,-1,0),True)],n_sphere=500,n_ring=16,ds_target=.025)
    K,M=s.operators();area,KT=s.local_stiffness()
    rng=np.random.default_rng(15);u=rng.normal(size=s.n_vertices);v=rng.normal(size=s.n_vertices)
    e=element_energy(s,area,KT,u,v);neck=s.tri_tube>=0
    expected=.5*(u@(K@u)+np.sum(M*(v*v+.25*u*u)))
    assert e[neck].sum()+e[~neck].sum()==pytest.approx(expected,rel=1e-13)
    assert e.min()>-1e-12


def test_fixed_glued_metric_does_not_trigger_default_contact():
    values,report=run(n_sphere=2000,glued=True)
    assert report['max_partition_residual']<1e-12
    assert report['max_relative_energy_drift']<1e-3
    assert report['max_neck_energy_fraction']>.65
    assert report['max_B_side_energy_fraction']>.5
    assert .35<report['waist_peak_over_reference']<.6
    assert not report['contact_detected']
    assert report['minimum_radius_in_rim_units']>report['contour_radius_threshold']+.1
    assert report['release_energy_injected']==0 and not report['topology_change_performed']

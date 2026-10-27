import numpy as np

from geometrodynamics.waves.deformable_mode_conversion import (
    DeformableDefect,Parameters,evolve,step_scattering,
)


def test_deck_map_preserves_potential_and_full_field_energy():
    model=DeformableDefect(Parameters(length=40.,dx=.2,twisted=True))
    rng=np.random.default_rng(3490)
    u=rng.normal(size=(3,model.n));v=rng.normal(size=u.shape)
    np.testing.assert_allclose(model.potential(model.deck(u)),model.potential(u),atol=1e-14)
    assert abs(model.energy(model.deck(u),model.deck(v))-model.energy(u,v))<1e-10
    np.testing.assert_allclose(model.acceleration(model.deck(u)),model.deck(model.acceleration(u)),atol=1e-13)


def test_acceleration_is_negative_gradient_of_discrete_energy():
    model=DeformableDefect(Parameters(length=40.,dx=.2,twisted=True))
    rng=np.random.default_rng(350)
    u=rng.normal(size=(3,model.n))*.3;v=np.zeros_like(u)
    direction=rng.normal(size=u.shape);h=2e-6
    derivative=(model.energy(u+h*direction,v)-model.energy(u-h*direction,v))/(2*h)
    expected=-model.par.dx*np.sum(model.acceleration(u)*direction)
    assert abs(derivative-expected)<2e-7


def test_scattering_ode_matches_tanh_formula_and_predicts_pull():
    for omega in [2.3,2.8,3.5]:
        r=step_scattering(omega)
        assert abs(r['reflection']-r['analytic_reflection'])<1e-10
        assert abs(r['flux_residual'])<2e-9
        assert r['impulse_per_energy']>0
    r=step_scattering(2.8,coupling=0.)
    assert abs(r['impulse_per_energy'])<2e-9


def test_whole_field_time_reversal_retraces_deformation():
    model=DeformableDefect(Parameters(length=80.,dx=.2,twisted=True))
    u0,v0=model.initial();model.add_packet(u0,v0,center=8.,width=2.,amplitude=.08)
    u,v,_=evolve(model,u0.copy(),v0.copy(),16.,.04,16.)
    recovered,reverse,_=evolve(model,u,-v,16.,.04,16.)
    np.testing.assert_allclose(recovered,u0,atol=2e-13)
    np.testing.assert_allclose(reverse,-v0,atol=2e-12)


def test_energy_error_improves_under_time_refinement():
    model=DeformableDefect(Parameters(length=80.,dx=.2,twisted=True))
    u0,v0=model.initial();model.add_packet(u0,v0,center=5.,width=1.5,amplitude=.1)
    errors=[]
    for dt in [.08,.04]:
        _,_,rows=evolve(model,u0.copy(),v0.copy(),12.,dt,.4)
        errors.append(max(abs(r['energy_error']) for r in rows))
    assert errors[1]<errors[0]/3.5

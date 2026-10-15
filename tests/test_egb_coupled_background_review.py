"""No-shell reference signature and physical density in the shrinking-mouth limit."""
import hashlib
import json
from pathlib import Path
import pytest
from scipy.special import eval_gegenbauer
from experiments.closure_ledger.egb_coupled_background_review import no_shell_mode,odd_mouth_control


def test_regular_sphere_spectrum_and_norm_by_quadrature():
    for ell in range(7):
        mode=no_shell_mode(ell)
        assert mode['full_sphere_L2']==pytest.approx(1.,abs=2e-12)
        assert mode['kinetic_norm_quadrature']==pytest.approx(mode['kinetic_norm'],abs=2e-10)
        assert mode['gauge']==(ell==1)
        # Gegenbauer differential equation evaluated using independent
        # polynomial derivatives in x=cos(u), not a shooting solution.
        for x in [-.7,.1,.6]:
            f=eval_gegenbauer(ell,1.5,x)
            df=0. if ell==0 else 3*eval_gegenbauer(ell-1,2.5,x)
            ddf=0. if ell<2 else 15*eval_gegenbauer(ell-2,3.5,x)
            assert (1-x*x)*ddf-4*x*df+ell*(ell+3)*f==pytest.approx(0.,abs=1e-10)
    assert no_shell_mode(0)['kinetic_norm']==pytest.approx(-9.12) # full-sphere L2=1
    assert no_shell_mode(1)['kinetic_norm']==0


def test_background_growth_is_universal_but_one_growing_mode_is_not():
    for C in [.1,.2,.76,2.]:
        assert no_shell_mode(0,C)['growing_physical']
    assert sum(no_shell_mode(l,.76)['growing_physical'] for l in range(7))==1
    assert no_shell_mode(2,.2)['growing_physical']
    assert no_shell_mode(2,.76)['omega2_R2']>0
    assert not no_shell_mode(1,.2)['growing_physical'] # negative Lambda but gauge
    assert no_shell_mode(2,R=2)['kinetic_norm']==pytest.approx(16*no_shell_mode(2)['kinetic_norm'])


def test_odd_mode_density_and_norm_fade_at_fixed_shape_normalization():
    rows=[odd_mouth_control(a) for a in [.3,.1,.02]]
    for row in rows:
        assert row['omega2_R2']<0 and row['h']>0
        assert row['density_contrast_quadrature']==pytest.approx(row['density_contrast_L2_half_bulk'],abs=1e-10)
        assert row['norm_gradient_quadrature']==pytest.approx(row['norm_half_bulk_L2_unit'],abs=2e-8)
    baseline=no_shell_mode(0)['omega2_R2']
    assert abs(rows[-1]['even_lowest_omega2_R2']-baseline)<.002
    assert abs(rows[-1]['even_lowest_omega2_R2']-baseline)<abs(rows[0]['even_lowest_omega2_R2']-baseline)
    assert rows[-1]['h']<.0061
    assert rows[-1]['density_contrast_L2_half_bulk']<rows[0]['density_contrast_L2_half_bulk']/200
    assert rows[-1]['growth_rate_R']>rows[0]['growth_rate_R'] # rate alone misleads


def test_review_archive_splits_reasons_without_changing_original_evidence():
    root=Path(__file__).resolve().parents[1]
    result=json.loads((root/'experiments/closure_ledger/runs/20261014_egb_coupled_review/result.json').read_text())
    assert result['label']=='BACKGROUND_AND_SHELL_FAILURES_SEPARATED'
    assert all(result['numerical_gates'].values())
    assert result['original_label']=='COUPLED_RADIAL_ADMISSION_FAILED'
    assert set(result['failure_reasons'])=={'background','shell'}
    assert len(result['failure_reasons']['shell'])==3
    c=result['comparison']
    assert c['with_shell_negative_even_directions']-c['no_shell_negative_even_directions']==1
    assert c['with_shell_growing_found']-c['no_shell_growing_physical']==2
    assert hashlib.sha256((root/result['original_archive']).read_bytes()).hexdigest()==result['original_archive_sha256']
    for path,digest in result['sources'].items():
        assert hashlib.sha256((root/path).read_bytes()).hexdigest()==digest

import copy
import json
import shutil
import pytest
from experiments.closure_ledger import r3_diagonal_obstruction_probe as probe
from experiments.closure_ledger import r3_diagonal_obstruction_replay as replay


def rows():return json.loads((probe.RUN/'scan.json').read_text())['rows']


def test_authenticated_score_and_provenance():
    out=replay.replay()
    assert out['replay']=='VERIFIED' and out['phases']==192
    assert out['label']=='NUMERICALLY_UNRESOLVED' and out['roots']==0
    assert out['supplemental_roots']==24


@pytest.mark.parametrize('damage',['truncate','duplicate','phase','slice','nonfinite','residual_flag','error_estimate'])
def test_modified_numerical_records_rejected(damage):
    r=copy.deepcopy(rows())
    if damage=='truncate':r.pop()
    elif damage=='duplicate':r[1]=r[0]
    elif damage=='phase':r[0]['phi']+=.01
    elif damage=='slice':r[0]['point']['c0'][2]+=.001
    elif damage=='nonfinite':r[0]['check']['nodes'][0][2]='nan'
    elif damage=='residual_flag':r[0]['check']['residual']='1e-10'
    elif damage=='error_estimate':r[0]['error_estimate']='1e-5'
    with pytest.raises(ValueError):replay.validate_scan(r)


def test_archive_byte_change_rejected(tmp_path):
    for name in replay.FILES+('manifest.json',):shutil.copy(probe.RUN/name,tmp_path/name)
    p=tmp_path/'scan.json';p.write_text(p.read_text()+' ')
    with pytest.raises(ValueError,match='fingerprint'):replay.replay(tmp_path)


def test_score_does_not_turn_failed_accuracy_into_a_null():
    r=copy.deepcopy(rows());r[0]['error_estimate']='1e-20'
    s=probe.scan_score(r,dict(ok=True));assert s['label']=='NUMERICALLY_UNRESOLVED'


@pytest.mark.parametrize('damage',['truncate','duplicate','isolation_flag','closure','slice_residual','nonfinite'])
def test_modified_root_records_rejected(damage):
    rr=json.loads((probe.RUN/'supplement.json').read_text())['rows']
    brackets=json.loads((probe.RUN/'scan_result.json').read_text())['brackets']
    if damage=='truncate':rr.pop()
    elif damage=='duplicate':rr[1]=copy.deepcopy(rr[0])
    elif damage=='isolation_flag':rr[0]['diagnostics']['ok']=not rr[0]['diagnostics']['ok']
    elif damage=='closure':rr[0]['diagnostics']['closure_residual']='1e-10'
    elif damage=='slice_residual':rr[0]['history'][-1]['residual']='1e-10'
    elif damage=='nonfinite':rr[0]['nodes'][0][0]='nan'
    with pytest.raises(ValueError):replay.validate_roots(rr,brackets)

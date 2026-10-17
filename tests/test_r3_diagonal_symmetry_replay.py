import copy
import json
import shutil
import numpy as np
import pytest
from experiments.closure_ledger import r3_diagonal_symmetry_probe as probe
from experiments.closure_ledger import r3_diagonal_symmetry_replay as replay


def records():
    return json.loads((probe.RUN/'actions.json').read_text())['records']


def test_authenticated_archive_rescores():
    assert replay.replay()['label']=='ORDER_TWELVE_SELECTION_SUPPORTED'


@pytest.mark.parametrize('damage',['truncate','duplicate','move_input','nonfinite','radau_grid'])
def test_structural_damage_is_rejected(damage):
    rows=copy.deepcopy(records())
    if damage=='truncate': rows.pop()
    elif damage=='duplicate': rows[1]=rows[0]
    elif damage=='move_input': rows[0]['node'][2]+=.001
    elif damage=='nonfinite': rows[0]['chain'][2][3]=float('nan')
    elif damage=='radau_grid': rows[1]['radau_half']=rows[0]['radau_half']
    with pytest.raises(ValueError):replay.validate(rows)


@pytest.mark.parametrize('key,label',[
    ('half_cycle','PROPOSED_GROUP_ACTION_FAILED'),
    ('matrix_half','NUMERICALLY_UNRESOLVED')])
def test_registered_label_can_fail(key,label):
    rows=copy.deepcopy(records());rows[0][key]['z'][2]+=.001
    result=probe.score(rows,probe.data()[:,:6])
    assert result['label']==label


def test_modified_archive_is_rejected(tmp_path):
    for name in ('manifest.json','actions.json','result.json'):
        shutil.copy(probe.RUN/name,tmp_path/name)
    path=tmp_path/'actions.json';path.write_text(path.read_text()+' ')
    with pytest.raises(ValueError,match='fingerprint'):replay.replay(tmp_path)

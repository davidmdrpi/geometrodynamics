"""Fast authenticity/decision regressions; CI separately reconstructs all fields."""
import json
import shutil
import numpy as np
import pytest
from experiments.closure_ledger import aperture_phase_probe as p

MANIFEST_SHA='02b8bd392282292cbf0f82a42dfdf17ba2263a84b709bebf38217fa394037a6d'


def test_phase_evidence_pins_inventory_sources_and_both_decisions():
    assert p.digest(p.RUN/'manifest.json')==MANIFEST_SHA
    manifest=json.loads((p.RUN/'manifest.json').read_text())
    assert set(manifest)=={n+'.npz.b64' for n,c in p.schedule()}|{'started.json','provenance.json','result.json'}
    for path,sha in manifest.items():assert p.digest(p.RUN/path)==sha
    provenance=json.loads((p.RUN/'provenance.json').read_text())
    assert provenance['sources']==p.sources()
    assert provenance['started_utc']<provenance['finished_utc']
    saved=json.loads((p.RUN/'result.json').read_text())
    scored=p.score(saved['diagnostics'],saved['comparisons'])
    assert scored==saved
    assert saved['numerical_validity'] and saved['phase_coverage']
    assert saved['retention_label']=='FAILED_IN_DECLARED_FAMILY'
    assert saved['inverse_phase_label']=='SUPPORTED_IN_DECLARED_FAMILY'


def test_ulp_archive_mutation_fails_authentication_before_numeric_tolerance(tmp_path):
    shutil.copyfile(p.RUN/'manifest.json',tmp_path/'manifest.json')
    filename=p.schedule()[0][0]+'.npz.b64'
    r=p.read(p.RUN/filename)
    index=np.argmax(abs(r['incoming'][:,0]));value=r['incoming'][index,0]
    r['incoming'][index,0]=np.nextafter(value,np.inf)
    assert abs(r['incoming'][index,0]-value)<1e-14
    p.save(tmp_path/filename,r)
    with pytest.raises(ValueError,match='archive mismatch'):p.replay(tmp_path,MANIFEST_SHA)


def test_phase_replay_rejects_missing_evidence(tmp_path):
    shutil.copyfile(p.RUN/'manifest.json',tmp_path/'manifest.json')
    with pytest.raises(FileNotFoundError):p.replay(tmp_path,MANIFEST_SHA)
    with pytest.raises(ValueError,match='manifest mismatch'):p.replay(tmp_path,'0'*64)

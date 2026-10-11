import shutil

import pytest

from experiments.closure_ledger import r3_half_clock_probe as probe
from experiments.closure_ledger import r3_half_clock_replay as replay


def test_authentication_and_complete_rescore():
    result = replay.replay()
    assert result['replay'] == 'VERIFIED'
    assert result['freeze_commit'] == replay.FREEZE


def test_changed_archive_bytes_fail_closed(tmp_path):
    for name in (*replay.FILES, 'manifest.json'):
        shutil.copy(probe.DIRECTORY/name, tmp_path/name)
    with (tmp_path/'scan.json').open('a') as stream:
        stream.write('\n')
    with pytest.raises(ValueError, match='archive fingerprint mismatch: scan.json'):
        replay.replay(tmp_path)


def test_changed_source_bindings_fail_closed(monkeypatch):
    original = probe.bindings()
    original['geometrodynamics/waves/r3_half_clock.py'] = '0'*64
    monkeypatch.setattr(probe, 'bindings', lambda: original)
    with pytest.raises(ValueError, match='source or input fingerprint mismatch'):
        replay.replay()

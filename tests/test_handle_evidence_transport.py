import hashlib
import shutil
from pathlib import Path
import pytest
from experiments.closure_ledger.restore_handle_evidence import DEFAULT, RAW_SHA256, RAW_SIZE, restore


@pytest.fixture
def transport(tmp_path):
    for source in DEFAULT.glob('evolution.*'):
        if source.name != 'evolution.json':
            shutil.copyfile(source, tmp_path/source.name)
    return tmp_path


def test_lossless_restore(transport):
    target = restore(transport)
    assert target.stat().st_size == RAW_SIZE
    assert hashlib.sha256(target.read_bytes()).hexdigest() == RAW_SHA256
    assert restore(transport) == target


def test_corrupted_part_is_rejected(transport):
    part = transport/'evolution.json.gz.b64.part000'
    content = bytearray(part.read_bytes())
    content[20] ^= 1
    part.write_bytes(content)
    with pytest.raises(ValueError, match='part hash mismatch'):
        restore(transport)
    assert not (transport/'evolution.json').exists()


def test_missing_part_is_rejected(transport):
    (transport/'evolution.json.gz.b64.part000').unlink()
    with pytest.raises(FileNotFoundError):
        restore(transport)


def test_existing_corrupt_evidence_is_not_silently_overwritten(transport):
    (transport/'evolution.json').write_text('{}')
    with pytest.raises(ValueError, match='existing evolution evidence hash mismatch'):
        restore(transport)

"""Admission checks for the prescribed contact geometry, not fluid dynamics."""
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from scripts.neck_wave_contact import ContactNeck, Parameters, H, END, X


@pytest.mark.parametrize("amplitude", [0., .08])
def test_null_controls_do_not_reconnect(amplitude):
    model = ContactNeck(Parameters(amplitude=amplitude))
    assert model.tc is None
    for t in np.linspace(0, END, 37):
        assert model.gap(t) > 0
        assert len(model.branches(t)) == 1
        if amplitude == 0:
            np.testing.assert_array_equal(model.radius(X, t), model.rest(X))


def test_first_contact_is_driven_by_wave_and_preserves_rims():
    model = ContactNeck()
    assert 0 < model.tc < END
    lo, hi = model.bracket
    assert lo <= model.tc <= hi
    assert model.gap(lo) > 0 >= model.gap(hi)
    assert abs(model.gap(model.tc)) < 1e-10
    assert model.gap(model.tc - 1e-3) > 0
    assert model.xc < 0  # incoming crest reaches contact before the waist center
    for t in np.linspace(0, model.tc, 101):
        assert np.min(model.radius(X, t)) > 0
        np.testing.assert_allclose(model.radius([-H, H], t), [1, 1], atol=1e-14)
    assert ContactNeck(Parameters(amplitude=.22)).tc < model.tc
    assert ContactNeck(Parameters(waist=.30)).tc is None


def test_reconnection_has_separated_closed_rounded_tips():
    model = ContactNeck()
    for t in np.linspace(model.tc, END, 31):
        left, right = model.branches(t, n=501)
        assert left[0][-1] < right[0][-1]
        for x, r, _, _ in [left, right]:
            assert r[0] == 1 and r[-1] == 0
            assert np.isfinite(r).all() and (r[:-1] > 0).all()
            # A rounded pole has r proportional to sqrt(distance), not distance.
            # Sample sufficiently close to the tip to resolve that asymptote.
            fine = model.branches(t, n=10001)[0 if x[0] < 0 else 1]
            distance = abs(fine[0] - fine[0][-1])
            exponent = np.log(fine[1][-3] / fine[1][-2]) / np.log(distance[-3] / distance[-2])
            assert .48 < exponent < .52


def test_recorded_contact_matches_current_model():
    root = Path(__file__).resolve().parents[1]
    recorded = json.loads((root / "docs/images/finite_neck/contact.json").read_text())
    model = ContactNeck()
    assert recorded["contact_time"] == pytest.approx(model.tc, abs=1e-12)


def test_interactive_model_and_controls():
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is needed for the HTML math/control checks")
    subprocess.run([node, str(Path(__file__).with_name("finite_neck_ui.cjs"))], check=True)

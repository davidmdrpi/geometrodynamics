"""CI run of the absolute-coupling capstone probe (PR #225).

The capstone is the arc's integration test: the canonical Hopf-KK
chain, the Einstein-frame radion, the modulus-rank audit, and the
alpha-dependent holdout over the committed run ledgers.  Running it
in CI keeps the ledger honest: any drift in the committed artifacts
or in the symbolic/numeric machinery fails the suite.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]
                       / "experiments" / "closure_ledger"))

import absolute_coupling_capstone_probe as probe  # noqa: E402


@pytest.fixture(scope="module")
def summary():
    return probe.run_probe()


def test_capstone_all_green(summary):
    failed = [t["name"] for t in summary["tests"] if not t["pass"]]
    assert not failed, f"capstone tests failed: {failed}"
    assert summary["n_passed"] == summary["n_total"] == 9


def test_capstone_verdict_established(summary):
    assert "INCONCLUSIVE" not in summary["verdict_class"]
    assert summary["verdict"].startswith("ESTABLISHED")


def test_alpha_holdout_one_common_alpha(summary):
    t7 = [t for t in summary["tests"] if t["name"] == "T7_holdout"][0]
    assert t7["n_checks"] == 12
    assert t7["alpha_inferred_relative_spread"] < 1e-12
    assert t7["worst_rel_dev"] < 1e-9


def test_rank_audit_cap_fixes_alpha_direction(summary):
    t5 = [t for t in summary["tests"] if t["name"] == "T5_answer"][0]
    ra = t5["rank_audit"]
    assert ra["rank_before"] == 1 and ra["flats_before"] == 4
    assert ra["grad_alpha_null_projection_before"] > 1.0
    assert ra["rank_after"] == 2 and ra["flats_after"] == 3
    assert ra["grad_alpha_null_projection_after"] < 1e-12

"""Smoke tests for lepton calibration CLI utilities."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_sweep_k_uplift_beta_cli_smoke():
    env = dict(os.environ)
    env["PYTHONPATH"] = str(REPO_ROOT)
    cmd = [
        sys.executable,
        "scripts/sweep_k_uplift_beta.py",
        "--beta-center", "72",
        "--beta-span", "4",
        "--beta-steps", "3",
        "--phase-steps", "3",
        "--transport-steps", "3",
        "--pinhole-steps", "3",
        "--n-points", "16",
        "--top-k", "1",
    ]
    proc = subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=True,
    )
    out = proc.stdout
    assert "Tau target" in out
    assert "action_base (locked)" in out
    assert "Betas with exact mu/e roots" in out
    if "Betas with exact mu/e roots = 0" not in out:
        assert "Best beta" in out


def test_refine_locked_tau_cli_smoke():
    env = dict(os.environ)
    env["PYTHONPATH"] = str(REPO_ROOT)
    cmd = [
        sys.executable,
        "scripts/refine_locked_tau.py",
        "--beta-min", "138",
        "--beta-max", "142",
        "--beta-step", "2",
        "--phase-steps", "3",
        "--transport-steps", "3",
        "--pinhole-steps", "3",
        "--n-points", "16",
        "--top-k", "1",
    ]
    proc = subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=True,
    )
    out = proc.stdout
    assert "geometric beta estimate" in out
    assert "Betas with exact mu/e roots in local window" in out
    assert "Integer-anchored geometric beta family check" in out


def test_lock_beta_50pi_probe_cli_smoke():
    env = dict(os.environ)
    env["PYTHONPATH"] = str(REPO_ROOT)
    cmd = [
        sys.executable,
        "scripts/lock_beta_50pi_probe.py",
        "--n-points", "16",
    ]
    proc = subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=True,
    )
    out = proc.stdout
    assert "Locked beta probe" in out
    assert "uplift_2pi_quanta" in out

"""Shared fixtures for the geometrodynamics test suite."""

import pytest


@pytest.fixture
def l_values():
    """Default angular momentum values for radial mode tests."""
    return [1, 3, 5]

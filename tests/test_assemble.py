"""Smoke tests for Assemble.pyx — verify module import and basic functionality."""

import pytest


def test_assemble_module_imports():
    """Verify the Assemble module imports successfully."""
    from LOCATE.Assemble import assemble_cluster
    assert assemble_cluster is not None


def test_assemble_module_has_expected_attributes():
    """Smoke test: check expected attributes exist on Assemble module."""
    import LOCATE.Assemble as mod
    attrs = dir(mod)
    assert 'assemble_cluster' in attrs

"""Integration-level smoke tests for the LOCATE pipeline."""

import pytest
import subprocess
import os


def test_module_imports():
    """Verify all extension modules import successfully."""
    from LOCATE import Main
    from LOCATE import FileIO
    from LOCATE import Cluster
    from LOCATE import Assemble
    from LOCATE import Annotate
    from LOCATE import ParallelTaskExecutor
    assert True

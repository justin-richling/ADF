# lib/test/conftest.py
import pytest
from lib.adf_diag import AdfDiag
from pathlib import Path
import shutil

@pytest.fixture
def diag_fixture(tmp_path):
    """
    Return a minimal AdfDiag object pointing to temporary directories
    so that diag.create_climo() writes to tmp_path instead of real paths.
    """
    config_yaml = Path("lib/test/test_files/config_simple_test.yaml")
    diag = AdfDiag(config_yaml, debug=True)

    # Override input/output locations to temp path
    diag.diag_cam_climo.cam_ts_loc = tmp_path / "ts"
    diag.diag_cam_climo.cam_climo_loc = tmp_path / "climo"
    diag.diag_cam_climo.cam_ts_loc.mkdir()
    diag.diag_cam_climo.cam_climo_loc.mkdir()

    # Copy tiny sample dataset
    sample_file = Path("lib/test/unit_tests/output/climo/baseline_climo.nc")
    #baseline_file_loc = Path("lib/test/unit_tests/output/climo/")
    shutil.copy(sample_file, diag.diag_cam_climo.cam_climo_loc)

    # Minimal climo years
    #diag.climo_yrs = {"syears": [2000], "eyears": [2001]}

    return diag

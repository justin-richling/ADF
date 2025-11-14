# lib/test/test_climo_pipeline.py
import xarray as xr
from pathlib import Path

def test_create_climo(diag_fixture):
    """
    Run diag.create_climo() on sample data and compare to baseline.
    """
    diag_fixture.create_climo()

    # Check that output files exist
    out_file = next(diag_fixture.diag_cam_climo.cam_climo_loc.glob("*.nc"))
    assert out_file.exists()

    # Compare NetCDF data to baseline
    baseline_file_loc = Path("lib/test/unit_tests/output/climo/")
    #Path("lib/test/unit_tests/output/climo/baseline_climo.nc")
    #Check if plot output directory exists, and if not, then create it:
    if not baseline_file_loc.is_dir():
        print(f"    {baseline_file_loc} not found, making new directory")
        baseline_file_loc.mkdir(parents=True)
    new_climo = xr.open_dataset(out_file)
    baseline = xr.open_dataset(baseline_file_loc / "baseline_climo.nc")

    xr.testing.assert_allclose(new_climo, baseline, rtol=1e-6, atol=1e-8)

# tests/climo_tests/test_climo_generation.py
import xarray as xr
from lib.adf_diag import create_climo  # replace with your function

def test_climo_generation(tmp_path):
    input_file = "tests/data/sample_input.nc"
    baseline_file = "tests/data/baseline_climo.nc"
    output_file = tmp_path / "climo_out.nc"

    # Generate climatology
    create_climo(input_file, output_file)

    # Compare with golden baseline
    new = xr.open_dataset(output_file)
    baseline = xr.open_dataset(baseline_file)

    xr.testing.assert_allclose(new, baseline, rtol=1e-6, atol=1e-8)

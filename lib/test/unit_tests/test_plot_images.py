import unittest
import os
import matplotlib.pyplot as plt
from pathlib import Path
from lib.plotting.global_latlon_map import global_latlon_map  # adjust to your actual import
from lib.test.unit_tests.image_utils import compare_images

class TestPlotImages(unittest.TestCase):
    BASE_DIR = Path(__file__).parent
    BASELINE = BASE_DIR / "baseline_images"
    OUTPUT = BASE_DIR / "output_images"

    def setUp(self):
        self.OUTPUT.mkdir(exist_ok=True)

    def test_global_latlon_plot(self):
        """Check that global_latlon_map plot matches baseline."""
        fig, ax = plt.subplots()
        # run your plotting function (replace with your own)
        global_latlon_map(ax=ax)
        out_path = self.OUTPUT / "test_global_latlon_map.png"
        fig.savefig(out_path, dpi=100)
        plt.close(fig)

        baseline_path = self.BASELINE / "test_global_latlon_map.png"
        self.assertTrue(
            compare_images(out_path, baseline_path, tol=5),
            msg="Global lat-lon plot differs from baseline image!"
        )

import numpy as np
from PIL import Image

def compare_images(img1_path, img2_path, tol=5):
    """
    Compare two images pixel-wise and return True if they are similar.
    `tol` is the mean absolute difference threshold (0–255 scale).
    """
    img1 = np.array(Image.open(img1_path).convert("RGB"), dtype=np.int16)
    img2 = np.array(Image.open(img2_path).convert("RGB"), dtype=np.int16)

    if img1.shape != img2.shape:
        return False

    diff = np.abs(img1 - img2)
    mean_diff = np.mean(diff)
    return mean_diff <= tol

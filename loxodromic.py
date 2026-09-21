"""Generate the loxodromic surface used in the accompanying writeup.

From this directory, install the adjacent mesh utility repository once:

    pip install -e ../mesh-utils

Then run `python loxodromic.py`.  The generated ASCII STL is written to
`meshes/` and includes this source file as a comment.
"""

from pathlib import Path

import numpy as np

from mesh_utils import save_solid


DPI = 80
XMAX, YMAX = 16.5, 35.5
XSTEPS, YSTEPS = int(XMAX * DPI), int(YMAX * DPI)
BOTTOM_Z = -1
ZOOM = 5
Z_SCALE = 0.15

# Sample more densely around the two singular points, w = ±i.
X_CENTERS = [(XMAX / 2, 2.0)]
Y_CENTERS = [(YMAX / 2 + ZOOM, 3.0), (YMAX / 2 - ZOOM, 3.0)]
FOCUS_STRENGTH = 30.0


def cayley_transform(w):
    """Map the upper half-plane coordinate `w` to the unit disk."""
    return (w + 1j) / (w - 1j)


GAMMA = 1.9 * np.exp(1j * np.pi / 3)
LOG_GAMMA = np.log(GAMMA)
A = np.real(1 / LOG_GAMMA)
B = np.imag(1 / LOG_GAMMA)


def inverse_loxodromic_coordinate(z):
    """Apply the logarithmic coordinate associated with `GAMMA`."""
    return np.log(z) / LOG_GAMMA


def surface(x, y):
    """Return surface heights and optional STL comment text for the x/y grid."""
    w = (x - XMAX / 2) / ZOOM + 1j * (y - YMAX / 2) / ZOOM
    coordinates = inverse_loxodromic_coordinate(cayley_transform(w))
    horizontal = np.real(coordinates)
    vertical = np.imag(coordinates)

    height = np.sin(8 * vertical / A) + 0.5 * np.sin(12 * horizontal / B)
    distance_to_positive_i = 3 * np.abs(w - 1j)
    distance_to_negative_i = 3 * np.abs(w + 1j)
    damping_positive_i = 1 - 0.95 / (1 + distance_to_positive_i**2)
    damping_negative_i = 1 - 0.95 / (1 + distance_to_negative_i**2)
    return Z_SCALE * height * damping_positive_i * damping_negative_i, ""


if __name__ == "__main__":
    save_solid(
        surface,
        XMAX,
        YMAX,
        XSTEPS,
        YSTEPS,
        BOTTOM_Z,
        Path(__file__).name,
        x_centers=X_CENTERS,
        y_centers=Y_CENTERS,
        focus_strength=FOCUS_STRENGTH,
    )

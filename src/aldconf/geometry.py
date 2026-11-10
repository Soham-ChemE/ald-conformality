"""Feature geometries. Lengths in metres.

Both shapes expose, as functions of the local opening (which shrinks as the film grows):
- hydraulic diameter d_h = 4 A / P (sets the Knudsen diffusivity),
- surface-to-volume ratio P / A (sets the wall sink per unit gas volume).
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class Slit:
    """Lateral channel of gap `gap` (height) and length `length`, width >> gap (PillarHall-type)."""
    gap: float
    length: float

    @property
    def aspect_ratio(self):
        return self.length / self.gap

    def hydraulic_diameter(self, opening):
        return 2.0 * np.asarray(opening)

    def surface_to_volume(self, opening):
        return 2.0 / np.asarray(opening)

    def walls(self):
        # number of growing walls that reduce the opening (top and bottom)
        return 2


@dataclass
class Hole:
    """Cylindrical hole of diameter `diameter` and depth `length` (3D NAND channel hole type)."""
    diameter: float
    length: float

    @property
    def gap(self):
        return self.diameter

    @property
    def aspect_ratio(self):
        return self.length / self.diameter

    def hydraulic_diameter(self, opening):
        return np.asarray(opening)

    def surface_to_volume(self, opening):
        return 4.0 / np.asarray(opening)

    def walls(self):
        # film grows radially inward: diameter shrinks by 2 x thickness
        return 2

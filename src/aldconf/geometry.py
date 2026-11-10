"""Feature geometries. Lengths in metres.

All shapes expose, as functions of the local opening (which shrinks as the film grows):
- hydraulic diameter d_h = 4 A / P (sets the Knudsen diffusivity),
- surface-to-volume ratio P / A (sets the wall sink per unit gas volume),
- cross-section A (per unit width for a slit), used by the solver when the opening varies along the feature,
- initial_opening(x): the opening before any film, uniform for Slit and Hole, linear in depth for TaperedHole.
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

    def cross_section(self, opening):
        # per unit width
        return np.asarray(opening, dtype=float)

    def initial_opening(self, x):
        return np.full(len(x), self.gap, dtype=float)

    @property
    def tapered(self):
        return False

    def wall_area(self):
        """Wall area per unit width (both faces), ignoring the closed end."""
        return 2.0 * self.length

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

    def cross_section(self, opening):
        return np.pi * np.asarray(opening, dtype=float) ** 2 / 4.0

    def initial_opening(self, x):
        return np.full(len(x), self.diameter, dtype=float)

    @property
    def tapered(self):
        return False

    def wall_area(self):
        """Side-wall area of the cylinder (the bottom disc is not counted)."""
        return np.pi * self.diameter * self.length

    def walls(self):
        # film grows radially inward: diameter shrinks by 2 x thickness
        return 2


@dataclass
class TaperedHole(Hole):
    """Hole whose diameter falls linearly from `diameter` at the mouth to `bottom_ratio * diameter` at the bottom
    (a frustum). bottom_ratio = 1 is the straight Hole; the aspect ratio is defined on the mouth diameter.
    The solver uses the area-weighted (frustum) form of the diffusion equation for this geometry."""
    bottom_ratio: float = 1.0

    def __post_init__(self):
        if not (0.0 < self.bottom_ratio <= 1.0):
            raise ValueError("bottom_ratio must be in (0, 1]")

    def diameter_at(self, x):
        x = np.asarray(x, dtype=float)
        return self.diameter * (1.0 - (1.0 - self.bottom_ratio) * x / self.length)

    def initial_opening(self, x):
        return self.diameter_at(x)

    @property
    def tapered(self):
        return self.bottom_ratio < 1.0

    def wall_area(self):
        """Lateral area of the frustum pi (r_top + r_bot) sqrt(L^2 + (r_top - r_bot)^2)."""
        r_t = self.diameter / 2.0
        r_b = self.bottom_ratio * self.diameter / 2.0
        return np.pi * (r_t + r_b) * np.sqrt(self.length ** 2 + (r_t - r_b) ** 2)

    def area_ratio_to_straight(self):
        """Wall area relative to the straight hole of the mouth diameter."""
        return self.wall_area() / (np.pi * self.diameter * self.length)

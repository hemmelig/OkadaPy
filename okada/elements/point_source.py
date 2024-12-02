"""

:copyright:
    2024, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

from dataclasses import dataclass


@dataclass
class PointShear:
    """Class to encapsulate point source capable of shearing."""

    x_start: float
    x_end: float
    y_start: float
    y_end: float
    z_start: float
    z_end: float
    dip_angle: float
    right_lateral_potency: float = 0.0
    reverse_potency: float = 0.0

    def __str__(self):
        """Returns a summary of the point shear source."""

        out = "Point source of shear:\n"
        out += f"   X: {(self.x_end + self.x_start) / 2:5.3f} / km\n"
        out += f"   Y: {(self.y_end + self.y_start) / 2:5.3f} / km\n"
        out += f"   Z: {(self.z_end + self.z_start) / 2:5.3f} / km\n"
        out += f" Dip: {self.dip_angle:5.3f}\n"
        out += "Components of shear:\n"
        out += f" Right-lateral: {self.right_lateral_potency:5.3f} / m\n"
        out += f"       Reverse: {self.reverse_potency:5.3f} / m\n"

        return out

    @property
    def raw_input(self) -> list:
        """Prepare the raw input as a 10 element wide list."""

        raw_input = [
            self.x_start,
            self.y_start,
            self.x_end,
            self.y_end,
            400,
            self.right_lateral_potency,
            self.reverse_potency,
            self.dip_angle,
            self.z_start,
            self.z_end,
        ]

        return raw_input


@dataclass
class PointInflation:
    """Class to encapsulate a point source capable of inflation/deflation."""

    x_start: float
    x_end: float
    y_start: float
    y_end: float
    z_start: float
    z_end: float
    dip_angle: float
    tensile_opening: float = 0.0
    point_opening: float = 0.0

    def __str__(self):
        """Returns a summary of the point inflation source."""

        out = "Point source of inflation:\n"
        out += f"   X: {(self.x_end + self.x_start) / 2:5.3f} / km\n"
        out += f"   Y: {(self.y_end + self.y_start) / 2:5.3f} / km\n"
        out += f"   Z: {(self.z_end + self.z_start) / 2:5.3f} / km\n"
        out += f" Dip: {self.dip_angle:5.3f}\n"
        out += "Components of opening:\n"
        out += f"   Tensile: {self.tensile_opening:5.3f} / m³\n"
        out += f" Inflation: {self.point_opening:5.3f} / m³\n"

        return out

    @property
    def raw_input(self) -> list:
        """Prepare the raw input as a 10 element wide list."""

        raw_input = [
            self.x_start,
            self.y_start,
            self.x_end,
            self.y_end,
            500,
            self.tensile_opening,
            self.point_opening,
            self.dip_angle,
            self.z_start,
            self.z_end,
        ]

        return raw_input

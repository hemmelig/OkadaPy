"""

:copyright:
    2024, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

from dataclasses import dataclass


@dataclass
class PlanarFault:
    """
    Class to encapsulate a planar fault source.

    """

    x_start: float
    x_end: float
    y_start: float
    y_end: float
    z_start: float
    z_end: float
    dip_angle: float
    right_lateral_slip: float = 0.0
    dip_slip: float = 0.0

    def __str__(self):
        """Returns a summary of the planar fault."""

        out = "Planar fault:\n"
        out += f"   X-start: {self.x_start:5.3f} / km; X-end: {self.x_end:5.3f} / km\n"
        out += f"   Y-start: {self.y_start:5.3f} / km; Y-end: {self.y_end:5.3f} / km\n"
        out += f"   Z-start: {self.z_start:5.3f} / km; Z-end: {self.z_end:5.3f} / km\n"
        out += f" Dip: {self.dip_angle:5.3f}\n"
        out += "Components of shear:\n"
        out += f" Right-lateral: {self.right_lateral_slip:5.3f} / m³\n"
        out += f"           Dip: {self.dip_slip:5.3f} / m³\n"

        return out

    @property
    def raw_input(self) -> list:
        """Prepare the raw input as a 10 element wide list."""

        raw_input = [
            self.x_start,
            self.y_start,
            self.x_end,
            self.y_end,
            100,
            self.right_lateral_slip,
            self.dip_slip,
            self.dip_angle,
            self.z_start,
            self.z_end,
        ]

        return raw_input

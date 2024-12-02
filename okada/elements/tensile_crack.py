"""

:copyright:
    2024, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

from dataclasses import dataclass


@dataclass
class TensileCrack:
    """
    Class to encapsulate a tensile crack, capable of tensile opening, strike-slip, and
    dip-slip.

    """

    x_start: float
    x_end: float
    y_start: float
    y_end: float
    z_start: float
    z_end: float
    dip_angle: float
    tensile_slip: float | None = None
    strike_slip: float | None = None
    dip_slip: float | None = None

    def __str__(self):
        """Returns a summary of the tensile crack."""

        out = "Point source of shear:\n"
        out += f"   X-start: {self.x_start:5.3f} / km; X-end: {self.x_end:5.3f} / km\n"
        out += f"   Y-start: {self.y_start:5.3f} / km; Y-end: {self.y_end:5.3f} / km\n"
        out += f"   Z-start: {self.z_start:5.3f} / km; Z-end: {self.z_end:5.3f} / km\n"
        out += f"       Dip: {self.dip_angle:5.3f}\n"
        out += "Components of opening and shear:\n"
        if self.tensile_slip is not None:
            out += f"     Tensile: {self.tensile_slip:5.3f} / m\n"
        if self.strike_slip is not None:
            out += f" Strike-slip: {self.strike_slip:5.3f} / m\n"
        if self.dip_slip is not None:
            out += f"    Dip-slip: {self.dip_slip:5.3f} / m\n"

        return out

    @property
    def raw_input(self) -> list:
        """Prepare the raw input as a 10 element wide list."""

        if self.dip_slip is None:
            strike_slip = 0.0 if self.strike_slip is None else self.strike_slip
            raw_input = [
                self.x_start,
                self.y_start,
                self.x_end,
                self.y_end,
                200,
                strike_slip,
                self.tensile_slip,
                self.dip_angle,
                self.z_start,
                self.z_end,
            ]
        elif self.strike_slip is None:
            dip_slip = 0.0 if self.dip_slip is None else self.dip_slip
            raw_input = [
                self.x_start,
                self.y_start,
                self.x_end,
                self.y_end,
                300,
                self.tensile_slip,
                dip_slip,
                self.dip_angle,
                self.z_start,
                self.z_end,
            ]

        return raw_input

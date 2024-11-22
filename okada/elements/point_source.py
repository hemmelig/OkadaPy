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
    """
    Class to encapsulate point source capable of shearing.

    """

    x_start: float
    x_end: float
    y_start: float
    y_end: float
    z_start: float
    z_end: float
    dip_angle: float
    right_lateral_potency: float = 0.0
    reverse_potency: float = 0.0


@dataclass
class PointInflation:
    """
    Class to encapsulate a point source capable of inflation/deflation.

    """

    x_start: float
    x_end: float
    y_start: float
    y_end: float
    z_start: float
    z_end: float
    dip_angle: float
    tensile_opening: float = 0.0
    point_opening: float = 0.0

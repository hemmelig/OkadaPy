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

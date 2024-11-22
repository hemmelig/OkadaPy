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
    tensile_slip: float = 0.0
    strike_slip: float = 0.0
    dip_slip: float = 0.0

"""

:copyright:
    2025, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

from dataclasses import dataclass

import numpy as np
from numpy import deg2rad as d2r


@dataclass
class ReceiverFault:
    """
    Class to encapsulate a receiver fault.

    Faults are defined following the definition of Aki & Richards.

        - Strike is measured from north
        - Dip is measured DOWN from the horizontal (that is, dip angle must be positive)
        - Rake is the angle between the strike vector and the slip vector, in the slip
          plane

    """

    strike: float
    dip: float
    rake: float

    def __post_init__(self) -> None:
        """Perform data validation."""

        while self.strike >= 180.0:
            self.strike -= 180.0
            self.dip *= -1  # Necessary? Dip should be always positive.

        self.rake -= 90.0
        while self.rake <= -180.0:
            self.rake += 360.0

    @property
    def strike_vector(self) -> np.ndarray:
        """Compose the vector parallel to fault strike."""

        return np.array(
            [
                np.cos(d2r(self.strike)),
                np.sin(d2r(self.strike)),
                0,
            ]
        )

    @property
    def dip_vector(self) -> np.ndarray:
        """Compose the vector parallel to the fault dip."""

        return np.array(
            [
                np.cos(d2r(self.dip)) * np.sin(d2r(self.strike)),
                -np.cos(d2r(self.dip)) * np.cos(d2r(self.strike)),
                -np.sin(d2r(self.dip)),
            ]
        )

    @property
    def normal_vector(self) -> np.ndarray:
        """Compose the vector normal to the fault plane from the strike and dip."""

        return np.cross(self.strike_vector, self.dip_vector)

    @property
    def slip_vector(self) -> np.ndarray:
        """Compute the slip vector from the strike and dip vectors and the rake."""

        along_strike = np.cos(d2r(self.rake)) * self.strike_vector
        along_dip = np.sin(d2r(self.rake)) * self.dip_vector

        return along_strike + along_dip

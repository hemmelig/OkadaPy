"""

:copyright:
    2024, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

from .planar_source import PlanarFault
from .point_source import PointInflation, PointShear
from .tensile_crack import TensileCrack


def coulomb2okadapy(element: list[float]):
    """"""

    kode = element[4]

    match kode:
        case 100:
            okadapy_element = PlanarFault(
                x_start=element[0],
                y_start=element[1],
                x_end=element[2],
                y_end=element[3],
                z_start=element[8],
                z_end=element[9],
                right_lateral_slip=element[5],
                dip_slip=element[6],
                dip_angle=element[7],
            )
        case 200:
            okadapy_element = TensileCrack(
                x_start=element[0],
                y_start=element[1],
                x_end=element[2],
                y_end=element[3],
                z_start=element[8],
                z_end=element[9],
                strike_slip=element[5],
                tensile_slip=element[6],
                dip_angle=element[7],
            )
        case 300:
            okadapy_element = TensileCrack(
                x_start=element[0],
                y_start=element[1],
                x_end=element[2],
                y_end=element[3],
                z_start=element[8],
                z_end=element[9],
                tensile_slip=element[5],
                dip_slip=element[6],
                dip_angle=element[7],
            )
        case 400:
            okadapy_element = PointShear(
                x_start=element[0],
                y_start=element[1],
                x_end=element[2],
                y_end=element[3],
                z_start=element[8],
                z_end=element[9],
                right_lateral_potency=element[5],
                reverse_potency=element[6],
                dip_angle=element[7],
            )
        case 500:
            okadapy_element = PointInflation(
                x_start=element[0],
                y_start=element[1],
                x_end=element[2],
                y_end=element[3],
                z_start=element[8],
                z_end=element[9],
                tensile_opening=element[5],
                point_opening=element[6],
                dip_angle=element[7],
            )

    return okadapy_element

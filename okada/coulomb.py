"""
Compute the Coulomb stress field for faults with a given orientation.

:copyright:
    2025, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import numpy as np

from okada.elements import ReceiverFault


def calculate_coulomb_stress(
    stress_tensor: np.ndarray, fault: ReceiverFault, friction_coefficient: float
) -> np.ndarray:
    """
    Resolve the shear and normal stresses on some arbitrarily oriented receiver fault
    (i.e., no slip occurs) with a given coefficient of friction.

    """

    # Traction vector on fault plane given by T = σ ⋅ n
    traction = np.einsum("ijk,jk->ik", stress_tensor, fault.normal_vector)

    # Normal stress given by σ_n = T ⋅ n
    normal = np.einsum("ik,ik->k", traction, fault.normal_vector)

    # Shear stress in slip direction given by τ = T ⋅ slip_vector
    shear = np.einsum("ik,ik->k", traction, fault.slip_vector)

    coulomb_stress = shear + friction_coefficient * normal

    return coulomb_stress, shear, normal

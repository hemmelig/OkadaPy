"""

:copyright:
    2024, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

from dataclasses import dataclass

import numpy as np


STRESS_COMPONENT_MAP = {
    "XX": 0,
    "YY": 1,
    "ZZ": 2,
    "YZ": 3,
    "ZY": 3,
    "XZ": 4,
    "ZX": 4,
    "XY": 5,
    "YX": 5,
}


@dataclass
class Result:
    """Class to encapsulate computed results for `OkadaPy`."""

    x_coords: np.ndarray
    y_coords: np.ndarray


@dataclass
class DisplacementResult(Result):
    """Class to encapsulate computed displacement results for `OkadaPy`."""

    displacement: np.ndarray


@dataclass
class StrainResult(Result):
    """Class to encapsulate computed strain results for `OkadaPy`."""

    strain: np.ndarray


@dataclass
class StressResult(Result):
    """Class to encapsulate computed stress results for `OkadaPy`."""

    stress: np.ndarray

    def shmax_vectors(self):
        shmax_vectors = []
        for x in self.stress:
            for stress_tensor in x:
                stress_tensor_2d = np.asarray(
                    [
                        stress_tensor[STRESS_COMPONENT_MAP["XX"]],
                        stress_tensor[STRESS_COMPONENT_MAP["XY"]],
                        stress_tensor[STRESS_COMPONENT_MAP["YX"]],
                        stress_tensor[STRESS_COMPONENT_MAP["YY"]],
                    ]
                ).reshape(2, 2)
                eigenvals, eigenvecs = np.linalg.eigh(stress_tensor_2d)
                max_idx = np.argmin(eigenvals)

                max_direction = eigenvecs[:, max_idx]
                max_magnitude = np.min(eigenvals)

                angle = np.arctan2(max_direction[1], max_direction[0])

                shmax_vectors.append((max_direction, max_magnitude, angle))

        return shmax_vectors

        # for stress_tensor in stress_tensors:
        #     stress_tensor = np.array(
        #         [
        #             [stress_tensor[0], stress_tensor[5], stress_tensor[4]],
        #             [stress_tensor[5], stress_tensor[1], stress_tensor[3]],
        #             [stress_tensor[4], stress_tensor[3], stress_tensor[2]],
        #         ]
        #     )

        #     eigenvalues, eigenvectors = eigh(stress_tensor)
        # s1n, s2n = eigenvectors[1][:2]
        # s1e, s2e = eigenvectors[0][:2]

        # R = (eigenvalues[1] - eigenvalues[0]) / (eigenvalues[2] - eigenvalues[0])
        # X = (s1n ** 2 - s1e ** 2) + (1 - R) * (s2n ** 2 - s2e ** 2)
        # Y = 2 * (s1n * s1e + (1 - R) * s2n * s2e)
        # tan2alpha = Y / X
        # alpha = np.arctan(tan2alpha) / 2
        # sl = np.sqrt(50) / 4
        # angle = np.rad2deg(alpha)
        # shx = sl * np.sin(alpha)
        # shy = sl * np.cos(alpha)

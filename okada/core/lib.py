"""
Bindings for the C library functions:

    - _evaluate_displacement
    - _evaluate_strain
    - _evaluate_stress

:copyright:
    2024, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import numpy as np
import numpy.ctypeslib as clib

from okada.core.libnames import _load_cdll
from okada.model import Model
from okada.results import DisplacementResult, Result, StrainResult, StressResult
from okada.utils import timeit


libokada = _load_cdll("libokada")

c_int32 = clib.ctypes.c_int32
c_dbl = clib.ctypes.c_double
c_dPt = clib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")


@timeit()
def evaluate_okada_model(
    model: Model, calculation_depth: float, mode: str, threads: int
) -> Result:
    """Evaluate the model using the analytical equations presented in Okada, 1992."""

    match mode:
        case "displacement":
            print("\tComputing displacement solution...")
            result = _evaluate_displacement(model, calculation_depth, threads)
        case "strain":
            print("\tComputing strain solution...")
            result = _evaluate_strain(model, calculation_depth, threads)
        case "stress":
            print("\tComputing stress solution...")
            result = _evaluate_stress(model, calculation_depth, threads)

    return result


libokada.compute_okada_displacement.argtypes = [
    c_dPt,
    c_dPt,
    c_int32,
    c_dPt,
    c_int32,
    c_dbl,
    c_dbl,
    c_dbl,
    c_dPt,
    c_int32,
]


def _evaluate_displacement(
    model: Model, calculation_depth: float, threads: int
) -> DisplacementResult:
    """Evaluate the model using the analytical equations presented in Okada, 1992."""

    n_coords = len(model.x_coords)
    nx, ny = len(set(model.x_coords)), len(set(model.y_coords))
    n_elements = int(len(model.raw_elements) / 10)

    displacement = np.zeros(n_coords * 12, dtype=np.float64, order="C")
    displacement = np.ascontiguousarray(displacement)

    libokada.compute_okada_displacement(
        model.x_coords,
        model.y_coords,
        c_int32(n_coords),
        model.raw_elements,
        c_int32(n_elements),
        model.youngs_modulus,
        model.poisson_ratio,
        calculation_depth,
        displacement,
        c_int32(threads),
    )

    kwargs = {
        "displacement": displacement.reshape(nx, ny, 12),
        "x_coords": model.x_coords,
        "y_coords": model.y_coords,
    }

    return DisplacementResult(**kwargs)


libokada.compute_okada_strain.argtypes = [
    c_dPt,
    c_dPt,
    c_int32,
    c_dPt,
    c_int32,
    c_dbl,
    c_dbl,
    c_dbl,
    c_dPt,
    c_int32,
]


def _evaluate_strain(
    model: Model,
    calculation_depth: float,
    threads: int,
) -> StrainResult:
    """Evaluate the model using the analytical equations presented in Okada, 1992."""

    n_coords = len(model.x_coords)
    nx, ny = len(set(model.x_coords)), len(set(model.y_coords))
    n_elements = int(len(model.raw_elements) / 10)

    strain = np.zeros(n_coords * 6, dtype=np.float64, order="C")
    strain = np.ascontiguousarray(strain)

    libokada.compute_okada_strain(
        model.x_coords,
        model.y_coords,
        c_int32(n_coords),
        model.raw_elements,
        c_int32(n_elements),
        model.youngs_modulus,
        model.poisson_ratio,
        calculation_depth,
        strain,
        c_int32(threads),
    )

    kwargs = {
        "strain": strain.reshape(nx, ny, 6),
        "x_coords": model.x_coords,
        "y_coords": model.y_coords,
    }

    return StrainResult(**kwargs)


libokada.compute_okada_stress.argtypes = [
    c_dPt,
    c_dPt,
    c_int32,
    c_dPt,
    c_int32,
    c_dbl,
    c_dbl,
    c_dbl,
    c_dPt,
    c_int32,
]


def _evaluate_stress(
    model: Model,
    calculation_depth: float,
    threads: int,
) -> StressResult:
    """Evaluate the model using the analytical equations presented in Okada, 1992."""

    n_coords = len(model.x_coords)
    nx, ny = len(set(model.x_coords)), len(set(model.y_coords))
    n_elements = int(len(model.raw_elements) / 10)

    stress = np.zeros(n_coords * 6, dtype=np.float64, order="C")
    stress = np.ascontiguousarray(stress)

    libokada.compute_okada_stress(
        model.x_coords,
        model.y_coords,
        c_int32(n_coords),
        model.raw_elements,
        c_int32(n_elements),
        model.youngs_modulus,
        model.poisson_ratio,
        calculation_depth,
        stress,
        c_int32(threads),
    )

    kwargs = {
        "stress": stress.reshape(nx, ny, 6),
        "x_coords": model.x_coords,
        "y_coords": model.y_coords,
    }

    return StressResult(**kwargs)

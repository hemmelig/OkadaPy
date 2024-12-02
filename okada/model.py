"""

:copyright:
    2024, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import itertools
from dataclasses import dataclass
import pathlib
import tomllib

import numpy as np
import pandas as pd
import pyproj
from pyproj.enums import TransformDirection

from okada.elements import coulomb2okadapy


@dataclass
class Model:
    """Class to encapsulate model information for `OkadaPy`."""

    poisson_ratio: float
    youngs_modulus: float
    friction_coefficient: float
    elements: list
    x_coords: np.ndarray
    y_coords: np.ndarray
    size: list[float] | None = None
    xsection: list[float] | None = None
    map_: list[float] | None = None
    transformer: pyproj.Transformer | None = None

    @property
    def grid_bounds_xy(self) -> list[tuple[float], tuple[float]]:
        """Get x and y coordinate bounds for the model grid."""

        bottom_left_corner = (min(self.x_coords), min(self.y_coords))
        top_right_corner = (max(self.x_coords), max(self.y_coords))

        return [bottom_left_corner, top_right_corner]

    @property
    def grid_bounds_coords(self) -> list[tuple[float], tuple[float]]:
        """Get lat/lon coordinate bounds for the model grid."""

        (min_lon, min_lat), (max_lon, max_lat) = [
            self.transformer.transform(
                x * 1000, y * 1000, direction=TransformDirection.INVERSE
            )
            for x, y in self.grid_bounds_xy
        ]
        return [min_lon, max_lon, min_lat, max_lat]

    @property
    def grid_xy(self) -> tuple[np.ndarray, np.ndarray]:
        """Get x and y coordinates of all nodes in model grid."""

        X, Y = np.mgrid[
            min(self.x_coords) : max(self.x_coords) : len(set(self.x_coords)) * 1j,
            min(self.y_coords) : max(self.y_coords) : len(set(self.y_coords)) * 1j,
        ]

        return X, Y

    @property
    def grid_coords(self) -> tuple[np.ndarray, np.ndarray]:
        """Get lat/lon coordinates of all nodes in model grid."""

        if self.transformer is None:
            print("No coordinate transformation defined for this model.")
            return

        grid_x, grid_y = self.grid_xy
        grid_lons, grid_lats = self.transformer.transform(
            grid_x * 1000,
            grid_y * 1000,
            direction=TransformDirection.INVERSE,
        )

        return grid_lons, grid_lats

    @property
    def raw_elements(self) -> np.ndarray:
        """
        Prepare all model elements for computation, which assumes parameters have
        specific positions in a flattened, contiguous array.

        """

        return np.asarray([element.raw_input for element in self.elements]).flatten()


def read(model_file: str, file_format: str = None) -> Model:
    """
    Parse a model file into the Model dataclass.

    Arguments
    ---------
    model_file: Path to a file containing the model to load.
    file_format: A string identifier for a specific file format e.g. "COULOMB"

    Returns
    -------
    model: The model in the `OkadaPy` Model format.

    Raises
    ------
    FileNotFoundError: If the specified model_file does not exist.
    Exception: If the model file format could not be detected.

    """

    print(f"\tReading in model from {model_file}...")

    model_file = pathlib.Path(model_file)
    if not model_file.is_file():
        raise FileNotFoundError

    # Read file. Will try to automatically detect the file format
    if model_file.suffix == ".toml" or file_format == "TOML":
        print("\t\t.toml file identified.")
        model = _read_toml(model_file)
    elif model_file.suffix == ".inp" or file_format == "COULOMB":
        print("\t\tCoulomb-style .inp file identified.")
        model = _read_coulomb(model_file)
    else:
        raise Exception("Could not detect file format—please specify.")
    print("\t...read successful.")

    return model


def _read_coulomb(model_file: pathlib.Path) -> Model:
    """
    Parse the model and model parameters from a file that conforms to the Coulomb3
    model file format.

    Arguments
    ---------
    model_file: Path to a Coulomb3 model file containing the model to load.

    Returns
    -------
    model: The model in the `OkadaPy` Model format.

    """

    with open(model_file, "r") as f:
        lines = f.readlines()

    header, body, *tail = [
        list(group)
        for k, group in itertools.groupby(lines, lambda x: x == "\n")
        if not k
    ]
    tail = [section for sublist in tail for section in sublist]

    # Parse information from model file header
    poisson_ratio = float(header[3].strip().split()[1])
    youngs_modulus = float(header[4].strip().split()[1])
    friction_coefficient = float(header[6].strip().split()[1])

    # Parse model elements from model file body
    raw_elements = (
        pd.DataFrame([line.strip().split() for line in body[2:]])
        .drop(columns=0)
        .values.astype(np.float64)
    )
    elements = [coulomb2okadapy(element) for element in raw_elements]

    # Parse grid information from model file tail
    grid, size, xsection, map_ = tail[:7], tail[7:11], tail[11:19], tail[19:]

    grid_keys = ["x_min", "y_min", "x_max", "y_max", "x_inc", "y_inc"]
    grid = {
        key: float(model_part.split("=")[1].strip())
        for key, model_part in zip(grid_keys, grid[1:])
    }

    x_coords = np.arange(grid["x_min"], grid["x_max"] + grid["x_inc"], grid["x_inc"])
    x_coords = np.ascontiguousarray(x_coords, dtype=np.float64)

    y_coords = np.arange(grid["y_min"], grid["y_max"] + grid["y_inc"], grid["y_inc"])
    y_coords = np.ascontiguousarray(y_coords, dtype=np.float64)

    x_coords, y_coords = np.meshgrid(x_coords, y_coords, indexing="ij")

    map_keys = ["lon_min", "lon_max", "lon_zero", "lat_min", "lat_max", "lat_zero"]
    map_ = {
        key: float(map_part.split("=")[1].strip())
        for key, map_part in zip(map_keys, map_[1:])
    }

    model = Model(
        poisson_ratio,
        youngs_modulus,
        friction_coefficient,
        elements,
        x_coords.flatten(),
        y_coords.flatten(),
        size,
        xsection,
        map_,
    )

    return model


def _read_toml(model_file: pathlib.Path) -> Model:
    """
    Parse the model and model parameters from a .toml file (defined by `OkadaPy`).

    Arguments
    ---------
    model_file: Path to an OkadaPy .toml model file containing the model to load.

    Returns
    -------
    model: The model in the `OkadaPy` Model format.

    """

    with model_file.open("rb") as f:
        tf = tomllib.load(f)

    elements = pd.read_csv(model_file.parent / tf["model"]["elements_file"])
    elements = elements.values.flatten()
    elements = elements.astype(np.float64)

    if "grid" in tf.keys():
        x_coords = np.arange(*tf["grid"]["x_range"])
        y_coords = np.arange(*tf["grid"]["y_range"])
        z_coords = np.arange(*tf["grid"]["z_range"])
    elif "coords" in tf.keys():
        x_coords = tf["coords"]["x"]
        y_coords = tf["coords"]["y"]
        z_coords = tf["coords"]["z"]

    x_coords = np.ascontiguousarray(x_coords, dtype=np.float64)
    y_coords = np.ascontiguousarray(y_coords, dtype=np.float64)
    z_coords = np.ascontiguousarray(z_coords, dtype=np.float64)

    x_coords, y_coords = np.meshgrid(x_coords, y_coords, indexing="ij")

    model = Model(
        tf["model"]["poisson_ratio"],
        float(tf["model"]["youngs_modulus"]),
        tf["model"]["friction_coefficient"],
        elements,
        x_coords.flatten(),
        y_coords.flatten(),
    )

    return model

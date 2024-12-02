"""
Plot a given component of the strain tensor field for a given deformation model.

:copyright:
    2024, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pygmt

from okada import Model
from okada.results import StressResult


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


def plot(*args, **kwargs) -> plt.Axes | pygmt.Figure | None:
    from okada.plot import PlottingBackendConfig

    match PlottingBackendConfig.requested:
        case "matplotlib":
            return _plot_mpl(*args, **kwargs)
        case "pygmt":
            return _plot_pygmt(*args, **kwargs)


def _plot_mpl(
    model: Model,
    stress: StressResult,
    ax: plt.Axes | None = None,
    save: bool = False,
    outfile: str = "stress_field.pdf",
    xy_files: list[pathlib.Path] | None = None,
    coordinate_space: str = "cartesian",
) -> plt.Axes | None:
    """
    Plot a map view of a given stress field, with maximum horizontal
    Plot the vertical and horizontal displacement fields for a given deformation model.

    The vertical displacement field is shown as a colour map, while the horizontal
    displacement field is represented by 2-D vectors.

    """

    pass


def _plot_pygmt(
    model: Model,
    displacement: np.ndarray,
    figure: pygmt.Figure | None = None,
    projection: str = "M18c",
    frame_params: list = ["WSne", "xa1f0.1", "ya1f0.1"],
    save: bool = False,
    outfile: str = "displacement_field.pdf",
    xy_files: list[pathlib.Path] | None = None,
    grdsample_params: dict = {"interpolation": "b", "spacing": "0.01k"},
) -> pygmt.Figure | None:
    pass


def plot_stress_vectors(
    model: Model,
    stress: StressResult,
    ax: plt.Axes | None = None,
    save: bool = False,
    coordinate_space: str = "cartesian",
) -> plt.Axes:
    """ """

    X, Y = model.grid_xy

    if coordinate_space == "geographic":
        if (transformer := model.transformer) is None:
            print(
                "No coordinate transformation specified, plotting in Cartesian space."
            )
        else:
            X, Y = model.grid_coords

    if coordinate_space == "cartesian":
        ax.set_aspect("equal")
    elif coordinate_space == "geographic" and transformer is not None:
        x_extent, y_extent = [max_ - min_ for min_, max_ in zip(*model.grid_bounds_xy)]
        min_lon, max_lon, min_lat, max_lat = model.grid_bounds_coords
        lon_extent, lat_extent = max_lon - min_lon, max_lat - min_lat
        aspect = (y_extent * lon_extent) / (x_extent * lat_extent)
        ax.set_aspect(aspect=aspect)

    max_directions, max_magnitudes, angles = zip(*stress.shmax_vectors())

    if ax is None:
        fig, ax = plt.subplots(1, figsize=(10, 5), constrained_layout=True)
    else:
        fig = ax.get_figure()

    normalized_magnitudes = max_magnitudes / np.max(max_magnitudes)

    Q = ax.quiver(
        X,
        Y,
        [d[0] for d in max_directions],
        [d[1] for d in max_directions],
        normalized_magnitudes,
        cmap="viridis",
        scale=10,
    )

    return Q, ax

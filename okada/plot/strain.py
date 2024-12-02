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
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import pygmt

from okada import Model
from okada.results import StrainResult


STRAIN_COMPONENT_MAP = {
    "XX": 0,
    "YY": 1,
    "ZZ": 2,
    "YZ": 3,
    "XZ": 4,
    "XY": 5,
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
    strain_result: StrainResult,
    strain_component: str,
    ax: plt.Axes | None = None,
    save: bool = False,
    outfile: str = "strain_field.pdf",
    xy_files: list[pathlib.Path] | None = None,
    coordinate_space: str = "cartesian",
) -> plt.Axes | None:
    """
    Plot the vertical and horizontal displacement fields for a given deformation model.

    The vertical displacement field is shown as a colour map, while the horizontal
    displacement field is represented by 2-D vectors.

    """

    X, Y = model.grid_xy

    if coordinate_space == "geographic":
        if (transformer := model.transformer) is None:
            print(
                "No coordinate transformation specified, plotting in Cartesian space."
            )
        else:
            X, Y = model.grid_coords

    strain = strain_result.strain
    if strain_component == "DIL":
        strain_to_plot = (strain[:, :, 0] + strain[:, :, 1] + strain[:, :, 2]) * 1e6
    else:
        strain_to_plot = strain[:, :, STRAIN_COMPONENT_MAP[strain_component]] * 1e6

    # Plot vertical displacement as colormap
    tick_limit = max(
        [
            abs(value)
            for value in [
                strain_to_plot[:, :].flatten().min(),
                strain_to_plot[:, :].flatten().max(),
            ]
        ]
    )
    levels = MaxNLocator(nbins=100).tick_values(
        -tick_limit,
        tick_limit,
    )
    cmap = plt.colormaps["RdBu_r"]

    if ax is None:
        fig, ax = plt.subplots(1, figsize=(10, 5), constrained_layout=True)
    else:
        fig = ax.get_figure()

    v_scale = max(abs(strain_to_plot[:, :].flatten()))
    cf = ax.contourf(
        X,
        Y,
        strain_to_plot[:, :],
        levels=levels,
        cmap=cmap,
        vmin=-v_scale,
        vmax=v_scale,
    )
    if coordinate_space == "cartesian":
        ax.set_aspect("equal")
    elif coordinate_space == "geographic" and transformer is not None:
        x_extent, y_extent = [max_ - min_ for min_, max_ in zip(*model.grid_bounds_xy)]
        min_lon, max_lon, min_lat, max_lat = model.grid_bounds_coords
        lon_extent, lat_extent = max_lon - min_lon, max_lat - min_lat
        aspect = (y_extent * lon_extent) / (x_extent * lat_extent)
        ax.set_aspect(aspect=aspect)

    # fig.colorbar(cf, ax=ax, shrink=0.5)

    if xy_files is not None:
        for xy_file in xy_files:
            xy_file = pd.read_csv(
                xy_file,
                names=["Longitude", "Latitude", "Z"],
                header=None,
                comment=">",
                sep="\s+",
            )
            if coordinate_space == "cartesian" and transformer is not None:
                x, y = (
                    transformer.transform(xy_file["Longitude"], xy_file["Latitude"])
                    / 1000
                )
            elif coordinate_space == "geographic":
                x, y = xy_file["Longitude"], xy_file["Latitude"]
            ax.plot(x, y, color="k")

    ax.set_xlim([min(X.flatten()), max(X.flatten())])
    ax.set_ylim([min(Y.flatten()), max(Y.flatten())])

    if save:
        plt.savefig(outfile)

    return cf, ax


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

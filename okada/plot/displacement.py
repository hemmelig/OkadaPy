"""
Plot the displacement field for a given deformation model.

:copyright:
    2024, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import pathlib
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import pygmt
from pyproj.enums import TransformDirection

from okada import Model


def plot(*args, **kwargs) -> plt.Axes | pygmt.Figure | None:
    from okada.plot import PlottingBackendConfig

    match PlottingBackendConfig.requested:
        case "matplotlib":
            return _plot_mpl(*args, **kwargs)
        case "pygmt":
            return _plot_pygmt(*args, **kwargs)


def _plot_mpl(
    model: Model,
    displacement,
    ax: plt.Axes | None = None,
    save: bool = False,
    outfile: str = "displacement_field.pdf",
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

    # Plot vertical displacement as colormap
    levels = MaxNLocator(nbins=100).tick_values(
        displacement[:, :, 2].flatten().min(), displacement[:, :, 2].flatten().max()
    )
    cmap = plt.colormaps["RdBu_r"]

    if ax is None:
        fig, ax = plt.subplots(1, figsize=(10, 5), constrained_layout=True)
    else:
        fig = ax.get_figure()

    v_scale = max(abs(displacement[:, :, 2].flatten()))
    cf = ax.contourf(
        X,
        Y,
        displacement[:, :, 2],
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

    fig.colorbar(cf, ax=ax, shrink=0.5)

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
                    transformer.transform(xy_file["Latitude"], xy_file["Longitude"])
                    / 1000
                )
            elif coordinate_space == "geographic":
                x, y = xy_file["Longitude"], xy_file["Latitude"]
            ax.plot(x, y, color="k")

    # Plot horizontal displacement as vector arrows
    ax.quiver(X, Y, displacement[:, :, 0], displacement[:, :, 1], scale=3)

    ax.set_xlim([min(X.flatten()), max(X.flatten())])
    ax.set_ylim([min(Y.flatten()), max(Y.flatten())])

    if save:
        plt.savefig(outfile)

    return ax


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
    """
    Plot the vertical and horizontal displacement fields for a given deformation model.

    The vertical displacement field is shown as a colour map, while the horizontal
    displacement field is represented by 2-D vectors.

    """

    if figure is None:
        figure = pygmt.Figure()
        pygmt.config(FORMAT_GEO_MAP="ddd.xx", MAP_FRAME_TYPE="plain")

    if model.transformer is None:
        print("Cannot use PyGMT backend without a defined coordinate transformation.")
        sys.exit(1)

    region = model.grid_bounds_coords
    grid_lons, grid_lats = model.grid_coords
    displacement_field_grd = pygmt.xyz2grd(
        x=grid_lons.flatten(),
        y=grid_lats.flatten(),
        z=displacement[:, :, 2].flatten(),
        region=region,
        spacing="1k",
    )
    smoothed_field_grd = pygmt.grdsample(
        displacement_field_grd,
        **grdsample_params,
        region=region,
    )

    # Vertical displacement field
    figure.grdimage(
        smoothed_field_grd,
        interpolation="b",
        projection=projection,
        cmap=True,
        frame=frame_params,
    )

    # Horizontal displacement field
    angles = (
        np.arctan2(displacement[:, :, 1].flatten(), displacement[:, :, 0].flatten())
        * 180
        / np.pi
    )
    amplitudes = [
        np.linalg.norm([x, y])
        for x, y in zip(
            displacement[:, :, 0].flatten(), displacement[:, :, 1].flatten()
        )
    ]

    figure.plot(
        region=region,
        projection=projection,
        x=grid_lons.flatten(),
        y=grid_lats.flatten(),
        style="v1.2c+e+h0.5",
        direction=(angles, amplitudes),
        pen="0.8p",
        fill="black",
    )

    if xy_files is not None:
        for xy_file in xy_files:
            figure.plot(xy_file, pen="1.2p,black", projection=projection)

    if save:
        figure.savefig(outfile)

    return figure

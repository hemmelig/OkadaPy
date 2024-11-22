"""
A collection of plotting utilities for OkadaPy.

:copyright:
    2024, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

from typing import Any, Literal

from .displacement import plot as displacement
from .strain import plot as strain
from .stress import plot_stress_vectors as shmax_vectors


class PlottingBackendConfig:
    """*internal* config dict"""

    requested = "matplotlib"  # default is matplotlib
    available = ("matplotlib", "pygmt")
    allow_fallback = False


def set_backend(
    backend: Literal["matplotlib", "pygmt"], force: bool = True, **kwargs: Any
) -> None:
    """
    Select the backend to be used when plotting visualisations on a map.

    Parameters
    ----------
    backend: Specify the backend to be used.
    force: Raises an ImportError when ``okada.plot.get_backend()`` is called and the
        requested backend can not be imported. force=True should be used in user code
        to ensure that the correct backend is being loaded.

    pyseabreeze only
    ----------------
    pyusb_backend: str
        either libusb1, libusb0 or openusb

    """

    if backend not in PlottingBackendConfig.available:
        raise ValueError(
            f"backend not in: {', '.join(PlottingBackendConfig.available)}"
        )

    PlottingBackendConfig.requested = backend
    PlottingBackendConfig.allow_fallback = not force


__all__ = [displacement, set_backend, strain, shmax_vectors]

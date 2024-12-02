"""
Utility module for loading compiled C libraries.

:copyright:
    2024, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import sysconfig
from ctypes import CDLL
from pathlib import Path


def _load_cdll(name: str) -> CDLL:
    """
    Helper function to load a compiled C library.

    Parameters
    ----------
    name: Name of library to load.

    Returns
    -------
    cdll: The C dynamically linked library, i.e. a shared library object.

    """

    lib = (Path(__file__).parent / "src" / name).with_suffix(
        sysconfig.get_config_var("EXT_SUFFIX")
    )
    try:
        cdll = CDLL(str(lib))
    except Exception as e:
        raise ImportError(
            f"Could not load extension library '{name}'.\n\n{e}\n\n"
            "If you have chosen to build and install from source, ensure the C library"
            "has been compiled correctly. See the documentation for more details."
        )

    return cdll

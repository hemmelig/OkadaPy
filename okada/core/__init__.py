"""
The :mod:`okada.core` module provides Python bindings for the library of
compiled C routines implementing the equations set out in Okada, 1992.

:copyright:
    2024, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

from .lib import evaluate_okada_model  # NOQA

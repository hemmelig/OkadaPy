"""
OkadaPy—a Python toolkit for evaluating analytical deformation models.

:copyright:
    2024, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import pkg_resources

from .model import Model, read as read_model
from .core import evaluate_okada_model as evaluate


__all__ = [evaluate, read_model, Model]
__version__ = pkg_resources.get_distribution("okada").version

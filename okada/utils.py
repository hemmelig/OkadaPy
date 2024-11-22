"""
Module of utilities for the `OkadaPy` package.

:copyright:
    2024, Conor A. Bacon.
:license:
    GNU General Public License, Version 3
    (https://www.gnu.org/licenses/gpl-3.0.html)

"""

import time
from functools import wraps


def timeit(*args_, **kwargs_):
    """Function wrapper that measures the time elapsed during its execution."""

    def inner_function(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            ts = time.time()
            result = func(*args, **kwargs)
            print(f"\t\tElapsed time: {time.time() - ts:6f} seconds.\n\t...complete.")
            return result

        return wrapper

    return inner_function

# -*- coding: utf-8 -*-

"""Common non-package specific utility functions.

Licensed under a 3-clause BSD style license - see LICENSE.rst

"""

__author__ = "Nathaniel Starkman"


__all__ = [
    "amuseify_array",
    "draw_unit_normal",
    # "store_function_input",
]


##############################################################################
# IMPORTS

# GENERAL

import numpy as np

# amuse
from amuse.units.core import unit as Unit

# typing
from typing import Optional, Sequence


# CUSTOM

from utilipy.decorators import store_function_input


##############################################################################
# TYPES

_GOOD_DECORATORS: bool


##############################################################################
# CODE
##############################################################################


def amuseify_array(
    arr: np.ndarray, to_unit: Optional[Unit] = None
) -> Sequence:
    """Convert an array of amuse quantities to an amuse quantity array.

    All the units in the array must be of compatible units.

    Paramaters
    ----------
    arr: ndarray
    to_unit: amuse unit

    Returns
    -------
    amuse_arr: amuse narray quantity

    Examples
    --------
    >>> amuseify_array([1 | amu.Myr, 0.002 | amu.Gyr], to_unit=amu.Myr)
    quantity<[1.0, 2.0] Myr>

    """
    if to_unit is None:
        to_unit = arr[0].unit
    amuse_arr = np.array([x.value_in(to_unit) for x in arr]) | to_unit

    return amuse_arr


# /def


# ------------------------------------------------------------------------


def draw_unit_normal(
    ndim, loc: float = 0.0, scale: float = 1.0, size: int = 1, random=0
) -> Sequence:
    """Draw ND-Normal Distribution.

    loc & scale, size apply equally to all dimensions
    use scipy.stats.multivariate_normal for greater customization

    Parameters
    ----------
    ndim: int
        number of dimensions
    loc: float
    scale: float
    size: int
    random: Seed

    Returns
    -------
    draws: (ndim, size) array
        multivariate normal

    """
    # if not isinstance(random, Generator):
    #     random = np.random.default_rng(random)
    if isinstance(random, int):
        random = np.random.RandomState(random)

    # draw random
    draws = np.array(
        [random.normal(loc=loc, scale=scale, size=size) for _ in range(ndim)]
    )
    draws /= np.linalg.norm(draws, axis=0)  # normalizing

    # handle when x=y=z=0
    while (draws == 0).sum() > 0:
        where_0 = draws == 0
        num_where_0 = where_0.sum()

        if num_where_0 > 0:
            draws[where_0] = np.array(
                [
                    random.normal(loc=loc, scale=scale, size=num_where_0)
                    for _ in range(ndim)
                ]
            )
            draws[where_0] /= np.linalg.norm(draws[where_0], axis=0)

    return draws.T


# /def


##############################################################################
# END

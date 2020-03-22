# -*- coding: utf-8 -*-

"""Common non-package specific utility functions.

Licensed under a 3-clause BSD style license - see LICENSE.rst

"""

__author__ = "Nathaniel Starkman"
__license__ = "BSD-3"

# __all__ = [
#     ""
# ]

##############################################################################
# IMPORTS

# GENERAL
import inspect
import functools
import numpy as np

# amuse
from amuse.units.core import unit as Unit

# typing
from typing import Optional


##############################################################################
# CODE
##############################################################################


def amuseify_array(arr: np.ndarray, to_unit: Optional[Unit] = None):
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
    >>> amuseify_array([1 | u.Myr, 0.002 | u.Gyr], to_unit=u.Myr)
    [1, 2] | u.Myr

    """
    if to_unit is None:
        to_unit = arr[0].unit
    amuse_arr = np.array([x.value_in(to_unit) for x in arr]) | to_unit

    return amuse_arr


# /def


def draw_unit_normal(
    ndim, loc: float = 0.0, scale: float = 1.0, size: int = 1, random=0
):
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


def store_function_input(
    function=None,
    *,
    store_inputs: bool = True,
    # _doc_style="numpy", _doc_fmt={}
):
    """Docstring for decorator.

    Description of this decorator

    Parameters
    ----------
    function : types.FunctionType or None, optional
        the function to be decoratored
        if None, then returns decorator to apply.
    _get_inputs : bool, optional
        whether to return all the inputs to the function as a dictionary

    Returns
    -------
    wrapper : types.FunctionType
        wrapper for function
        does a few things
        includes the original function in a method `.__wrapped__`

    Other Parameters
    ----------------
    _doc_style: str or formatter, optional
        default 'numpy'
        parameter to `astroPHD.wraps`
    _doc_fmt: dict, optional
        default None
        parameter to `astroPHD.wraps`

    """
    if function is None:  # allowing for optional arguments
        return functools.partial(
            store_function_input,
            store_inputs=store_inputs,
            # _doc_style=_doc_style,
            # _doc_fmt=_doc_fmt,
        )

    sig = inspect.signature(function)
    # _doc_fmt['wrapped_function'] = function.__name__

    @functools.wraps(
        function,
        # _doc_style=_doc_style, _doc_fmt=_doc_fmt
    )
    def wrapper(*args, store_inputs: bool = store_inputs, **kw):
        """Wrapper docstring.

        Parameters
        ----------
        store_inputs: bool
            whether to store function inputs in a BoundArguments instance
            default {store_inputs}

        Returns
        -------
        inputs: BoundArguments
            the inputs to ``{wrapped_function}``
            only returned if `store_inputs` is True
            other returned values are in now in a tuple

        """
        return_ = function(*args, **kw)
        if store_inputs:
            inputs = sig.bind_partial(*args, **kw)
            inputs.apply_defaults()
            return return_, inputs
        else:
            return return_

    # /def

    return wrapper


# /def


##############################################################################
# END

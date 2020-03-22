#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Amuse-Astropy unit stuff.

Not all of the Amuse & Astropy quantities match exactly

"""

__author__ = "Nathaniel Starkman"

__all__ = ["to_astropy", "to_amuse"]


##############################################################################
# IMPORTS

# GENERAL

import functools

# amuse
from amuse.units.quantities import Quantity


# PROJECT-SPECIFIC

from . import astropy_units as apu
from . import amuse_units as amu


###############################################################################
# CODE
###############################################################################


def to_astropy(quantity):
    """Convert AMUSE quantity to astropy quantity.

    only amuse quantities need to be converted.
    floats and astropy quantities are left as is.

    Parameters
    ----------
    quantity: AMUSE quantity or array-like

    Returns
    -------
    quantity: astropy quantity or array-like

    Notes
    -----
    requires that amuse units are string represented in astropy
    this is why astroPHD is used.

    """
    if isinstance(quantity, Quantity):
        return quantity.value_in(quantity.unit) * apu.Unit(str(quantity.unit))

    else:
        return quantity


# /def


def to_astropy_decorator(function=None, *, arguments=[]):
    """Function decorator to convert inputs to Astropy quantities.

    Parameters
    ----------
    function : types.FunctionType or None, optional
        the function to be decoratored
        if None, then returns decorator to apply.
    arguments : list, optional
        arguments to convert
        integers are indices into `arguments`
        strings are names of `kw` arguments

    Returns
    -------
    wrapper : types.FunctionType
        wrapper for function
        does a few things
        includes the original function in a method `.__wrapped__`

    """
    if function is None:  # allowing for optional arguments
        return functools.partial(to_astropy_decorator, arguments=arguments)

    @functools.wraps(function)
    def wrapper(*args, **kw):
        """Wrapper docstring."""
        for itm in arguments:
            if isinstance(itm, int):
                args[itm] = to_astropy(args[itm])
            elif isinstance(itm, str):
                kw[itm] = to_astropy(kw[itm])
            else:
                raise TypeError("elements of `args` must be int or str")

        return function(*args, **kw)

    # /def

    return wrapper


# /def


# ----------------------------------------------------------------------------


def to_amuse(quantity):
    """Convert AMUSE quantity to astropy quantity.

    only amuse quantities need to be converted.
    floats and astropy quantities are left as is.

    """
    if isinstance(quantity, apu.Quantity):
        return quantity.to_value(quantity.unit) | getattr(
            amu, (str(quantity.unit))
        )

    else:
        return quantity


# /def


def to_amuse_decorator(function=None, *, arguments=[]):
    """Function decorator to convert inputs to AMUSE quantities.

    Parameters
    ----------
    function : types.FunctionType or None, optional
        the function to be decoratored
        if None, then returns decorator to apply.
    arguments : list, optional
        arguments to convert
        integers are indices into `arguments`
        strings are names of `kw` arguments

    Returns
    -------
    wrapper : types.FunctionType
        wrapper for function
        does a few things
        includes the original function in a method `.__wrapped__`

    """
    if function is None:  # allowing for optional arguments
        return functools.partial(to_amuse_decorator, arguments=arguments)

    @functools.wraps(function)
    def wrapper(*args, **kw):
        """Wrapper docstring."""
        for itm in arguments:
            if isinstance(itm, int):
                args[itm] = to_amuse(args[itm])
            elif isinstance(itm, str):
                kw[itm] = to_amuse(kw[itm])
            else:
                raise TypeError("elements of `args` must be int or str")

        return function(*args, **kw)

    # /def

    return wrapper


# /def


# ----------------------------------------------------------------------------


def convert_units_decorator(
    function=None, *, to_astropy_args=[], to_amuse_args=[]
):
    """Function decorator to convert inputs to AMUSE quantities.

    Parameters
    ----------
    function : types.FunctionType or None, optional
        the function to be decoratored
        if None, then returns decorator to apply.
    to_astropy_args : list, optional
        arguments to convert to astropy units
        integers are indices into `args`
        strings are names of `kw` arguments
    to_amuse_args : list, optional
        arguments to convert to amuse units
        integers are indices into `args`
        strings are names of `kw` arguments

    Returns
    -------
    wrapper : types.FunctionType
        wrapper for function
        does a few things
        includes the original function in a method `.__wrapped__`

    """
    if function is None:  # allowing for optional arguments
        return functools.partial(
            convert_units_decorator,
            to_astropy_args=to_astropy_args,
            to_amuse_args=to_amuse_args,
        )

    @functools.wraps(function)
    def wrapper(*args, **kw):
        """Wrapper docstring."""
        for itm in to_astropy_args:
            if isinstance(itm, int):
                args[itm] = to_astropy(args[itm])
            elif isinstance(itm, str):
                kw[itm] = to_astropy(kw[itm])
            else:
                raise TypeError(
                    "elements of `to_astropy_args` must be int or str"
                )
        for itm in to_amuse_args:
            if isinstance(itm, int):
                args[itm] = to_amuse(args[itm])
            elif isinstance(itm, str):
                kw[itm] = to_amuse(kw[itm])
            else:
                raise TypeError(
                    "elements of `to_amuse_args` must be int or str"
                )

        return function(*args, **kw)

    # /def

    return wrapper


# /def

##############################################################################
# END

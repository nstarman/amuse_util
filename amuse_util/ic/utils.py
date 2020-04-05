# -*- coding: utf-8 -*-

"""Utility Functions for Initial Conditions.

Routine Listings
----------------
`num_particles_from_mtot_given_mass_func`
`imf_number_of_particles_decorator`

"""

__author__ = "Nathaniel Starkman"


__all__ = [
    "num_particles_from_mtot_given_mass_func",
    "imf_number_of_particles_decorator",
]


##############################################################################
# IMPORTS

# GENERAL

from typing import Optional, Union, Callable, Tuple, Any

import numpy as np

from amuse.units import units as u


# CUSTOM

from utilipy.utils.functools import wraps, partial


##############################################################################
# PARAMETERS


##############################################################################
# CODE
##############################################################################


def num_particles_from_mtot_given_mass_func(
    target_mass: u.MSun,
    imf_func: Callable,
    imf_args: list = [],
    imf_kwargs: dict = {},
    tolerance: float = 0.01,
    random=0,
) -> Tuple[int, Any, Any]:  # TODO correct Any
    """Finds the number of particles in a cluster of mass `target_mass`.

    Parameters
    ----------
    target_mass: amuse units.MSun
        target mass of cluster
    imf_func: function
        the function that generates the IMF
        signature: ``imf_func(number_of_particles, *args, **kwargs)``
    imf_args: list
        list of arguments for `imf_func`
        should NOT include number_of_particles
    imf_kwargs: dict
        keyword arguments for `imf_func`
        should NOT include `random`
    tolerance: float
        fractional error between `target_mass` and mass(imf(N))
    random: number
        the random seed
        for reproducibility

    Returns
    -------
    N: int
        number of stars in cluster
    M: amuse units.MSun
        mass of stars in cluster
    err: float
        fractional error between M and `target_mass`

    """
    assert target_mass.value_in(u.MSun) >= 0

    def _frac_err(M):
        return np.abs(target_mass - M) / target_mass

    # figuring out good guess
    # by getting the approximate answer from the total / mean mass
    mean = imf_func(1000, *imf_args, random=random, **imf_kwargs).mean()
    N = int(target_mass / mean)
    M = imf_func(N, *imf_args, **imf_kwargs).sum()

    # lower bound, just N / 2
    N_low = int(N / 2)
    M_low = imf_func(N_low, *imf_args, random=random, **imf_kwargs).sum()
    assert M_low < target_mass

    # upper bound, just 2 N
    N_up = int(2 * N)
    M_up = imf_func(N_up, *imf_args, random=random, **imf_kwargs).sum()
    assert M_up > target_mass

    # iterating to find correct N
    # stop when reach desired tolerance
    i: int = 0  # prevent infinite loop
    while (_frac_err(M) > tolerance) & (i < 50):
        i += 1

        if M < target_mass:  # set lower bounds = current
            N_low = N
            M_low = M

        if M > target_mass:  # set upper bounds = current
            N_up = N
            M_up = M

        # set guess
        N = int((N_low + N_up) / 2)
        M = imf_func(N, *imf_args, random=random, **imf_kwargs).sum()

        # catching edge cases
        if N < N_low:
            N = N_low
            M = M_low
            break
        elif N == N_low:
            M_Np1 = imf_func(
                N_low + 1, *imf_args, random=random, **imf_kwargs
            ).sum()
            if _frac_err(M_Np1) < _frac_err(M_low):
                N = N_low + 1
                M = M_Np1
            else:
                N = N_low
                M = M_low
            break
        elif N == N_up:
            M_Nm1 = imf_func(
                N_up - 1, *imf_args, random=random, **imf_kwargs
            ).sum()
            if _frac_err(M_Nm1) < _frac_err(M_up):
                N = N_up - 1
                M = M_Nm1
            else:
                N = N_up
                M = M_up
            break
        elif N > N_up:
            N = N_up
            M = M_up
            break

    return N, M, _frac_err(M)


# /def


# ------------------------------------------------------------------------


def imf_number_of_particles_decorator(
    function: Optional[Callable] = None, *, tolerance: float = 1e-5
) -> Callable:
    """For IMF functions, allows mass argument for `number_of_particles`.

    `number_of_particles` must be a float. If a mass-unit argument is given,
    it is treated as the total mass and the IMF function is used to
    infer the float-value for `number_of_particles`.

    Parameters
    ----------
    function : Callable or None, optional
        the function to be decorated
        if None, then returns decorator to apply.
    tolerance : float, optional
        mass to number_of_particles conversion tolerance
        key-word only argument
        sets the wrapper's default value.

    Returns
    -------
    wrapper : Callable
        wrapper for function
        does a few things
        includes the original function in a method `.__wrapped__`

    """
    if function is None:  # allowing for optional arguments
        return partial(imf_number_of_particles_decorator, tolerance=tolerance)

    # make decorator, specifying docstring style
    @wraps(function, _doc_style="numpy")
    def wrapper(
        number_of_particles: Union[int, Any],  # TODO correct Any
        *func_args: Any,
        m_to_n_tol: float = tolerance,
        **func_kwargs: Any,
    ) -> Any:
        """Wrapper docstring.

        Parameters
        ----------
        m_to_n_tol: float
            fractional error between `target_mass` and mass(imf(N))

        Notes
        -----
        if number_of_particles has units of mass, treated as Mtot, and will
        find number_of_particles via Newton-Raphson iteration.

        """
        try:
            number_of_particles.unit
        except AttributeError:  # it's actually the number of stars
            if not isinstance(number_of_particles, int):
                raise TypeError("number_of_particles (Mtot) is not an integer")
            random = func_kwargs.pop("random", 0)
            pass
        else:  # need to figure out number of stars
            random = func_kwargs.pop("random", 0)
            N, _, _ = num_particles_from_mtot_given_mass_func(
                number_of_particles,
                function,
                tolerance=m_to_n_tol,
                random=random,
                imf_kwargs=func_kwargs,
            )
            number_of_particles = N
            # print(number_of_particles)

        return function(
            number_of_particles, *func_args, random=random, **func_kwargs
        )

    # /def

    return wrapper


# /def


##############################################################################
# END

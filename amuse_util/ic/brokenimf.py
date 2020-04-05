# -*- coding: utf-8 -*-

"""Broken IMF codes.

Routine Listings
----------------
new_kroupa_mass_distribution

"""

__author__ = "Nathaniel Starkman"


__all__ = [
    "new_kroupa_mass_distribution",
]


##############################################################################
# IMPORTS

# GENERAL

import numpy as np

# amuse
from amuse.ic.brokenimf import new_broken_power_law_mass_distribution
from amuse.units import units


# PROJECT-SPECIFIC

from .utils import imf_number_of_particles_decorator


##############################################################################
# CODE
##############################################################################


@imf_number_of_particles_decorator(tolerance=1e-7)
def new_kroupa_mass_distribution(
    number_of_particles: int,
    mass_min: units.MSun = 0.01 | units.MSun,
    mass_max: units.MSun = 100.0 | units.MSun,
    random: bool = True,
):
    """Kroupa (2001) mass distribution in SI units with custom minimum mass.

    Modified version of amuse.ic.brokenimf.new_kroupa_mass_distribution
    that allows for a minimum mass, as well as maximum mass
    The default mass range is [0.01, 0.08, 0.5, 100.0] MSun
    and power-law exponents of each mass range: [-0.3, -1.3, -2.3]

    Parameters
    ----------
    number_of_particles: int
        the number of particles
        {number_of_particles}
    mass_min: float quantity
        the minimum mass, will modify the default minimum mass range
    mass_max: float quantity
        the cut-off mass (defaults to 100.0 MSun)
    random: int, random number generator
        ex: np.random.default_rng(seed=0)
        will default to random seed

    Returns
    -------
    kroupa: (N, ) list quantity
        list of masses
        length `number_of_particles`

    """
    mass_boundaries = [0.01, 0.08, 0.5, 100.0] | units.MSun
    alphas = [-0.3, -1.3, -2.3]

    if mass_min is not None:
        ind = (
            min(np.where(mass_boundaries > mass_min)[0]) - 1
        )  # find where to adj
        mass_boundaries = mass_boundaries[ind:]  # truncate boundaries
        alphas = alphas[ind:]
        if mass_boundaries[0] < mass_min:  # adjust value
            mass_boundaries[0] = mass_min

    kroupa = new_broken_power_law_mass_distribution(
        number_of_particles,
        mass_boundaries=mass_boundaries,
        mass_max=mass_max,
        alphas=alphas,
        random=random,
    )

    return kroupa


# /def


##############################################################################
# END

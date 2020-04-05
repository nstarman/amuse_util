# -*- coding: utf-8 -*-
# see LICENSE.rst

"""Amuse Utilities with improved Astropy compatibility.

See Also
--------
The actual AMUSE docs (https://amuse.readthedocs.io/en/latest/index.html)

References
----------
The AMUSE papers [1]_, [2]_, [3]_, [4]_, [5]_.

.. [1] Portegies Zwart, S. & McMillan, S.L.W., 2018, “Astrophysical Recipes:
    the art of AMUSE”, AAS IOP Astronomy publishing (411 pages)
.. [2] Portegies Zwart, S. et al., 2013, Multi-physics Simulations Using a
    Hierarchical Interchangeable Software Interface,
    Computer Physics Communications 183, 456-468 [2013CoPhC.183..456P]
.. [3] Pelupessy, F. I. et al., 2013, The Astrophysical Multipurpose Software
    Environment, Astronomy and Astrophysics 557, 84 [2013A&A…557A..84P]
.. [4] Portegies Zwart, S. et al., 2009, A multiphysics and multiscale
    software environment for modeling astrophysical systems,
    New Astronomy, Volume 14, Issue 4, 369-378 [2009NewA…14..369P]
.. [5] Portegies Zwart, Simon, van Elteren, Arjen, Pelupessy, Inti, McMillan,
    Steve, Rieder, Steven, de Vries, Nathan, … Altay, Gabriel. (2019, June 28).
    AMUSE: the Astrophysical Multipurpose Software Environment
    (Version v12.0.0). Zenodo. http://doi.org/10.5281/zenodo.3260650

"""

__author__ = "Nathaniel Starkman"
__license__ = "BSD-3"


__all__ = []


##############################################################################
# IMPORTS

# Packages may add whatever they like to this file, but
# should keep this content at the top.
from ._astropy_init import *  # noqa

# GENERAL

from amuse.units import constants, nbody_system

# PROJECT-SPECIFIC

from .datamodel import System, Systems
from .ic import (
    initialize_system,
    recreate_system,
)
from .units import (
    # amuse_units as u,
    amuse_units as amu,
    astropy_units as apu,
    to_astropy,
    to_amuse,
)

# from .utils import amuseify_array

# Import top-level modules for __all__
from . import data, utils, datamodel, ic, simulation, units


##############################################################################
# __ALL__

# Then you can be explicit to control what ends up in the namespace,
# __all__ += ['do_primes']   # noqa
# or you can keep everything from the subpackage with the following instead
# __all__ += example_mod.__all__

__all__ += [
    # GENERAL
    "constants",
    "nbody_system",
    # PROJECT-SPECIFIC
    # datamodel
    "datamodel",
    "System",
    "Systems",
    # ic
    "initialize_system",
    "recreate_system",
    # units
    "amu",
    "apu",
    "to_astropy",
    "to_amuse",
]


##############################################################################
# END

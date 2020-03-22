# -*- coding: utf-8 -*-

# Docstring and Metadata
"""utility file for ic."""

__author__ = "Nathaniel Starkman"

##############################################################################
# IMPORTS

# GENERAL

# PROJECT-SPECIFIC
from ._initialize_system import initialize_system, recreate_system
from .brokenimf import new_kroupa_mass_distribution

from .util import (
    num_particles_from_mtot_given_mass_func,
    imf_number_of_particles_decorator,
)


##############################################################################
# END

# -*- coding: utf-8 -*-

# Docstring
"""Astropy Units.

Merge existing astropy units with AMUSE-compatible additions.

"""

__author__ = "Nathaniel Starkman"
__credits__ = ["astropy"]

__all__ = []


###############################################################################
# IMPORTS

# GENERAL

from astropy import units
from astropy.units import *


# PROJECT-SPECIFIC

from . import _add_to_astropy_units
from ._add_to_astropy_units import *


###############################################################################
# __ALL__

__all__ += units.__all__
__all__ += _add_to_astropy_units.__all__


###############################################################################
# END

# -*- coding: utf-8 -*-

"""Amuse-Astropy unit stuff.

Not all of the Amuse & Astropy quantities match exactly

"""

__author__ = "Nathaniel Starkman"

__all__ = []


##############################################################################
# IMPORTS

# GENERAL
import amuse.units.units as _amu
import amuse.units.core as _core

from amuse.units.units import *  # so this module is a drop-in for amuse units


###############################################################################
# UNITS

mas = _core.named_unit("milliarcsecond", "mas", _amu.arcsec / 1000)


###############################################################################
# __ALL__

# __all__ += _amu.__all__
__all__ += ["mas"]

##############################################################################
# END

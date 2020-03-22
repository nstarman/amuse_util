# -*- coding: utf-8 -*-

"""Amuse Utilities with improved Astropy compatibility.

Licensed under a 3-clause BSD style license - see LICENSE.rst

"""

__author__ = "Nathaniel Starkman"
__license__ = "BSD-3"


##############################################################################
# IMPORTS

# PROJECT-SPECIFIC

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *   # noqa
# ----------------------------------------------------------------------------

__all__ = []
from .example_mod import *   # noqa
# Then you can be explicit to control what ends up in the namespace,
__all__ += ['do_primes']   # noqa
# or you can keep everything from the subpackage with the following instead
# __all__ += example_mod.__all__


##############################################################################
# END

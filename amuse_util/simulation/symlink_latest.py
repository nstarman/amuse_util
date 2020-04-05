# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   :
# AUTHOR  :
# PROJECT :
#
# ----------------------------------------------------------------------------

# Docstring
"""**DOCSTRING**.

TODO, make this a general function

Routing Listings
----------------

"""

__author__ = "Nathaniel Starkman"
# __copyright__ = "Copyright 2019, "
# __credits__ = [""]
# __license__ = "MIT"
# __version__ = "0.0.0"
# __maintainer__ = ""
# __email__ = ""
# __status__ = "Production"

# __all__ = [
#     ""
# ]


###############################################################################
# IMPORTS

# GENERAL

import string
import numpy as np
import os
import os.path

# CUSTOM

# PROJECT-SPECIFIC

###############################################################################
# CODE
###############################################################################


def sorted_date_folders(contents):
    """Return sorted folders that start with a digit."""
    # sort out files
    isfile = np.array([os.path.isfile(c) for c in contents])
    folders = contents[~isfile]

    # filter to date-start files
    folders = [f for f in folders if f[0] in string.digits]

    return np.sort(folders)


# /def


# --------------------------------------------------------------------------


def make_symlink(drct):
    """Make a symlink called `latest` to folder with most recent date name.

    Parameters
    ----------
    drct: str
        the directory in which to make the symbolic link

    """
    old_dir = os.getcwd()
    os.chdir(drct)

    contents = np.array(os.listdir("./"))
    folders = sorted_date_folders(contents)

    if len(folders) > 0:  # not empty
        try:
            os.rmdir("latest")
        except OSError:
            pass
        try:
            os.unlink("latest")
        except OSError:
            pass

        os.symlink("./" + folders[-1], "./latest")

    os.chdir(old_dir)


###############################################################################
# END

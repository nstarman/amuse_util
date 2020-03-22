# -*- coding: utf-8 -*-

# Docstring and Metadata
"""System(s)."""

__author__ = "Nathaniel Starkman"

###############################################################################
# IMPORTS

# GENERAL
import itertools
from dataclasses import dataclass as _dataclass

# amuse
from amuse.units.quantities import ScalarQuantity
from amuse.couple import bridge
from amuse.datamodel.particles import Particles

# typing
from typing import Any

# PROJECT-SPECIFIC
from ._container import AmuseContainer


###############################################################################
# CODE
###############################################################################


def _system_reconstructor(args, kwargs):
    """A reconstructor function for pickling a System.

    Parameters
    ----------
    args: list or tuple
    kwargs: dict

    Returns
    -------
    system

    """
    return System(*args, **kwargs)


# /def


@_dataclass
class System:
    """Class for grouping AMUSE code about a single system.

    Parameters
    ----------
    objects: Particles
        the basic component.
    evolution: evolution code
    gravity: gravity code
    channels

    Returns
    -------
    system

    """

    particles: Particles = None
    evolution: Any = None
    gravity: Any = None

    # util
    converter: Any = None
    channel_attrs: Any = None

    def __getitem__(self, key):
        """Getitem via getattr."""
        return getattr(self, key)

    # /def

    def __setitem__(self, key, value):
        """Setitem via setattr."""
        setattr(self, key, value)

    # /def

    def __post_init__(self):
        """Dataclass post initialization.

        ensures the particles, evolution, and gravity codes are wrapped
        in AmuseContainer

        makes channels between the all the codes

        """
        names = ["particles", "evolution", "gravity"]

        for n in names:
            if (self[n] is not None) and not isinstance(
                self[n], AmuseContainer
            ):
                self[n] = AmuseContainer(self[n], n)

        for n1, n2 in itertools.permutations(names, 2):
            # print("Making channels:")
            if (self[n1] is not None) and (self[n2] is not None):
                print("\t", n1, "<->", n2)
                self[n1].add_channel_to(
                    self[n2], attributes=self.channel_attrs
                )
                self[n1].add_channel_from(
                    self[n2], attributes=self.channel_attrs
                )

    # /def

    def __reduce_ex__(self, protocol):
        """Reduce for serialization.

        TODO FIX
        """
        redx = (
            _system_reconstructor,
            (
                [],
                {
                    "particles": self.particles,
                    "evolution": self.evolution,
                    "gravity": self.gravity,
                    "converter": self.converter,
                    "channel_attrs": self.channel_attrs,
                    # TODO any user-defined attribute
                },
            ),
            None,
            None,
            None,
        )
        return redx

    # /def

    def plot_particles(self, d1="x", d2="y"):
        """In progress."""
        pass

    # /def


# /class


###############################################################################


class Systems:
    """Group of System instances.

    Parameters
    ----------
    use_threading: bool (default False)
        whether to use threading
    **systems: {name: System} pairs

    Returns
    -------
    systems: Systems instance

    """

    def __init__(self, *, use_threading=False, **systems):
        """Initialize a Systems instance."""
        super().__init__()
        self.gravity = bridge.Bridge(use_threading=use_threading)
        self.system_list = list(systems.keys())

        for k, v in systems.items():
            self[k] = v

    # /def

    def __getitem__(self, key):
        """Getitem via getattr."""
        return getattr(self, key)

    # /def

    def __setitem__(self, key, value):
        """Setitem via setattr."""
        setattr(self, key, value)

    # /def

    # def bridge_systems(
    #     self, *systems, timestep=None, dependencies: dict = {0: 1, 1: 0}
    # ):
    #     """Bridge gravity systems.

    #     works on gravity only

    #     Parameters
    #     ----------
    #     system0: str, System, AMUSE thing
    #         the first system to bridge together
    #         if string, use self.`system0`.gravity
    #         if System, use system0.gravity
    #         else use system0

    #     system1: str, System, AMUSE thing
    #         the second system to bridge together
    #         if string, use self.`system1`.gravity
    #         if System, use system1.gravity
    #         else use system1

    #     timestep: None, amuse ScalarQuantity, str, System, AMUSE thing
    #         if None, assume already have set a time step
    #         if number, use
    #         if string, use self.`system1`.gravity
    #         if System, use system1.gravity
    #         else use system1

    #     dependencies: d.get(i, '')i"t
    #         gravitational dependencies of the systems
    #         integer keys, referring to the other system indices
    #         ex :
    #         {0: 1, 1: 0} means sys 0 deps on sys 1 and vice versa

    #         {0: 1, 1: None} means sys 0 deps on sys 1
    #                         but sys 1 does not dep on sys 0

    #     Returns
    #     -------
    #     self.gravity:
    #         for chaining

    #     @TODO
    #     support more than 2 systems

    #     """
    #     # timestep
    #     set_timestep = True  # whether need to set the time-step
    #     if timestep is None:  # assume already set
    #         set_timestep = False
    #     elif isinstance(timestep, ScalarQuantity):  # a number
    #         pass
    #     elif isinstance(timestep, str):  # from a System in self
    #         timestep = self[timestep].gravity.parameters.timestep / 2.0
    #     elif isinstance(timestep, System):  # from an external system
    #         timestep = timestep.gravity.parameters.timestep / 2.0
    #     else:  # from the gravity code of an external system
    #         timestep = timestep.parameters.timestep / 2.0

    #     # construct list of systems
    #     syss = []  # blank list of systems
    #     for system in systems:  # iterate thru systems

    #         if isinstance(system, str):  # system in self
    #             syss.append(self[system].gravity)
    #         elif isinstance(system, System):  # external system
    #             syss.append(system.gravity)
    #         else:  # gravity code
    #             syss.append(system)

    #     # add systems and gravitational dependencies to bridge codes
    #     for i, system in enumerate(syss):
    #         # print(f"adding {system}: {dependencies.get(i, 'no dependencies')}")

    #         if i not in dependencies.keys():
    #             deps = set()
    #         elif dependencies[i] is None:
    #             deps = set()
    #         else:
    #             deps = syss[dependencies[i]]

    #             if not isinstance(deps, list):  # catch single values
    #                 deps = set([deps])
    #             else:
    #                 deps = set(deps)

    #         self.gravity.add_system(system, partners=deps)
    #     # /def

    #     # Set how often to update external potential
    #     if set_timestep:
    #         print("set timestep")
    #         self.gravity.timestep = timestep

    #     return self.gravity

    # # /def


# /class


##############################################################################
# END

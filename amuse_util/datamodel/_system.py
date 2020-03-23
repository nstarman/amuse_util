# -*- coding: utf-8 -*-

"""System(s).

Routine Listings
----------------
System
Systems

"""

__author__ = "Nathaniel Starkman"

__all__ = [
    "System",
    "Systems",
]


###############################################################################
# IMPORTS

# GENERAL

import itertools
from dataclasses import dataclass as _dataclass
from typing import Any, Optional, Union

# amuse
from amuse.couple import bridge
from amuse.datamodel.particles import Particles
from amuse.units.quantities import ScalarQuantity

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


# ------------------------------------------------------------------------


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

    def __init__(
        self,
        particles: Optional[Particles] = None,
        evolution: Optional[Any] = None,
        gravity: Optional[Any] = None,
        converter: Optional[Any] = None,
        channel_attrs: Optional[list] = None,
        **kw,
    ):
        """Dataclass post initialization.

        ensures the particles, evolution, and gravity codes are wrapped
        in AmuseContainer

        makes channels between the all the codes

        """
        self.particles = particles
        self.evolution = evolution
        self.gravity = gravity

        # util
        self.converter = converter
        self.channel_attrs = channel_attrs

        for k, v in kw.items():
            setattr(self, k, v)

        names = ["particles", "evolution", "gravity"]

        for n in names:
            if (self[n] is not None) and not isinstance(
                self[n], AmuseContainer
            ):
                self[n] = AmuseContainer(self[n], n)

        for n1, n2 in itertools.permutations(names, 2):
            # print("Making channels:")
            if (self[n1] is not None) and (self[n2] is not None):
                # print("\t", n1, "<->", n2)
                self[n1].add_channel_to(
                    self[n2], attributes=self.channel_attrs
                )
                self[n1].add_channel_from(
                    self[n2], attributes=self.channel_attrs
                )

    # /def

    # ---------------------------------------------------------------
    # Get / Set

    def __getitem__(self, key):
        """Getitem via getattr.

        supports key.subkey where `key` is in self and `subkey`
        is an item (dictionary get) of `self.key`

        """
        k0, *ks = key.split(".")
        if not ks:  # ks is empty
            return getattr(self, k0)
        else:  # pass to AmuseContainer
            return getattr(self, k0)[".".join(ks)]

    # /def

    def __setitem__(self, key, value):
        """Setitem via setattr."""
        setattr(self, key, value)

    # /def

    # ---------------------------------------------------------------
    # Representations

    def __repr__(self):
        """Modified representation.

        Returns
        -------
        str
            standard representation
            with particles, evolution, and gravity fields.

        """
        _repr = "\n".join(
            (
                super().__repr__(),
                f"\tparticles: {self.particles.__repr__()}",
                f"\tevolution: {self.evolution.__repr__()}",
                f"\tgravity: {self.gravity.__repr__()}",
            )
        )
        return _repr

    # /def

    # ---------------------------------------------------------------
    # Serialize

    def __reduce_ex__(self, protocol):
        """Reduce for serialization.

        TODO FIX
        """
        redx = (
            _system_reconstructor,
            ([], self.__dict__,),
            None,
            None,
            None,
        )
        return redx

    # /def

    # ---------------------------------------------------------------
    # Plot

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
    internal_bridges: dict:
        `bridges` arguments into ``self.bridge_internal_systems``
    timestep:
        Bridge timestep
    use_threading: bool (default False)
        whether to use threading
    **systems: {name: System} pairs

    Returns
    -------
    systems: Systems instance

    TODO
    ----
    make internal_bridges update when add more bridges,
        currently static on initialization

    """

    def __init__(
        self,
        *,
        internal_bridges={},
        timestep=None,
        use_threading=True,
        **systems,
    ):
        """Initialize a Systems instance."""
        super().__init__()
        # save inputs
        self._use_threading = use_threading
        self._init_internal_bridges = internal_bridges
        self._init_timestep = timestep

        self.gravity = bridge.Bridge(use_threading=use_threading)
        self.system_list = set()  # system list set

        for k, v in systems.items():
            self[k] = v

        # now do bridges. needed to create systems first
        self.bridge_internal_systems(
            bridges=internal_bridges, timestep=timestep
        )

    # /def

    # ---------------------------------------------------------------
    # Representations

    def __repr__(self):
        """Modified representation.

        Returns
        -------
        str
            standard representation
            with gravity and all user-defined fields.
            TODO full field detail

        """
        _repr = "\n".join(
            (
                super().__repr__(),
                f"\t{self.gravity}",
                *[f"\t{k}" for k in self.system_list],
            )
        )
        return _repr

    # /def

    # ---------------------------------------------------------------
    # Serialize

    def __getitem__(self, key):
        """Getitem via getattr.

        supports key.subkey where `key` is in self and `subkey`
        is an item (dictionary get) of `self.key`

        """
        k0, *ks = key.split(".")
        if not ks:  # ks is empty
            return getattr(self, k0)
        else:  # pass to System
            return getattr(self, k0)[".".join(ks)]

    # /def

    def __setitem__(self, key, value):
        """Setitem via setattr."""
        setattr(self, key, value)
        self.system_list.add(key)

    # /def

    # ---------------------------------------------------------------
    # Bridges

    def bridge_internal_systems(
        self, bridges: dict, timestep: Union[ScalarQuantity, str, None] = None
    ):
        """Bridge gravity within Systems.

        Works on gravity only.

        Parameters
        ----------
        bridges: dict
            ex: "cluster.gravity": ["galaxy.gravity", "cdf_code"]
        timestep: None or amuse ScalarQuantity or str
            if None, assume already have set a time step
            if str, self[timestep]

        Returns
        -------
        self.gravity:
            for chaining

        """
        # add bridges
        # iterate through dictionary of bridges {'system': []}
        for name, deps in bridges.items():
            if not isinstance(deps, (list, tuple)):
                raise TypeError(f"args for {name} must be a list/tuple")
            self.gravity.add_system(
                self[name], partners=[self[dep] for dep in deps]
            )

        if timestep is not None:
            if isinstance(timestep, str):
                timestep = self[timestep]
            self.gravity.timestep = timestep  # set gravity
            self._init_timestep = timestep  # also (re)save timestep

        return self.gravity

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

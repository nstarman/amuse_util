# -*- coding: utf-8 -*-

"""Container for Amuse Particles objects.

Routine Listings
----------------
AttrDict
proxy_reconstructor
amuse_socket_reconstructor
AmuseContainer


"""

__author__ = "Nathaniel Starkman"

__all__ = [
    "AmuseContainer",
    "AttrDict",
    "proxy_reconstructor",
    "amuse_socket_reconstructor",
]


###############################################################################
# IMPORTS

# GENERAL

import warnings
import wrapt
from inspect import BoundArguments

# typing
from typing import Optional, Any

# amuse
from amuse.datamodel.particles import Particles


###############################################################################
# CODE
###############################################################################


class AttrDict(dict):
    """A dictionary with getattr and setattr methods."""

    def __getattr__(self, name: str):
        """Redirects to __getitem__."""
        return super().__getitem__(name)

    # /def

    def __setattr__(self, name: str, value: Any):
        """Redirects to __setitem__."""
        return super().__setitem__(name, value)

    # /def


# /class


def proxy_reconstructor(cls, wrapped_reduce: tuple, name: str):
    """Reconstruct a proxy object.

    Parameters
    ----------
    cls: type
        the ObjectProxy class object
    wrapped_reduce: tuple
        the reduced form of the proxied oject.
        form (reconstructor, (cls, *args), {state})
    name: str
        the name of the object
        generally 'particles', 'evolution', or 'gravity'

    Returns
    -------
    container: wrapt.ObjectProxy
        the reconstruced proxied object

    """
    # break apart wrapped_reduce
    if len(wrapped_reduce) == 2:
        proxied_obj_reconstructor, args = wrapped_reduce
        state = None
        app = None
        kv = None
        objstate = None
    else:
        proxied_obj_reconstructor, args, *rest = wrapped_reduce

        state = None
        app = None
        kv = None
        objstate = None

        if len(rest) > 0:
            state = rest[0]
        if len(rest) > 1:
            app = rest[1]
        if len(rest) > 2:
            kv = rest[2]
        if len(rest) > 3:
            objstate = rest[3]  # TODO

    # make wrapped object
    wrapped = proxied_obj_reconstructor(*args)

    if state is not None:
        try:
            wrapped.__setstate__
        except AttributeError:
            wrapped.__dict__.update(state)
        else:
            wrapped.__setstate__(state)

    if app is not None:
        for item in app:
            wrapped.append(item)

    if kv is not None:
        for k, v in kv:
            wrapped[k] = v

    # put into container
    container = cls(wrapped, name)

    return container


# /def


def amuse_socket_reconstructor(
    cls,
    ba: BoundArguments,
    particles: Particles,
    # can_make_blank: Optional[int] = None,
):
    """Reconstruct an AMUSE Socket Object.

    Like an SSE or BHTree object.
    Must have particles.

    Parameters
    ----------
    cls: type
        the function or class
        like SSE() or BHTree()
    ba: BoundArguments
        BoundArguments for `cls`
    particles: Particles
        the particles. will be assigned to obj.

    """
    obj = cls(*ba.args, **ba.kwargs)

    obj.particles.remove_particles(obj.particles)
    obj.particles.add_particles(particles)

    return obj


# /def


##########################################################################


class AmuseContainer(wrapt.ObjectProxy):
    """AmuseContainer for object that can hold channels.

    AmuseContainer is built with ``wrapt.ObjectProxy``, so looks like
    the object it holds, with a few added attributes and methods.

    Attributes
    ----------
    name: str
        the name of the object
    channel_to: AttrDict
        dictionary of channels to other objects
    channel_from: AttrDict
        dictionary of channels from other objects
    _inputs: BoundArguments or None
        arguments used to create the contained object
        needed for reconstructing socket objects like BHTree
    # _min_p_recon: int, optional
    #     whether can make a blank version of the wrapped object
    #     if None (default), make full
    #     if int, the minimum number of particles to create before replacing

    """

    def __init__(
        self, wrapped: Any, name: str, _inputs: Optional[BoundArguments] = None
    ):
        """Initialize container, adding name and channel info."""
        # look like wrap
        super().__init__(wrapped)
        # add name
        self._self_name = name
        # add channel hooks
        self._self_channel_to = AttrDict()
        self._self_channel_from = AttrDict()
        # add wrapped object's bound-argument inputs
        self._self__inputs = _inputs

    # /def

    def __getattr__(self, name: str):
        """Redefine __getattr__ so that also checks container."""
        # first check original object
        try:
            return super().__getattr__(name)
        # now try getting from container
        except AttributeError:
            try:  # try container
                return object.__getattribute__(self, "_self_" + name)
            except AttributeError:  # not in container either
                raise AttributeError(
                    "neither container nor object has attribute"
                )

    # /def

    # ---------------------------------------------------------------
    # Get / Set

    def __getitem__(self, key):
        """Getitem via getattr.

        supports key.subkey where `key` is in self and `subkey`
        is an attribute of `self.key`

        """
        try:
            ks = key.split(".")
            if len(ks) == 1:  # only one key
                return getattr(self, ks[0])
            else:  # pass to AmuseContainer
                return getattr(getattr(self, ks[0]), ks[1])
        except AttributeError:  # it's not a string
            return getattr(self, key)

    # /def

    # ---------------------------------------------------------------
    # Serialize

    def __reduce_ex__(self, protocol: int):
        """Reduce method for pickling.

        Needed b/c Proxies can't normally pickle.
        works using `proxy_reconstructor`, which remakes the wrapped
        object and then applies the proxy.

        Parameters
        ----------
        protocol
            passed to wrapped object's ``__reduce_ex__``

        Returns
        -------
        redx: tuple
            reduced form of self
            uses ``proxy_reconstructor``

        Notes
        -----
        first tries the internal serialization technique (__reduce_ex__)
        if that fails, which it will for socket objects like BHTree,
        fall back to my socket reconstructor function
        this requires that the AmuseContainer has the input information
        to the wrapped object as the serialization works by making the object
        then adding the particles.
        The particles need to be copied. This can be memory intensive but also
        means that any links MUST BE RE-ESTABLISHED.

        See Also
        --------
        proxy_reconstructor
        amuse_socket_reconstructor

        TODO
        ----
        pickle large items with hdf5 in a folder format

        """
        # Serialize the wrapped object
        if self._inputs is None:
            internal_redx = self.__wrapped__.__reduce_ex__(protocol)

        else:
            internal_redx = (
                amuse_socket_reconstructor,  # reconstructor
                (
                    self.__class__,
                    self._inputs,
                    self.particles.copy(),
                ),  # (cls, *args)
                None,  # state
            )

        redx = (
            proxy_reconstructor,
            (
                AmuseContainer,  # proxy
                internal_redx,  # internal pickler
                self._self_name,  # name
            ),
            None,
            None,
            None,
        )
        return redx

    # /def

    # ---------------------------------------------------------------
    # channels

    def add_channel_to(
        self,
        going_to: Any,
        name: Optional[str] = None,
        attributes: Optional[list] = None,
    ):
        """Channels to things.

        Parameters
        ----------
        going_to: object
            something that  supports AMUSE channels
        name: str, optional
            needed if `going_to` is not an AmuseContainer with a name
            overrides the name in `going_to`, if `going_to` has a ``.name``
        attributes: list, optional
            attributes in the channel
            see AMUSE docs

        Returns
        -------
        channel
            see AMUSE docs

        Raises
        ------
        ValueError
            if `going_to` is not an AmuseContainer and `name` is None
            or if self.name == `name`
        warning
            if `name` (or `going_to.name`) is already in ``.channel_to``

        """
        # 1st check if thing is my object
        if name is None:
            if hasattr(going_to, "name"):  # System has name
                name = going_to.name
            else:
                raise ValueError("need to pass name")

        # next check that NOT making channel to same thing
        if self.name == name:
            raise ValueError("cannot make channel to self")
        elif name in self.channel_to.keys():
            warnings.warn(f"{name} already exists in channels. overwriting.")

        # making channel
        if hasattr(self, "particles") and hasattr(going_to, "particles"):
            self.channel_to[name] = self.particles.new_channel_to(
                going_to.particles, attributes=attributes
            )
        elif hasattr(self, "particles"):  # (so going_to does not)
            self.channel_to[name] = self.particles.new_channel_to(
                going_to, attributes=attributes
            )
        elif hasattr(going_to, "particles"):  # (so self does not)
            self.channel_to[name] = self.new_channel_to(
                going_to.particles, attributes=attributes
            )
        else:  # (neither self nor going_to has .particles)
            self.channel_to[name] = self.new_channel_to(
                going_to, attributes=attributes
            )

        return self.channel_to[name]

    # /def

    def add_channel_from(
        self,
        going_from: Any,
        name: Optional[str] = None,
        attributes: Optional[list] = None,
    ):
        """Channels from things.

        Parameters
        ----------
        going_from: object
            something that  supports AMUSE channels
        name: str, optional
            needed if `going_from` is not an AmuseContainer with a name
            overrides the name in `going_from`, if `going_from` has a ``.name``
        attributes: list, optional
            attributes in the channel
            see AMUSE docs

        Returns
        -------
        channel
            see AMUSE docs

        Raises
        ------
        ValueError
            if `going_from` is not an AmuseContainer and `name` is None
            or if self.name == `name`
        warning
            if `name` (or `going_from.name`) is already in ``.channel_from``

        """
        # 1st check if thing is my object
        if name is None:
            if hasattr(going_from, "name"):
                name = going_from.name
            else:
                raise ValueError("need to pass name")

        # next check that not making channel to same thing
        if self.name == name:
            raise ValueError("cannot make channel to self")
        elif name in self.channel_from.keys():
            warnings.warn(f"{name} already exists in channels. overwriting.")

        # making channels
        if hasattr(self, "particles") and hasattr(going_from, "particles"):
            self.channel_from[name] = going_from.particles.new_channel_to(
                self.particles, attributes=attributes
            )
        elif hasattr(self, "particles"):  # (so going_from does not)
            self.channel_from[name] = going_from.new_channel_to(
                self.particles, attributes=attributes
            )
        elif hasattr(going_from, "particles"):  # (so self does not)
            self.channel_from[name] = going_from.particles.new_channel_to(
                self, attributes=attributes
            )
        else:  # (neither self nor going_from has .particles)
            self.channel_from[name] = going_from.new_channel_to(
                self, attributes=attributes
            )

        return self.channel_from[name]

    # /def


##############################################################################
# END

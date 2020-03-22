# -*- coding: utf-8 -*-

# Docstring and Metadata
"""Initialize system method.

initialize_system returns

system: datamodel.System object
    a dataclass object with parameters
        - particles
        - evolution
        - gravity
        - converter
        - channel_attrs
    will try to automatically make channels to/from all things
    the amuse particles / evolution / gravity classes are proxied
    in a datamodel.Container that adds .name, .channel_to/from

"""

__author__ = "Nathaniel Starkman"

__all__ = ["initialize_system", "recreate_system"]

##############################################################################
# IMPORTS

# GENERAAL
import copy
import inspect
import numpy as np
from types import FunctionType

# amuse
from amuse.units import units as u, nbody_system
from amuse.datamodel import Particles

# CUSTOM
# from astroPHD import LogFile
# from astroPHD.decorators import store_function_input


# PROJECT-SPECIFIC
from ..datamodel import AmuseContainer, System

from ..utils import store_function_input


##############################################################################
# PARAMETERS

# typing
func_or_cls = (FunctionType, type)


##############################################################################
# CODE
##############################################################################


@store_function_input(store_inputs=True)
def initialize_particles(
    number_of_particles: int,
    *,  # must use kwargs
    # for IMF
    imf_func: (bool, FunctionType) = True,
    imf_args: list = [],
    imf_kwargs: dict = {},
    # for distribution function
    distr_func: (bool, FunctionType) = True,
    distr_args: list = [],
    distr_kwargs: dict = {},
    # object properties
    Rvirial: u.parsec = 10 | u.parsec,
    position: u.kpc = [0, 0, 0] | u.kpc,
    velocity: u.kms = [0, 0, 0] | u.kms,
    obj_radius: u.AU = 0 | u.AU,
    # util
    converter: (nbody_system.nbody_to_si, None) = None,
    random: (bool, int, np.random.RandomState) = True,
    _scale_to_standard: (bool, dict) = True,
    # logging
    # logger: LogFile = _LOGFILE,
    # verbose: (int, None) = None,
):
    """Create objs in a system.

    Parameters
    ----------
    number_of_particles: int or Particles
        if int, number of particles in the system
        if Particles instance, then the `imf_`, `distr_,
        and kwargs before `evln_func` are ignored

    imf_func: function
        function for initial mass function
        ex) new_kroupa_mass_distribution
        signature of function should be
        func(number_of_particles, *imf_args, random=random, **imf_kwargs)
    imf_args: list, optional
        the arguments for `imf_func`
        (default [])
    imf_kwargs: dictionary, optional
        the kwargs for `imf_func`
        (default {})

    distr_func: function
        function for object spatial distribution
        ex) new_plummer_model
        signature of function should be
        func(number_of_particles, *distr_args,
             convert_nbody=converter, **distr_kwargs)
    distr_args: list, optional
        the arguments for `distr_func`
        (default [])
    distr_kwargs: dictionary, optional
        the kwargs for `distr_func`
        (default {})

    Rvirial: distance quantity, optional
        the virial radius of the system
        (default 10 pc)
    position: distance 3-list quantity, optional
        the position of the system in GC coordinates
        (default [0,0,0] pc)
    velocity: velocity 3-list quantity, optional
        the velocity of the system in GC coordinates
        (default [0,0,0] kms)
    obj_radius: distance quantity, optional
        the radius of the individual objects
        (default 0 AU)

    _num_particles_reconstruct: int, optional
        minimum number of particles that can be used to reconstruct
        these particles if in a system. Proxy object like BHTree
        cannot be pickled and need to be manually reconstructed.
        If the reconstruction is expensive, making a minimum copy
        and replacing the pickled particles is much faster.

    Returns
    -------
    Particles
    converter

    References
    ----------
    modelled from:
    https://github.com/amusecode/amuse/blob/master/examples/textbook/solar_cluster_in_live_galaxy.py
    https://github.com/webbjj/galpy_profile/blob/master/test_cluster.py

    """
    # -------------------------------------------
    # checking types

    if not isinstance(imf_func, FunctionType):
        raise ValueError("Need to provide an IMF function")
    if not isinstance(distr_func, FunctionType):
        raise ValueError("Need to provide a distribution function")

    # -------------------------------------------

    # logger.report(
    #     "\nInitialize Particles Arguments:\n",
    #     f"number_of_particles: {number_of_particles}",
    #     # for IMF
    #     f"imf_func: {imf_func}",
    #     f"imf_args: {imf_args}",
    #     f"imf_kwargs: {imf_kwargs}",
    #     # for distribution function
    #     f"distr_func: {distr_func}",
    #     f"distr_args: {distr_args}",
    #     f"distr_kwargs: {distr_kwargs}",
    #     # object properties
    #     f"Rvirial: {Rvirial}",
    #     f"position: {position}",
    #     f"velocity: {velocity}",
    #     f"obj_radius: {obj_radius}",
    #     sep="\n    ",
    #     verbose=verbose,
    # )

    # ------------------------------------------------------------------------
    # Objects

    # ------------------------------------
    # create mass distribution

    # call initial mass function
    # accepts any arguments, normally :convert_nbody:
    # accepts any keyword arguments, except :random:
    masses = imf_func(
        number_of_particles, *imf_args, random=random, **imf_kwargs
    )

    # logger.report(
    #     "made masses",
    #     f"made masses, mean: {np.mean(masses)}, sum: {masses.sum()}",
    #     verbose=verbose,
    # )

    # ------------------------------------
    # converter
    # scaled to the system's total mass and virial radius
    if converter is None:
        converter = nbody_system.nbody_to_si(masses.sum(), Rvirial)

        # logger.report(
        #     f"made converter: {masses.sum(), Rvirial}", verbose=verbose,
        # )

    # ------------------------------------
    # Create system with origin at <0,0,0> and no net velocity
    # the system will be spatially distributed according to `distr_func`
    # also scale system to `converter` from unitless size

    # call distribution function
    objs = distr_func(
        number_of_particles,
        *distr_args,
        random=random,
        convert_nbody=converter,
        **distr_kwargs,
    )  # make objects from distribution

    # scale to standard
    if _scale_to_standard is not False:
        kw = {} if _scale_to_standard is True else _scale_to_standard
        objs.scale_to_standard(
            convert_nbody=converter, **kw,
        )

    # logger.report(
    #     "made objects",
    #     f"system has virial radius: {objs.virial_radius().in_(u.parsec)}",
    #     verbose=verbose,
    # )

    # add in the properties
    objs.radius = obj_radius  # object sizes (generally 0, pt size)
    objs.mass = masses  # masses

    # Place system in Galactocentric position
    objs.position += position
    objs.velocity += velocity

    # logger.report("added mean position & velocity to system", verbose=verbose)

    return objs, converter


# /def


# ----------------------------------------------------------------------------


@store_function_input(store_inputs=True)
def initialize_system(
    number_of_particles: (int, Particles),
    *,  # must use kwargs
    # for IMF
    imf_func: (bool, FunctionType) = True,
    imf_args: list = [],
    imf_kwargs: dict = {},
    # for distribution function
    distr_func: (bool, FunctionType) = True,
    distr_args: list = [],
    distr_kwargs: dict = {},
    # object properties
    Rvirial: u.parsec = 10 | u.parsec,
    position: u.kpc = [0, 0, 0] | u.kpc,
    velocity: u.kms = [0, 0, 0] | u.kms,
    obj_radius: u.AU = 0 | u.AU,  # size of each object
    # for evolution
    evln_func: (bool, FunctionType) = False,
    evln_kwargs: dict = {},
    # for gravity
    gravity_func: (bool, FunctionType) = False,
    gravity_args: list = [],
    gravity_kwargs: dict = {},
    smoothing_length: u.parsec = 0.0 | u.parsec,
    opening_angle: float = 0.6,
    number_of_workers: int = 8,
    use_self_gravity: int = 1,
    # util
    converter: (nbody_system.nbody_to_si, None) = None,
    timestep: u.Myr = 1.0 | u.Myr,
    random=True,
    _num_particles_reconstruct: (int, None) = None,
    # channel_attrs=None,
    # logging
    # logger: LogFile = _LOGFILE,
    # verbose: (int, None) = None,
    # debugging
    _scale_to_standard: bool = True,
):
    """Create objs in a system.

    modelled from:
    https://github.com/amusecode/amuse/blob/master/examples/textbook/solar_cluster_in_live_galaxy.py
    https://github.com/webbjj/galpy_profile/blob/master/test_cluster.py

    Parameters
    ----------
    number_of_particles: int or Particles
        if int, number of particles in the system
        if Particles instance, then the `imf_`, `distr_,
        and kwargs before `evln_func` are ignored

    imf_func: function
        function for initial mass function
        ex) new_kroupa_mass_distribution
        signature of function should be
        func(number_of_particles, *imf_args, random=random, **imf_kwargs)
    imf_args: list, optional (default [])
        the arguments for `imf_func`
    imf_kwargs: dictionary, optional (default {})
        the kwargs for `imf_func`

    distr_func: function
        function for object spatial distribution
        ex) new_plummer_model
        signature of function should be
        func(number_of_particles, *distr_args,
             convert_nbody=converter, **distr_kwargs)
    distr_args: list, optional (default [])
        the arguments for `distr_func`
    distr_kwargs: dictionary, optional (default {})
        the kwargs for `distr_func`

    Rvirial: distance quantity, optional (default 10 pc)
        the virial radius of the system
    position: distance 3-list quantity, optional (default [0,0,0] pc)
        the position of the system in GC coordinates
    velocity: velocity 3-list quantity, optional (default [0,0,0] kms)
        the velocity of the system in GC coordinates
    obj_radius: distance quantity, optional (default 0 AU)
        the radius of the individual objects

    evln_func: function or False, optional  (default False)
        object evolution function
    imf_kwargs: dictionary, optional
        the kwargs for `evln_func`

    gravity_func: function or False, optional  (default False)
        gravity code
        signature of function should be
        func(converter, *`gravity_args`,
             number_of_workers=`number_of_workers`,
             **`gravity_kwargs`)
    gravity_args: list, optional
        the arguments for `gravity_func`
    gravity_kwargs: dictionary, optional
        the kwargs for `gravity_func`
    smoothing_length: distance quantity, optional (default 0 pc)
        the smoothing length used in scaling and the gravity
    opening_angle: float, optional  (default 0.6)
    number_of_workers: int, optional  (default 1)
        number of gravity workers
    use_self_gravity: int
        flag for usage of self gravity, 1 or 0 (true or false)
    converter: nbody_to_si, optional
        nbody converter, takes total mass and virial radius
        calculated if not provided
    timestep: time quantity, optional  (default 1 Myr)
        the timestep for evolving the gravity and stellar evolution
    random: True or random number generator, optional (default True)
        ex: np.random.default_rng(seed=0)
        will default to random seed
    # channel_attrs: list or None, optional (default None)

    Returns
    -------
    system: datamodel.System
        a dataclass object with parameters
            - particles
            - evolution
            - gravity
            - converter
            - channel_attrs
        will try to automatically make channels to/from all things
        the amuse particles / evolution / gravity classes are proxied
        in a datamodel.Container that adds .name, .channel_to/from

    """
    # -------------------------------------------
    # checking types

    if evln_func is True:
        raise ValueError("Need to specify an evolution code")
    elif (evln_func is False) | isinstance(evln_func, func_or_cls):
        pass
    else:
        raise ValueError("Need to provide an evolution code")

    if gravity_func is True:
        raise ValueError("Need to specify a gravity code")
    elif (gravity_func is False) | isinstance(gravity_func, func_or_cls):
        pass
    else:
        raise ValueError("Need to provide a gravity code")

    # -------------------------------------------

    # logger.newsection("Initialize System")
    # logger.report(
    #     "\nInitialize System Arguments:\n",
    #     f"number_of_particles: {number_of_particles}",
    #     # for IMF
    #     f"imf_func: {imf_func}",
    #     f"imf_args: {imf_args}",
    #     f"imf_kwargs: {imf_kwargs}",
    #     # for distribution function
    #     f"distr_func: {distr_func}",
    #     f"distr_args: {distr_args}",
    #     f"distr_kwargs: {distr_kwargs}",
    #     # object properties
    #     f"Rvirial: {Rvirial}",
    #     f"position: {position}",
    #     f"velocity: {velocity}",
    #     f"obj_radius: {obj_radius}",
    #     # for evolution
    #     f"evln_func: {evln_func}",
    #     f"evln_kwargs: {evln_kwargs}",
    #     # for gravity
    #     f"gravity_func: {gravity_func}",
    #     f"gravity_args: {gravity_args}",
    #     f"gravity_kwargs: {gravity_kwargs}",
    #     f"smoothing_length: {smoothing_length}",
    #     f"opening_angle: {opening_angle}",
    #     f"number_of_workers: {number_of_workers}",
    #     f"use_self_gravity: {use_self_gravity}",
    #     # util
    #     f"converter: {converter}",
    #     f"timestep: {timestep}",
    #     f"random: {random}",
    #     f"channel_attrs: {'deprecated'}",
    #     # debugging
    #     f"_scale_to_standard: {_scale_to_standard}",
    #     sep="\n    ",
    #     verbose=verbose,
    # )

    # ------------------------------------------------------------------------
    # Objects

    if isinstance(number_of_particles, Particles):

        objs = number_of_particles

        number_of_particles = len(objs)

        if converter is None:
            raise ValueError(
                "need converter if `number_of_particles` is Particles"
            )

    elif isinstance(number_of_particles, int):

        objs, converter = initialize_particles(
            number_of_particles,
            # for IMF
            imf_func=imf_func,
            imf_args=imf_args,
            imf_kwargs=imf_kwargs,
            # for distribution function
            distr_func=distr_func,
            distr_args=distr_args,
            distr_kwargs=distr_kwargs,
            # object properties
            Rvirial=Rvirial,
            position=position,
            velocity=velocity,
            obj_radius=obj_radius,
            # util
            converter=converter,
            random=random,
            _scale_to_standard=_scale_to_standard,
            # logging
            # logger=logger,
            # verbose=verbose,
            # decorator
            store_inputs=False,
        )

    else:
        raise TypeError("number_of_particles is not int nor Particles")

    # ------------------------------------------------------------------------
    # evln evolution

    if evln_func is not False:  # edge cases considered above

        # call evolution function
        # only uses key-word arguments
        _sig = inspect.signature(evln_func)
        ba = _sig.bind_partial(**evln_kwargs)

        _evln = evln_func(*ba.args, **ba.kwargs)

        if _num_particles_reconstruct is not None:
            ba.arguments["number_of_particles"] = _num_particles_reconstruct
        evln = AmuseContainer(_evln, "evolution", _inputs=ba)

        # add particles to evolution code
        evln.particles.add_particles(objs)

        # logger.report(
        #     "added stellar evolution",
        #     f"added stellar evolution from function {evln_func}",
        #     verbose=verbose,
        # )

        # particles should be the master list
        # the evolution, if it exists needs to update the particles properties
        channel = evln.particles.new_channel_to(objs)
        channel.copy()

    else:
        evln = None

        # logger.report("no evolution", verbose=verbose)

    # ------------------------------------------------------------------------
    # gravity

    if gravity_func is not False:  # edge cases considered above

        # create gravity code
        _sig = inspect.signature(gravity_func)
        ba = _sig.bind_partial(
            converter,
            *gravity_args,
            number_of_workers=number_of_workers,
            **gravity_kwargs,
        )

        _gravity = gravity_func(*ba.args, **ba.kwargs)

        if _num_particles_reconstruct is not None:
            ba.arguments["number_of_particles"] = _num_particles_reconstruct
        gravity = AmuseContainer(_gravity, "gravity", _inputs=ba)

        # gravity = gravity_func(
        #     converter,
        #     *gravity_args,
        #     number_of_workers=number_of_workers,
        #     **gravity_kwargs,
        # )

        # add in gravity parameters
        gravity.parameters.use_self_gravity = use_self_gravity
        gravity.parameters.epsilon_squared = smoothing_length ** 2
        gravity.parameters.opening_angle = opening_angle
        gravity.parameters.timestep = timestep

        # add particles to gravity code
        gravity.particles.add_particles(objs)

        # BHtree, not sure if affects others
        gravity.commit_particles()  # needed for BHTree. doesn't affect others?

        # logger.report(f"added gravity code {gravity_func}", verbose=verbose)

    else:

        gravity = None

        # logger.report("no no gravity", verbose=verbose)

    # -------------------------------------------------------------------------
    # make a system & store inputs for reproduction

    system = System(
        particles=objs, evolution=evln, gravity=gravity, converter=converter,
    )

    return system


# /def


##############################################################################


def recreate_system(
    ba: inspect.BoundArguments,
    particles: Particles,
    converter: (nbody_system.nbody_to_si, None) = None,
):
    """Reinitialize a System made with `initialize_system`.

    Parameters
    ----------
    ba: BoundArguments
    particles: Particles
        should have everything, including for evolution
    converter:
        needed if 'converter' not in `ba.kwargs`
        if provided, will overwrite built-in, if exists

    Notes
    -----
    A good way to do this is to have `particles` be the set with all the
    up-to-date information. Then there is no need to pass `evolution` or
    `gravity`.

    """
    # change number_of_particles arg to the `particles` set
    # this will then initialize the system off the existing particles
    # skipping the IMF and distribution function steps
    ba = copy.copy(ba)

    ba.arguments["number_of_particles"] = particles

    if converter is not None:
        ba.arguments["converter"] = converter

    # make the system
    system = initialize_system(*ba.args, store_inputs=False, **ba.kwargs)

    return system


# /def


##############################################################################
# END

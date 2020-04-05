# -*- coding: utf-8 -*-

"""Star Cluster Code.

Routine Listings
----------------
initialize_star_cluster
separate_bound_unbound

_initialize_star_cluster_bound
_initialize_star_cluster_unbound


"""

__author__ = "Nathaniel Starkman"

__all__ = [
    "SC_Tuple",
    "initialize_star_cluster",
    "separate_bound_unbound",
    "_initialize_star_cluster_bound",
    "_initialize_star_cluster_unbound",
]


###############################################################################
# IMPORTS

# GENERAL

import copy
from collections import namedtuple
from typing import Tuple, Any, Callable

import numpy as np

# amuse
import amuse.units.units as amu
from amuse.datamodel import Particles
from amuse.units import nbody_system

# typing
# from typing import Optional

# CUSTOM

try:
    import utilipy
except ImportError:
    from ..utils import store_function_input
    from ..utils._logging import LogFile
else:
    from utilipy import LogFile
    from utilipy.decorators import store_function_input

# PROJECT-SPECIFIC

from . import initialize_system
from ..units import to_amuse_decorator


###############################################################################
# PARAMETERS

SC_Tuple = namedtuple(
    "separate_bound_unbound",
    ["bound_system", "bound_inputs", "unbound_system", "unbound_inputs"],
)  # star cluster bound / unbound

_LOGFILE = LogFile(header=False)


###############################################################################
# CODE
###############################################################################


@store_function_input(store_inputs=True)
@to_amuse_decorator(  # ensure inputs are in AMUSE units
    arguments=[
        "Rvirial",
        "position",
        "velocity",
        "obj_radius",
        "smoothing_length",
        "timestep",
    ]
)
def initialize_star_cluster(
    number_of_particles: (int, Particles),
    *,  # must use kwargs
    # for IMF
    imf_func: (bool, Callable) = True,
    imf_args: list = [],
    imf_kwargs: dict = {},
    # for distribution function
    distr_func: (bool, Callable) = True,
    distr_args: list = [],
    distr_kwargs: dict = {},
    # object properties
    Rvirial: amu.parsec = 10 | amu.parsec,
    position: amu.kpc = [0, 0, 0] | amu.kpc,
    velocity: amu.kms = [0, 0, 0] | amu.kms,
    obj_radius: amu.AU = 0 | amu.AU,  # size of each object
    # for evolution
    evln_func: (bool, Callable) = False,
    evln_kwargs: dict = {},
    # for gravity
    gravity_func: (bool, Callable) = False,
    gravity_args: list = [],
    gravity_kwargs: dict = {},
    smoothing_length: amu.parsec = 0.0 | amu.parsec,
    opening_angle: float = 0.6,
    number_of_workers: int = 8,
    use_self_gravity: int = 1,
    # util
    converter: (nbody_system.nbody_to_si, None) = None,
    timestep: amu.Myr = 1.0 | amu.Myr,
    random=True,
    make_unbound=True,
    _num_particles_reconstruct: (int, None) = None,
    # logging
    logger: LogFile = _LOGFILE,
    verbose: (int, None) = None,
    # debugging
    _scale_to_standard: bool = True,
) -> SC_Tuple:
    """Docstring.

    Parameters
    ----------
    number_of_particles: int or Particles
        if int, number of particles in the system
        if Particles instance, then the `imf_`, `distr_,
        and kwargs before `evln_func` are ignored

    imf_func: Callable
        function for initial mass function
        ex) new_kroupa_mass_distribution
        signature of function should be
        func(number_of_particles, *imf_args, random=random, **imf_kwargs)
    imf_args: list, optional
        the arguments for `imf_func`
    imf_kwargs: dict, optional
        the kwargs for `imf_func`

    distr_func: Callable
        function for object spatial distribution
        ex) new_plummer_model
        signature of function should be
        func(number_of_particles, *distr_args,
             convert_nbody=converter, **distr_kwargs)
    distr_args: list, optional
        the arguments for `distr_func`
    distr_kwargs: dict, optional
        the kwargs for `distr_func`

    Rvirial: Quantity, optional
        the virial radius of the system, default 10 pc
    position: Quantity, optional
        the position of the system in GC coordinates
        default [0,0,0] pc
    velocity: Quantity, optional
        the velocity of the system in GC coordinates
        (default [0,0,0] kms)
    obj_radius: Quantity, optional
        the radius of the individual objects,
        (default 0 AU)

    evln_func: Callable or False, optional
        (default False)
        object evolution function
    imf_kwargs: dict, optional
        the kwargs for `evln_func`

    gravity_func: Callable or False, optional
        gravity code
        signature of function should be
        func(converter, *`gravity_args`,
             number_of_workers=`number_of_workers`,
             **`gravity_kwargs`)
    gravity_args: list, optional
        the arguments for `gravity_func`
    gravity_kwargs: dict, optional
        the kwargs for `gravity_func`
    smoothing_length: Quantity, optional
        the smoothing length used in scaling and the gravity
        (default 0 pc)
    opening_angle: float, optional
        (default 0.6)
    number_of_workers: int, optional
        number of gravity workers
        (default 1)
    use_self_gravity: int
        flag for usage of self gravity, 1 or 0 (True or False)
    converter: nbody_to_si, optional
        nbody converter, takes total mass and virial radius
        calculated if not provided
    timestep: time quantity, optional
        the timestep for evolving the gravity and stellar evolution
        (default 1 Myr)

    random: True or Generator, optional
        (default True)
        ex: np.random.default_rng(seed=0)
        will default to random seed
    _scale_to_standard: bool
        whether to call ``scale_to_standard``

    Returns
    -------
    bound_system : System
    bound_inputs : BoundArguments
    unbound_system : System
    unbound_inputs : BoundArguments

    """
    # call helper function _initialize_star_cluster_bound
    # returns system, inputs to ``initialize_system``,
    # and inputs to _initialize_star_cluster_bound
    (system, inputs), _inputs = _initialize_star_cluster_bound(
        number_of_particles=number_of_particles,
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
        # for evolution
        evln_func=evln_func,
        evln_kwargs=evln_kwargs,
        # for gravity
        gravity_func=gravity_func,
        gravity_args=gravity_args,
        gravity_kwargs=gravity_kwargs,
        smoothing_length=smoothing_length,
        opening_angle=opening_angle,
        number_of_workers=number_of_workers,
        use_self_gravity=use_self_gravity,
        # util
        converter=converter,
        timestep=timestep,
        random=random,
        _num_particles_reconstruct=_num_particles_reconstruct,
        # logging
        logger=logger,
        verbose=verbose,
        # debugging
        _scale_to_standard=_scale_to_standard,
    )

    # Make Unbound
    if make_unbound:

        (unbound_system, unbound_inputs,) = _initialize_star_cluster_unbound(
            system, _inputs, bound_system_func=_initialize_star_cluster_bound
        )

    else:

        unbound_system, unbound_inputs = None, None

    return SC_Tuple(
        bound_system=system,
        bound_inputs=inputs,
        unbound_system=unbound_system,
        unbound_inputs=unbound_inputs,
    )


# /def

# -------------------------------------------------------------------


@store_function_input(store_inputs=True)
def _initialize_star_cluster_bound(
    number_of_particles: (int, Particles),
    *,  # must use kwargs
    # for IMF
    imf_func: (bool, Callable) = True,
    imf_args: list = [],
    imf_kwargs: dict = {},
    # for distribution function
    distr_func: (bool, Callable) = True,
    distr_args: list = [],
    distr_kwargs: dict = {},
    # object properties
    Rvirial: amu.parsec = 10 | amu.parsec,
    position: amu.kpc = [0, 0, 0] | amu.kpc,
    velocity: amu.kms = [0, 0, 0] | amu.kms,
    obj_radius: amu.AU = 0 | amu.AU,  # size of each object
    # for evolution
    evln_func: (bool, Callable) = False,
    evln_kwargs: dict = {},
    # for gravity
    gravity_func: (bool, Callable) = False,
    gravity_args: list = [],
    gravity_kwargs: dict = {},
    smoothing_length: amu.parsec = 0.0 | amu.parsec,
    opening_angle: float = 0.6,
    number_of_workers: int = 8,
    use_self_gravity: int = 1,
    # util
    converter: (nbody_system.nbody_to_si, None) = None,
    timestep: amu.Myr = 1.0 | amu.Myr,
    random=True,
    _num_particles_reconstruct: (int, None) = None,
    # logging
    logger: LogFile = _LOGFILE,
    verbose: (int, None) = None,
    # debugging
    _scale_to_standard: bool = True,
) -> tuple:
    """Bound Star Cluster Helper Function.

    Parameters
    ----------
    number_of_particles: int or Particles
        if int, number of particles in the system
        if Particles instance, then the `imf_`, `distr_,
        and kwargs before `evln_func` are ignored

    imf_func: Callable
        function for initial mass function
        ex) new_kroupa_mass_distribution
        signature of function should be
        func(number_of_particles, *imf_args, random=random, **imf_kwargs)
    imf_args: list, optional
        the arguments for `imf_func`
    imf_kwargs: dict, optional
        the kwargs for `imf_func`

    distr_func: Callable
        function for object spatial distribution
        ex) new_plummer_model
        signature of function should be
        func(number_of_particles, *distr_args,
             convert_nbody=converter, **distr_kwargs)
    distr_args: list, optional
        the arguments for `distr_func`
    distr_kwargs: dict, optional
        the kwargs for `distr_func`

    Rvirial: Quantity, optional
        the virial radius of the system, default 10 pc
    position: Quantity, optional
        the position of the system in GC coordinates
        default [0,0,0] pc
    velocity: Quantity, optional
        the velocity of the system in GC coordinates
        (default [0,0,0] kms)
    obj_radius: Quantity, optional
        the radius of the individual objects,
        (default 0 AU)

    evln_func: Callable or False, optional
        (default False)
        object evolution function
    imf_kwargs: dict, optional
        the kwargs for `evln_func`

    gravity_func: Callable or False, optional
        gravity code
        signature of function should be
        func(converter, *`gravity_args`,
             number_of_workers=`number_of_workers`,
             **`gravity_kwargs`)
    gravity_args: list, optional
        the arguments for `gravity_func`
    gravity_kwargs: dict, optional
        the kwargs for `gravity_func`
    smoothing_length: Quantity, optional
        the smoothing length used in scaling and the gravity
        (default 0 pc)
    opening_angle: float, optional
        (default 0.6)
    number_of_workers: int, optional
        number of gravity workers
        (default 1)
    use_self_gravity: int
        flag for usage of self gravity, 1 or 0 (True or False)
    converter: nbody_to_si, optional
        nbody converter, takes total mass and virial radius
        calculated if not provided
    timestep: time quantity, optional
        the timestep for evolving the gravity and stellar evolution
        (default 1 Myr)
    random: True or Generator, optional
        (default True)
        ex: np.random.default_rng(seed=0)
        will default to random seed
    _scale_to_standard: bool
        whether to call ``scale_to_standard``

    Returns
    -------
    system: datamodel.System
        a dataclass object with parameters
            - particles
            - evolution
            - gravity
            - converter
        will try to automatically make channels to/from all things
        the amuse particles / evolution / gravity classes are proxied
        in a datamodel.Container that adds .name, .channel_to/from

    inputs: BoundArguments
        the inputs to initialize_system

    """
    system, inputs = initialize_system(
        number_of_particles=number_of_particles,
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
        # for evolution
        evln_func=evln_func,
        evln_kwargs=evln_kwargs,
        # for gravity
        gravity_func=gravity_func,
        gravity_args=gravity_args,
        gravity_kwargs=gravity_kwargs,
        smoothing_length=smoothing_length,
        opening_angle=opening_angle,
        number_of_workers=number_of_workers,
        use_self_gravity=use_self_gravity,
        # util
        converter=converter,
        timestep=timestep,
        random=random,
        _num_particles_reconstruct=_num_particles_reconstruct,
        # logging
        logger=logger,
        verbose=verbose,
        # debugging
        _scale_to_standard=_scale_to_standard,
        # ensure get inputs
        store_inputs=True,
    )

    return system, inputs


# /def


# -------------------------------------------------------------------


def _initialize_star_cluster_unbound(
    bound_cluster,
    bound_inputs,
    bound_system_func: Callable,
    bound_radius=None,  # TODO should None be allowable?
    init_arg: Tuple[str, Any] = ("number_of_particles", 2),
) -> tuple:
    """Initialize unbound system from Contenta (2018) star cluster.

    function can be used separately from `_initialize_star_cluster_bound`,
    which often needs to be overwritten, while this function does not.



    Parameters
    ----------
    bound_cluster : System
    bound_inputs : BoundArguments
    bound_system_func: Callable
        function to construct bound cluster
        defaults to ``_initialize_star_cluster_bound``
    bound_radius : Any, optional
        the bound radius beyond which anything is considered unbound
    init_arg : Tuple[str, Any], optional
        the arguments needed to make the smallest version of the star cluster
        default is ("number_of_particles", 2)

    Returns
    -------
    unbound_cluster : System
    unbound_inputs : BoundArguments

    See Also
    --------
    ``initialize_star_cluster``
    ``_initialize_star_cluster_bound``

    """

    ba = copy.copy(bound_inputs)

    ba.arguments[init_arg[0]] = init_arg[1]  # need to initialize, then remove

    (unbound_cluster, unbound_inputs,) = bound_system_func(
        *ba.args, store_inputs=False, **ba.kwargs
    )

    # remove initialization particles
    unbound_cluster.particles.remove_particles(unbound_cluster.particles)
    if unbound_cluster.gravity is not None:
        unbound_cluster.gravity.particles.synchronize_to(
            unbound_cluster.particles
        )
    if unbound_cluster.evolution is not None:
        unbound_cluster.evolution.particles.synchronize_to(
            unbound_cluster.particles
        )

    # find bound & unbound particles
    if bound_radius is None:
        bound_radius = bound_cluster.bound_radius_cutoff
    bound_cluster, unbound_cluster, _ = separate_bound_unbound(
        bound_cluster, unbound_cluster, bound_cluster.bound_radius_cutoff
    )

    return unbound_cluster, unbound_inputs


# /def


# -------------------------------------------------------------------


def separate_bound_unbound(
    bound_cluster, unbound_cluster, bound_radius, cdf_code=None, converter=None
):
    """Separate Unbound particles from Bound Particles by a simple radius cut from the center of mass.

    Also reset cdf_code parameters, if provided.

    Parameters
    ----------
    bound_cluster : System
    unbound_cluster : System
    cdf_code : DynamicalFrictionClass, optional
        reset the total_mass
    converter: nbody_to_si
        needed if cdf_code is not None and bound_cluster
        does not have a method `.converter`

    """
    if converter is None:
        if hasattr(bound_cluster, "converter"):
            converter = bound_cluster.converter
        else:
            raise ValueError("need a `converter` if `cdf_code` is not None")
    else:
        pass  # assume in correct format

    # find density center
    # use instead of center of mass b/c less sensitive to filling the orbit donut.
    (
        densitycentre,
        coreradius,
        coredens,
    ) = bound_cluster.particles.densitycentre_coreradius_coredens(
        unit_converter=converter
    )
    rcenter = np.linalg.norm(densitycentre)

    # find bound & unbound particles
    # TODO plot what is and isn't bound
    # subset = bound_cluster.particles.bound_subset(
    #     unit_converter=bound_cluster.converter,
    # )
    unbound = bound_cluster.particles.select_array(
        lambda p: np.linalg.norm(p.value_in(amu.pc), axis=-1)
        > rcenter.value_in(amu.pc) + bound_radius.value_in(amu.pc),
        ["position"],
    )

    # move particles that left the cluster to `unbound`
    unbound_cluster.particles.add_particles(unbound)
    unbound_cluster.particles.channel_to.gravity.copy()

    # remove unbound particles from cluster so don't do dynfric
    bound_cluster.particles.remove_particles(unbound)
    try:
        bound_cluster.particles.synchronize_to(bound_cluster.gravity.particles)
    except AttributeError:  # don't have gravity
        pass
    try:
        bound_cluster.particles.synchronize_to(
            bound_cluster.evolution.particles
        )
    except AttributeError:  # don't have gravity
        pass

    # update mass of cluster in dynamical friction code
    if cdf_code is not None:

        rhm = bound_cluster.particles.LagrangianRadii(
            mf=[0.5], unit_converter=converter
        )[0]
        cdf_code.reset_parameters(
            GMs=bound_cluster.particles.total_mass(), rhm=rhm,  # half-mass
        )

    return bound_cluster, unbound_cluster, cdf_code


# /def


# -------------------------------------------------------------------

###############################################################################
# END

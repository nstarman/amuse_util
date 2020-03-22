# -*- coding: utf-8 -*-

# Docstring and Metadata
"""Model generator from CDF distribution."""

__author__ = "Nathaniel Starkman"

__all__ = [
    "MakeSphericalMassEnclModel",
    "new_sphericalmassencl_model",
    "new_sphericalmassencl_sphere",
]

##############################################################################
# IMPORTS

# General
import numpy as np
import numpy.random

from types import FunctionType

from amuse import datamodel
from amuse.units import units as u, nbody_system

# Project-Specific
from ..util import amuseify_array


###############################################################################
# CODE
###############################################################################


class MakeSphericalMassEnclModel(object):
    """Make Mass Model from Enclosed Mass Function.

    https://github.com/peterewills/itsample/blob/master/itsample.py
    https://codereview.stackexchange.com/questions/196286/inverse-transform-sampling
    https://usmanwardag.github.io/python/astronomy/2016/07/10/inverse-transform-sampling-with-python.html
    https://en.wikipedia.org/wiki/Inverse_transform_sampling
    """

    def __init__(
        self,
        number_of_particles: int,
        encl_mass_func: FunctionType,
        vel_potential,
        convert_nbody=None,
        radius_cutoff: u.kpc = 100 | u.kpc,
        # mass_cutoff=0.999,
        do_scale: bool = False,
        _vel_adj: float = 1.0,
        random_state=None,
        random=None,
    ):
        """Instantiate spatial / velocity distribution from enclosed mass.

        Parameters
        ----------
        number_of_particles: int
        encl_mass_func:
            signature:: encl_mass_func(R)
        vel_potential: amuse potential
            potential from which to sample for the velocities
            signature:: vel_potential(position)

        """
        super().__init__()

        self.number_of_particles = number_of_particles
        self.encl_mass_func = encl_mass_func
        self.vel_potential = vel_potential

        self.convert_nbody = convert_nbody
        self.radius_cutoff = radius_cutoff

        # _tmp = self.calculate_mass_cuttof_from_radius_cutoff(radius_cutoff))
        # self.mass_cutoff = min(mass_cutoff, _tmp
        self.do_scale = do_scale

        self._vel_adj = 1.0

        if random_state is not None:
            print("DO NOT USE RANDOM STATE")

        self.random_state = None

        if random is None:
            self.random = np.random
        else:
            self.random = random

        return

    # /def

    def calculate_radius_invcdf_distribution(self):
        """Calculate Radius from Inverted CDF sampling.

        Returns
        -------
        rs: ndarray
            radii for spatial distribution
            units of kpc

        Notes
        -----
        calls encl_mass_func(R), where R has units, and returns units

        """
        # radii, finely sampled to the cutoff
        R_unit = self.radius_cutoff.unit
        R = np.logspace(
            -5,
            np.log10(self.radius_cutoff.value_in(R_unit)),
            num=max(10 * self.number_of_particles, 1e6),
        )

        # get enclosed masses for inverse sampling
        sample_encl_mass = self.encl_mass_func(R | R_unit)
        m_unit = sample_encl_mass.unit  # the mass unit

        # inverse CDF transform  # TODO as separate function
        cdf = sample_encl_mass.value_in(m_unit)  # cdf
        us = self.random.uniform(cdf[0], cdf[-1], self.number_of_particles)
        # get sample from `us`
        rs = (
            np.array(
                [
                    R[np.argwhere(cdf == (cdf[(cdf - u) >= 0])[0])[0][0]]
                    for u in us
                ]
            )
            | R_unit
        )

        return rs

    # /def

    def new_positions_spherical_coordinates(self):
        """Create positions in spherical coordinates.

        Returns
        -------
        radius: distance quantity array
        theta, phi: ndarray
            units of radiuas

        """
        radius = self.calculate_radius_invcdf_distribution()

        num_objs = self.number_of_particles
        theta = np.arccos(self.random.uniform(-1.0, 1.0, size=num_objs))
        phi = self.random.uniform(0.0, np.pi ** 2, size=num_objs)

        return radius, theta, phi

    # /def

    def new_velocities_spherical_coordinates(self, x, y, z):
        """Create velocities in spherical coordinates.

        Parameters
        ----------
        x, y, z: array
            unit-ed

        Returns
        -------
        velocity : ndarray
            velocity
        theta, phi: ndarray
            units of radii

        Notes
        -----
        adjusts velocity be 1 / sqrt(_vel_adj), default of 1.

        """
        num_objs = self.number_of_particles

        # potential = km^2/s^2
        pot_at_pt = self.vel_potential.get_potential_at_point(0, x, y, z)
        # velocity
        velocity = amuseify_array(
            np.sqrt(np.abs(pot_at_pt) / self._vel_adj), to_unit=u.km / u.s
        )

        # random directions
        theta = np.arccos(self.random.uniform(-1.0, 1.0, size=num_objs))
        phi = self.random.uniform(0.0, np.pi ** 2, size=num_objs)

        return velocity, theta, phi

    # /def

    def coordinates_from_spherical(self, radius: u.kpc, theta, phi):
        """Convert Coordinates to Cartesian from Spherical Coords.

        Returns
        -------
        x, y, z: distance quantity ndarray
            units from `radius`

        """
        R_unit = radius.unit
        x = (radius * np.sin(theta) * np.cos(phi)).in_(R_unit)
        y = (radius * np.sin(theta) * np.sin(phi)).in_(R_unit)
        z = (radius * np.cos(theta)).in_(R_unit)

        return x, y, z

    # /def

    def new_model(self):
        """Make New Model.

        Returns
        -------
        m: ndarray
            masses
        position: ndarray
        velocity: ndarray

        """
        num_objs = self.number_of_particles

        m = np.zeros((num_objs, 1)) + (1.0 / num_objs)

        # position
        radius, theta, phi = self.new_positions_spherical_coordinates()
        x, y, z = self.coordinates_from_spherical(radius, theta, phi)
        position = np.hstack(x, y, z)

        # velocity
        radius, theta, phi = self.new_velocities_spherical_coordinates(x, y, z)
        x, y, z = self.coordinates_from_spherical(radius, theta, phi)
        velocity = np.hstack(x, y, z)

        return m, position, velocity

    # /def

    @property
    def result(self):
        """Result.

        Returns
        -------
        result: Particles
            Particles datamodel

        """
        num_objs = self.number_of_particles

        masses = np.ones(num_objs) / num_objs
        print("made masses")

        radius, theta, phi = self.new_positions_spherical_coordinates()
        x, y, z = self.coordinates_from_spherical(radius, theta, phi)
        print("made coordinates")

        speed, theta, phi = self.new_velocities_spherical_coordinates(x, y, z)
        vx, vy, vz = self.coordinates_from_spherical(speed, theta, phi)
        print("made velocities")

        # ---------------
        # build Particles
        result = datamodel.Particles(num_objs)
        # result.mass = nbody_system.mass.new_quantity(masses)
        result.mass = masses | u.MSun  # TODO FIX
        # spatial
        # result.x = nbody_system.length.new_quantity(x.reshape(num_objs))
        # result.y = nbody_system.length.new_quantity(y.reshape(num_objs))
        # result.z = nbody_system.length.new_quantity(z.reshape(num_objs))
        result.x = x.reshape(num_objs)
        result.y = y.reshape(num_objs)
        result.z = z.reshape(num_objs)
        # velocity
        # result.vx = nbody_system.speed.new_quantity(vx.reshape(num_objs))
        # result.vy = nbody_system.speed.new_quantity(vy.reshape(num_objs))
        # result.vz = nbody_system.speed.new_quantity(vz.reshape(num_objs))
        result.vx = vx.reshape(num_objs)
        result.vy = vy.reshape(num_objs)
        result.vz = vz.reshape(num_objs)
        # radius
        # result.radius = 0 | nbody_system.length
        result.radius = 0 | u.AU

        # ---------------

        result.move_to_center()

        if self.do_scale:
            result.scale_to_standard()

        if self.convert_nbody is not None:
            converter = self.convert_nbody.as_converter_from_si_to_generic()
            result = datamodel.ParticlesWithUnitsConverted(result, converter)
            result = result.copy()

        return result

    # /def


# /class


# --------------------------------------------------------------------------


def new_sphericalmassencl_model(
    number_of_particles: int,
    encl_mass_func: FunctionType,
    vel_potential,
    *list_arguments,
    **keyword_arguments
):
    """Create a sphere with the given number of particles.

    Returns a set of stars with equal mass and positions and velocities
    distributed to fit an enclosed mass function. The model is centered
    around the origin. Positions and velocities are optionally scaled such
    that the kinetic and potential energies are 0.25 and -0.5 in nbody-units,
    respectively.

    Parameters
    ----------
    number_of_particles: int
        Number of particles to include in the plummer sphere
    encl_mass_func: function
    vel_potential: amuse potential

    convert_nbody:  When given will convert the resulting set to SI units
    radius_cutoff: Cutoff value for the radius (defaults to 22.8042468)
    mass_cutoff: Mass percentage inside radius of 1
    do_scale: scale the result to exact nbody units (M=1, K=0.25, U=-0.5)

    """
    uc = MakeSphericalMassEnclModel(
        number_of_particles,
        encl_mass_func,
        vel_potential,
        *list_arguments,
        **keyword_arguments
    )
    return uc.result


new_sphericalmassencl_sphere = new_sphericalmassencl_model


##############################################################################
# END

#!/usr/bin/env python
# vim: set fileencoding=utf8 :
"""
skyvec2ins JWST Coronagraph Visibility Calculator

Developed by Chris Stark (cstark@stsci.edu), translated to Python
from IDL by Joseph Long (jlong@stsci.edu)

The allowed pointing of JWST leads to target visibility that depends on
ecliptic latitude, and the range of roll angles allowed depends on
solar elongation. The allowed PAs for a target can thus be a
complicated function of time. As a result, it can be difficult to 1)
understand the possible orientations of a given target on the detector,
2) determine the ideal roll angle offsets for multi-roll observations,
and 3) determine a group of targets that are simultaneously visible.
The JWST Coronagraph Visibility Calculator (CVC) was created to address
these issues and assist with creating APT programs and diagnosing
scheduling errors.

We stress that the CVC is designed to provide quick illustrations of
the possible observable orientations for a given target. As such, the
CVC rapidly approximates JWST’s pointing restrictions and does not
query the official JWST Proposal Constraint Generator (PCG). The CVC
does not include detailed pointing restrictions like Earth and Moon
avoidance, etc. Additionally, results may differ from official
constraints by a degree or so. Users should treat the results as close
approximations.
"""
from __future__ import print_function, division

import datetime

import numpy as np


def _wrap_to_2pi(scalar_or_arr):
    """Offsets angles outside 0 <= x <= 2 * pi to lie within the interval"""
    return np.asarray(scalar_or_arr) % (2 * np.pi)

def sun_ecliptic_longitude(start_date):
    """Compute ecliptic longitude of sun on start_date
    using equations from http://aa.usno.navy.mil/faq/docs/SunApprox.php
    """
    n_days = (start_date - datetime.datetime(2000, 1, 1, 12, 00)).days
    mean_longitude = 280.459 + 0.98564736 * n_days
    mean_anomaly = 357.529 + 0.98560028 * n_days
    mean_longitude %= 360.
    mean_anomaly %= 360.
    lambda_sun = (mean_longitude +
                  1.915 * np.sin(np.deg2rad(mean_anomaly)) +
                  0.020 * np.sin(2 * np.deg2rad(mean_anomaly)))
    return lambda_sun

def ad2lb(alpha_rad, delta_rad):
    """
    Converts celestial coordinates (ra, dec), i.e. (alpha, delta)
    to ecliptic coordinates (lambda, beta). All angles in radians.

    See Eq 3 in Leinert et al. 1998
    """
    obliq = _tenv(23, 26, 21.45)  # J2000 obliquity of Earth in degrees
    obliq *= np.pi / 180.0

    beta_rad = np.arcsin(np.sin(delta_rad) * np.cos(obliq) - np.cos(delta_rad) * np.sin(obliq) * np.sin(alpha_rad))
    coslambda = np.cos(alpha_rad) * np.cos(delta_rad) / np.cos(beta_rad)
    sinlambda = (np.sin(delta_rad) * np.sin(obliq) +
                 np.cos(delta_rad) * np.cos(obliq) * np.sin(alpha_rad)) / np.cos(beta_rad)
    lambda_rad = np.arctan2(sinlambda, coslambda)  # make sure we get the right quadrant
    lambda_rad = _wrap_to_2pi(lambda_rad)
    return lambda_rad, beta_rad


def lb2ei(lmlsun, beta):
    """Convert ecliptic coordinates (lambda-lambda_sun, beta) to
    alternative ecliptic coordinates (epsilon, i). All angles in radians.

    See Eq 11 in Leinert et al. 1998
    """
    # convert to an elongation in radians (see Eq 11 in Leinert et al. 1998)
    elong = np.arccos(np.cos(lmlsun) * np.cos(beta))
    cosinc = np.cos(beta) * np.sin(lmlsun) / np.sin(elong)
    sininc = np.sin(beta) / np.sin(elong)
    inc = np.arctan2(sininc, cosinc)  # make sure we get the right quadrant
    inc = _wrap_to_2pi(inc)
    # make inc lie between 0 - 2*pi
    # j = np.where(inc < 0.)
    # if len(j[0]) != 0:
    #     inc[j] += 2*np.pi
    # j = np.where(inc > 2*np.pi)
    # if len(j[0]) != 0:
    #     inc[j] -= 2*np.pi
    return elong, inc


def ei2lb(elong, inc):
    """Convert alternative ecliptic coordinates (epsilon, i) to
    ecliptic coordinates (lambda-lambda_sun, beta). All angles in radians.

    See Eq 12 in Leinert et al. 1998
    """
    beta = np.arcsin(np.sin(inc) * np.sin(elong))
    coslmlsun = np.cos(elong) / np.cos(beta)
    sinlmlsun = np.cos(inc) * np.sin(elong) / np.cos(beta)
    lmlsun = np.arctan2(sinlmlsun, coslmlsun)
    return lmlsun, beta


def lb2ad(lambda_rad, beta_rad):
    """Converts ecliptic coordinates (lambda, beta) to
    celestial coordinates (ra, dec), i.e. (alpha, delta).
    All angles in radians.

    See Eq 4 in Leinert et al. 1998
    """

    obliq = _tenv(23, 26, 21.45)  # J2000 obliquity of Earth in degrees
    obliq *= np.pi / 180.

    delta = np.arcsin(np.sin(beta_rad) * np.cos(obliq) + np.cos(beta_rad) * np.sin(obliq) * np.sin(lambda_rad))
    cosalpha = np.cos(lambda_rad) * np.cos(beta_rad) / np.cos(delta)
    sinalpha = (-np.sin(beta_rad) * np.sin(obliq) +
                np.cos(beta_rad) * np.cos(obliq) * np.sin(lambda_rad)) / np.cos(delta)
    alpha = np.arctan2(sinalpha, cosalpha)
    alpha = _wrap_to_2pi(alpha)

    return alpha, delta


def _tenv(dd, mm, ss):
    sgn, dd_mag = dd / dd, np.abs(dd)
    return sgn * (dd_mag + np.abs(mm) / 60.0 + np.abs(ss) / 3600.0)


def skyvec2ins(ra, dec,
               pa1, pa2, pa3,
               separation_as1, separation_as2, separation_as3,
               aper, start_date,
               npoints=360, nrolls=15, maxvroll=7.0):
    """
    Parameters
    ----------
    ra : float
        right ascension of target in decimal degrees (0-360)
    dec : float
        declination of target in decimal degrees (-90, 90)
    pa1, pa2, pa3 : float
        position angles of companions in degrees east of north
    separation_as1, separation_as2, separation_as3 : float
        separations of companions in arcseconds
    aper : jwxml.Aperture object
        Aperture as loaded from the instrument SIAF
    start_date : datetime.datetime
        Start date of the year-long interval evaluated by skyvec2ins
    npoints : int
        number of points to sample in the year-long interval
        to find observable dates (default: 360)
    nrolls : int
        number of roll angles in the allowed roll angle range to
        sample at each date (default: 15)
    maxvroll : float
        maximum number of degrees positive or negative roll around
        the boresight to allow (as designed: 7.0)

    .. note::

        `lambda_rad0` is the longitude of quadrature at
        day 0 of the code, so it should be 90 deg W of the
        solar longitude.

    Returns
    -------
    x : numpy.ndarray
        float array of length `npoints` containing days from starting
        date
    observable : numpy.ndarray
        uint8 array of shape (`nrolls`, `npoints`) that is 1 where
        the target is observable and 0 otherwise
    elongation_rad : numpy.ndarray
        float array of length `npoints` containing elongation of the
        observatory in radians
    roll_rad : numpy.ndarray
        float array of shape (`nrolls`, `npoints`) containing V3 PA
        in radians
    c1_x, c1_y, c2_x, c2_y, c3_x, c3_y : numpy.ndarray
        float array of shape (`nrolls`, `npoints`) containing the
        location of the companions in "Idl" (ideal) frame coordinates
    n_x, n_y, e_x, e_y : numpy.ndarray
        float array of shape (`nrolls`, `npoints`) containing the location
        of a reference "north" vector and "east" vector from the
        center in "Idl" (ideal) frame coordinates
    """
    # Per Chris Stark:
    # > lambda_rad0 is commented as the longitude of quadrature at day 0 of the code.
    # > So it should be 90 deg W of the solar longitude.
    # West is negative, so subtract 90 from the angle (in deg) and convert to radians.
    lambda_sun = sun_ecliptic_longitude(start_date)
    lambda_rad0 = np.deg2rad(lambda_sun - 90)

    # Conversions
    ra_rad = np.deg2rad(ra)
    dec_rad = np.deg2rad(dec)

    # We want to know the PA of the V3 axis. A simple way to do this
    # would be to form the 3 telescope axes as normal unit vectors in Cartesian
    # coordinates (v1, v2, v3) w/ v1 pointed at RA = 0, Dec = 0.
    # Then rotate those three by the angles required to point v1 at the star.
    # Then take the dot product of the v3 and v2 unit vectors with a N unit
    # vector to get orientation & PA.  I essentially do this below, but go
    # about the first step in a much more convoluted manner.

    # Now, let's approximate the pointing to get the v2,v3 celestial coordinates
    # b/c the actual pointing direction depends on the roll angle, which depends
    # on the pointing vector, etc...to solve the problem correctly we'd need
    # a root finder.  To speed things up, we simply assume the PA of the
    # telescope when pointing at the star is the same as when the
    # coronagraphmask is on the star--this is an approximation! The
    # approximation is valid if the coronagraphs are close to
    # the optical axis, which isn't too bad of an assumption.
    pointing_rad = [ra_rad, dec_rad]  # this is our approximation

    # Convert pointing from equatorial coords to alternative ecliptic coords
    # (alpha, delta) -> (lambda, beta)
    lambda_rad, beta_rad = ad2lb(pointing_rad[0], pointing_rad[1])  # first, convert to ecliptic coords

    # the ecliptic latitude of the direction of quadrature vs. time in steps of 1 degree
    quadrature_lambda_rad = np.arange(npoints) / npoints * 2 * np.pi + lambda_rad0
    j = np.where(quadrature_lambda_rad < 0)
    if len(j[0]) > 0:
        quadrature_lambda_rad[j] += 2 * np.pi

    j = np.where(quadrature_lambda_rad > 2 * np.pi)
    if len(j[0]) > 0:
        quadrature_lambda_rad[j] -= 2 * np.pi

    lmlsun_rad = (np.pi / 2) + (quadrature_lambda_rad - lambda_rad)  # the longitude relative to the sun
    elongation_rad, inc_rad = lb2ei(lmlsun_rad, beta_rad)  # now convert to alternative ecliptic coords
    assert elongation_rad.shape == (npoints,)
    assert inc_rad.shape == (npoints,)

    # Calculate celestial coordinates of V2 & V3 axis
    # First, calculate solar elongation & inclination

    v3_elongation_rad = elongation_rad + (np.pi / 2)
    v3_inc_rad = inc_rad.copy()  # explicit copy -- does mutating this affect things later??
    j = np.where(v3_elongation_rad > np.pi)  # Make sure the solar elongation is between 0 - 180 degrees
    if len(j[0]) > 0:
        # make sure all our angles follow standard definitions
        v3_elongation_rad[j] = (2 * np.pi) - v3_elongation_rad[j]
        v3_inc_rad[j] = np.pi + inc_rad[j]
        v3_inc_rad %= 2 * np.pi

    # Now transform to ecliptic coordinates (lambda, beta)
    v3_lmlsun_rad, v3_beta_rad = ei2lb(v3_elongation_rad, v3_inc_rad)
    # Finally, get the celestial coordinates (alpha, beta), i.e. (ra, dec)
    v3_lambda_rad = (np.pi / 2) + quadrature_lambda_rad - v3_lmlsun_rad
    v3_lambda_rad = _wrap_to_2pi(v3_lambda_rad)
    v3_alpha_rad, v3_delta_rad = lb2ad(v3_lambda_rad, v3_beta_rad)

    # At this point, we have the approximate celestial coordinates of the v2, v3 axes
    # Now we want to know the V3 PA. To get that, we create unit vectors and
    # eventually take a dot product...
    # We're simply taking spherical coordinates and making a cartesian unit vector
    # First, the V1 vector (the pointing vector)...
    alpha_rad = ra_rad
    delta_rad = dec_rad

    # a single unit vector
    unit_vec = np.array([
        np.cos(delta_rad) * np.cos(alpha_rad),
        np.cos(delta_rad) * np.sin(alpha_rad),
        np.sin(delta_rad)
    ])

    ux, uy, uz = unit_vec
    # alias for clarity below:
    v1_unit_vec = unit_vec

    # Here are the V3 unit vectors for each elongation
    v3_unit_vec = np.array([
        np.cos(v3_delta_rad) * np.cos(v3_alpha_rad),
        np.cos(v3_delta_rad) * np.sin(v3_alpha_rad),
        np.sin(v3_delta_rad)
    ])
    # Take the cross product to get v2 (v2 = v3 x v1)
    v2_unit_vec = np.zeros((3, len(v3_delta_rad)))
    for i in range(len(v3_delta_rad)):
        v2_unit_vec[:, i] = np.cross(v3_unit_vec[:, i], v1_unit_vec)

    # Now make unit vector arrays including all vehicle roll angles
    # Make the rotation matrix about the v1 axis
    vroll = np.linspace(-maxvroll, maxvroll, nrolls)
    cosvroll = np.cos(np.deg2rad(vroll))
    sinvroll = np.sin(np.deg2rad(vroll))

    # Find the v1 axis in cartesian space
    utensu = np.array([
        [ux * ux, uy * ux, uz * ux],
        [ux * uy, uy * uy, uz * uy],
        [ux * uz, uy * uz, uz * uz]
    ])
    ucpm = np.array([
        [0, -uz, uy],
        [uz, 0, -ux],
        [-uy, ux, 0]
    ])
    ident = np.array([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
    ])

    v2_uva = np.zeros((3, nrolls, npoints))
    v3_uva = np.zeros((3, nrolls, npoints))
    for i in range(npoints):
        for j in range(nrolls):
            # a rotation matrix about the v1 axis...
            rotation_matrix = cosvroll[j] * ident + sinvroll[j] * ucpm + (1 - cosvroll[j]) * utensu
            # rotate those unit vectors
            # (dotting rotation_matrix . v1_unit_vec should reproduce
            # v1_unit_vec, so don't bother)
            v2_uva[:, j, i] = np.dot(rotation_matrix, v2_unit_vec[:, i])
            v3_uva[:, j, i] = np.dot(rotation_matrix, v3_unit_vec[:, i])

    # Now that we have the unit vectors at all vehicle roll angles and elongations,
    # we use the dot product method of finding the PA
    # Let's add a tiny displacement in the north and east directions and
    # make unit vectors out of them...
    ddelta_rad = 1. / (1000. * 60. * 60. * 180. / np.pi)  # 1 milliarcsec in rad
    dalpha_rad = 1. / (1000. * 60. * 60. * 180. / np.pi)
    # First, get vectors from the origin to the tiny displacements
    #  ... a vector in the E direction
    unit_vec2 = np.array([
        np.cos(delta_rad + ddelta_rad) * np.cos(alpha_rad),
        np.cos(delta_rad + ddelta_rad) * np.sin(alpha_rad),
        np.sin(delta_rad + ddelta_rad)
    ])
    #  ... a vector in the N direction
    unit_vec3 = np.array([
        np.cos(delta_rad) * np.cos(alpha_rad + dalpha_rad),
        np.cos(delta_rad) * np.sin(alpha_rad + dalpha_rad),
        np.sin(delta_rad)
    ])
    # Now subtract off the unit vector pointing to the star to determine
    # the direction of the N vector at the star
    north_vec = unit_vec2 - unit_vec  # a single unit vector
    north_vec /= np.sqrt(np.sum(north_vec * north_vec))  # normalize it to make it a unit vector
    assert north_vec.shape == (3,)
    north_vec_all_rolls = np.zeros((3, nrolls, npoints))  # make it a 3 x nrolls x npoints
    # TODO refactor this to remove loops
    for i in range(npoints):
        for j in range(nrolls):
            # i, j transposed from IDL bc of indexing
            # just for making this the right shape to go with npoints and nrolls
            north_vec_all_rolls[:, j, i] = north_vec
    # Subtract off the unit vector pointing to the star to determine the direction of the E vector at the star
    east_vec = unit_vec3 - unit_vec
    east_vec /= np.sqrt(np.sum(east_vec * east_vec))  # make it a unit vector
    east_vec_all_rolls = np.zeros((3, nrolls, npoints))  # make it a npoints x nrolls x 3 array
    for i in range(npoints):
        for j in range(nrolls):
            east_vec_all_rolls[:, j, i] = east_vec

    # Take the dot products of the V3 vectors with the N unit vectors

    north_theta = np.arccos(np.sum(v3_uva * north_vec_all_rolls, axis=0))
    # DIY dot product, both are unit vectors so acos(sum(elemwise mult, 3)) gives angle in rad...
    # in the IDL, the first dim in this `shape` is the 3 elements of the vec
    # so we sum over axis 0
    assert north_theta.shape == (nrolls, npoints)
    # Take the dot products of the V3 vectors with the E unit vectors
    east_theta = np.arccos(np.sum(v3_uva * east_vec_all_rolls, axis=0))
    roll_rad = north_theta  # this is the PA of the V3 axis, not a telescope roll !!!!!
    j = np.where(east_theta > np.pi / 2)
    if len(j[0]) > 0:
        roll_rad[j] = 2 * np.pi - roll_rad[j]
    # We're not quite done...we need to know the orientation...
    # Now translate the v2 vector into ecliptic coords to determine quadrant of PA
    # ~~~ this might have to do with when JWST has to flip to reach things?
    v2_delta_rad = np.arcsin(v2_uva[:, :, 2])
    v2_sinalpha = v2_uva[:, :, 1] / np.cos(v2_delta_rad)
    v2_cosalpha = v2_uva[:, :, 0] / np.cos(v2_delta_rad)
    v2_alpha_rad = np.arctan2(v2_sinalpha, v2_cosalpha)
    j = np.where(v2_alpha_rad < 0)
    if len(j[0]) > 0:
        v2_alpha_rad[j] += 2 * np.pi

    # Finally, we have the approximate PA of the V3 axis!
    # ---------------------

    # ---------------------
    # Now that we have the approximate PA, we can calculate the
    # approximate (V2,V3) coordinates of the target and PA reference points.
    # The target & PA reference points are specified relative to North.
    # So do this, we rotate the points by the the V3 PA, which is measured
    # relative to North. Then we shift them to the (V2,V3) location of the coronagraph center.

    # Determine when the target is observable
    # we're going to explore other roll angles now, so
    # these are going to be 2D arrays (roll angle, elongation)
    vroll_rad = np.zeros((nrolls, npoints))
    vroll_rad[:, :] = np.deg2rad(vroll)[:, np.newaxis]
    # vroll_rad is an array: nrolls x npoints
    elongation_rad_arr = np.zeros((nrolls, npoints))
    elongation_rad_arr[:, :] = elongation_rad[np.newaxis, :]
    elongation_rad = elongation_rad_arr
    # elongation_rad is now an array of npoints x nrolls
    vpitch_rad = np.pi / 2 - elongation_rad  # definition change elong -> pitch
    # JWST pointing restrictions are in terms of solar pitch and solar roll
    sroll_rad = np.arcsin(np.sin(vroll_rad) * np.cos(vpitch_rad))
    assert sroll_rad.shape == (nrolls, npoints)
    spitch_rad = np.arctan(np.tan(vpitch_rad) / np.cos(vroll_rad))
    assert spitch_rad.shape == (nrolls, npoints)
    sroll = np.rad2deg(sroll_rad)
    spitch = np.rad2deg(spitch_rad)

    observable = np.zeros((nrolls, npoints), dtype=bool)

    # sun pitch constrained to -45 to +5.2
    mask = np.where(
        (spitch < 2.5) &
        (spitch >= -45.) &
        (np.abs(sroll) <= 5.2)
    )
    observable[mask] = True
    # personal communication Wayne Kinzel to Stark: (stark double checked with OBS mission req)
    # there is a linearized region
    mask2 = np.where(
        (spitch >= 2.5) &
        (spitch <= 5.2) &
        (np.abs(sroll) <= ((3.5 - 5.2) / (5.2 - 2.5) * (spitch - 2.5) + 5.2))
    )
    observable[mask2] = True

    x = np.arange(npoints) * (365.25 / npoints)  # days since launch

    pa_sep_pairs = [
        (pa1, separation_as1),
        (pa2, separation_as2),
        (pa3, separation_as3),
        (0., 0.1),  # 0.1 arcsec at North for North vector overlay
        (90., 0.1),  # 0.1 arcsec at East for East vector overlay
    ]
    results = []
    for pa, sep in pa_sep_pairs:
        idl_x, idl_y = detector_transform(
            nrolls, npoints, roll_rad,
            pa,
            sep,
            aper
        )
        results.extend((idl_x, idl_y))
    return [x, observable.astype(np.uint8), elongation_rad, roll_rad] + results


def detector_transform(nrolls, npoints, roll_rad, pa, separation_as, aper):
    pa_rad = np.deg2rad(pa)  # companion 1
    # Calculate the (V2,V3) coordinates of the coronagraph center
    # That's where we want to stick the target
    # The centers of the coronagraphic masks correspond to the XDetRef & YDetRef
    # locations on the detector (according to Colin Cox)
    cortelcoords = aper.Det2Tel(aper.XDetRef, aper.YDetRef)
    # convert arcseconds to radians
    cortelcoords_rad = np.asarray(cortelcoords) / 206264.806247
    # At this point, we have the coronagraph mask location in
    # science, telescope, and celestial coordinates

    # BEGIN DETECTOR COORDS TRANSFORM SECTION

    # Calculate approximate celestial ([alpha, delta], i.e. [ra, dec]) coordinates of
    # companions relative to stars. This can be added to [ra,dec] of
    # stars to get absolute positions (approximately)
    dcoords1_rad = (separation_as / 206264.806247) * np.array([np.sin(pa_rad), np.cos(pa_rad)])

    # First, the rotation...
    refceloffset1_rad = dcoords1_rad  # PA reference point
    refteloffset1_rad = np.zeros((nrolls, npoints, 2))
    refteloffset1_rad[:, :, 0] = np.cos(roll_rad) * refceloffset1_rad[0] - np.sin(roll_rad) * refceloffset1_rad[1]
    refteloffset1_rad[:, :, 1] = np.sin(roll_rad) * refceloffset1_rad[0] + np.cos(roll_rad) * refceloffset1_rad[1]

    reftelcoords1_rad = np.zeros((nrolls, npoints, 2))
    reftelcoords1_rad[:, :, 0] = refteloffset1_rad[:, :, 0] + cortelcoords_rad[0]
    reftelcoords1_rad[:, :, 1] = refteloffset1_rad[:, :, 1] + cortelcoords_rad[1]

    # Transform from (v2,v3) telescope coordinates to science coordinates
    refidlcoords1 = np.zeros((nrolls, npoints, 2))

    for i in range(npoints):
        for j in range(nrolls):
            temptelcoords1 = np.array([
                reftelcoords1_rad[j, i, 0],
                reftelcoords1_rad[j, i, 1]
            ]) * 206264.806247  # arcseconds
            tempidlcoords1 = aper.Tel2Idl(temptelcoords1[0], temptelcoords1[1])
            refidlcoords1[j, i] = tempidlcoords1

    # Detector coordinates of star and companion
    c1_x = refidlcoords1[:, :, 0]
    c1_y = refidlcoords1[:, :, 1]

    # END OF DETECTOR POS SECTION
    return c1_x, c1_y

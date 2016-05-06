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
from os.path import join, exists, dirname
import numpy as np
from .jwxml import SIAF

def _wrap_to_2pi(scalar_or_arr):
    """Offsets angles outside 0 <= x <= 2 * pi to lie within the interval"""
    return np.asarray(scalar_or_arr) % (2 * np.pi)

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
    sinlambda = (np.sin(delta_rad) * np.sin(obliq) + np.cos(delta_rad) * np.cos(obliq) * np.sin(alpha_rad)) / np.cos(beta_rad)
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

    obliq = _tenv(23, 26, 21.45)     # J2000 obliquity of Earth in degrees
    obliq *= np.pi / 180.

    delta = np.arcsin(np.sin(beta_rad) * np.cos(obliq) + np.cos(beta_rad) * np.sin(obliq) * np.sin(lambda_rad))
    cosalpha = np.cos(lambda_rad) * np.cos(beta_rad) / np.cos(delta)
    sinalpha = (-np.sin(beta_rad) * np.sin(obliq) + np.cos(beta_rad) * np.cos(obliq) * np.sin(lambda_rad)) / np.cos(delta)
    alpha = np.arctan2(sinalpha, cosalpha)
    alpha = _wrap_to_2pi(alpha)

    return alpha, delta

def _tenv(dd, mm, ss):
    sgn, dd_mag = dd / dd, np.abs(dd)
    return sgn * (dd_mag + np.abs(mm) / 60.0 + np.abs(ss) / 3600.0)

def skyvec2ins(ra, dec, pa1, pa2, pa3, separation_as1, separation_as2,
               separation_as3, instrname, apername, lambda_rad0, npoints=360, nrolls=14, maxvroll=7.0):
    """
    ;---------------------
    ; INPUTS
    ;---------------------
    ;ra = right ascension of target in decimal degrees
    ;dec = declination of target in decimal degrees
    ;pa1 - pa3 = position angle of companions 1 - 3 in degrees
    ;separation_as1 - separation_as3 = separation of companions 1 - 3 in arcseconds
    ;instrname = name of instrument
    ;apername = name of SIAF aperture
    lambda_rad0 = ecliptic longitude of quadrature with the sun, in radians, at the beginning of the year-long interval sampled (indirectly, this specifies the start date)
    (via Chris Stark: "lambda_rad0 is commented as the longitude of quadrature at day 0 of the code.  So it should be 90 deg W of the solar longitude")
    ;npoints = # of elongations calculated, default 360
    ;nrolls = # of roll angles calculated, default 20
    maxvroll = degrees max V roll to consider (not max allowable), default 7.0

    ;---------------------
    ; OUTPUTS
    ;---------------------
    ;x = 1D vector, days
    ;observable = 2D array (0 = point unobservable, 1 = observable)
    ;elongation_rad = 1D vector, elongation in radians
    ;roll_rad = 2D array, v3 PA in radians
    ;s_x, s_y = 2D array, stellar location on detector in pix
    ;c#_x, c#_y = 2D array, companion locations on detector in pix

    ;scisize = dimension of science aperture in pix
    ;sciscale = scale of science aperture in arcsec / pix
    ;n_x, n_y = north vector in science frame
    ;e_x, e_y = east vector in science frame
    """

    # Constants
    hours2deg = 360. / 24.
    deg2rad = np.pi / 180.
    obliq = _tenv(23, 26, 21.45)  # J2000 obliquity of Earth in degrees

    # Conversions
    ra_rad = ra * deg2rad
    dec_rad = dec * deg2rad
    obliq_rad = obliq * deg2rad
    paN_rad = 0.  # companion marking North
    paE_rad = 90. * deg2rad  # companion marking East
    pa1_rad = pa1 * deg2rad  # companion 1
    pa2_rad = pa2 * deg2rad  # companion 2
    pa3_rad = pa3 * deg2rad  # companion 3

    # Calculations
    # -- Given the desired instrument, load its parameters
    siaf_path = join(dirname(__file__), 'data', '{}_SIAF.xml'.format(instrname))
    assert exists(siaf_path), 'no SIAF for {}'.format(instrname)
    siaf = SIAF(instr=instrname, filename=siaf_path)
    aper = siaf[apername]

    scisize = [aper.XSciSize, aper.YSciSize]  # save this info to an output
    sciscale = [aper.XSciScale, aper.YSciScale]  # save this info to an output
    sciyangle = aper.V3SciYAngle

    # Calculate the (V2,V3) coordinates of the coronagraph center
    # That's where we want to stick the target
    # The centers of the coronagraphic masks correspond to the XDetRef & YDetRef
    # locations on the detector (according to Colin Cox)
    cordetcoords = np.array([aper.XDetRef, aper.YDetRef])
    corscicoords = aper.Det2Sci(cordetcoords[0], cordetcoords[1])
    coridlcoords = aper.Sci2Idl(corscicoords[0], corscicoords[1])
    cortelcoords = aper.Idl2Tel(coridlcoords[0], coridlcoords[1])
    # convert arcseconds to radians
    cortelcoords_rad = np.asarray(cortelcoords) / 206264.806247
    # At this point, we have the coronagraph mask location in
    # science, telescope, and celestial coordinates

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

    quadrature_lambda_rad = np.arange(npoints) / npoints * 2 * np.pi + lambda_rad0  # the ecliptic latitude of the direction of quadrature vs. time in steps of 1 degree
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

    # DIRECT TRANSLATION FROM IDL
    # v3_elongation_rad = elongation_rad + (np.pi/2)
    # v3_inc_rad = inc_rad
    # j = np.where(v3_elongation_rad > np.pi)  # Make sure the solar elongation is between 0 - 180 degrees
    # if len(j[0]) > 0:
    #    v3_elongation_rad[j] = (2 * np.pi) - v3_elongation_rad[j]  # make sure all our angles follow standard definitions
    #    v3_inc_rad[j] = np.pi + inc_rad[j]
    #    j = np.where(v3_inc_rad < 0)
    #    if len(j) > 0:
    #        v3_inc_rad[j] += 2 * np.pi
    #    j = np.where(v3_inc_rad > 2 * np.pi)
    #    if len(j) > 0:
    #        v3_inc_rad[j] -= 2 * np.pi

    # SHORTER VERSION THAT HANDLES CASE WHERE v3_inc_rad == 2pi DIFFERENTLY
    v3_elongation_rad = elongation_rad + (np.pi/2)
    v3_inc_rad = inc_rad.copy()  # explicit copy -- does mutating this affect things later??
    j = np.where(v3_elongation_rad > np.pi)  # Make sure the solar elongation is between 0 - 180 degrees
    if len(j[0]) > 0:
        v3_elongation_rad[j] = (2 * np.pi) - v3_elongation_rad[j]  # make sure all our angles follow standard definitions
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
        [np.cos(delta_rad) * np.cos(alpha_rad)],
        [np.cos(delta_rad) * np.sin(alpha_rad)],
        [np.sin(delta_rad)]
    ])
    unit_vec = np.array([
        np.cos(delta_rad) * np.cos(alpha_rad),
        np.cos(delta_rad) * np.sin(alpha_rad),
        np.sin(delta_rad)
    ])

    ux = unit_vec[0]
    uy = unit_vec[1]
    uz = unit_vec[2]

    v1_unit_vec = np.zeros((3, npoints))
    # a bunch of v1 unit vectors, all the same, for each elongation
    v1_unit_vec[:, :] = unit_vec[:, np.newaxis]
    assert v1_unit_vec.shape == (3, 360)
    # Here are the V3 unit vectors for each elongation
    v3_unit_vec = np.array([
        np.cos(v3_delta_rad) * np.cos(v3_alpha_rad),
        np.cos(v3_delta_rad) * np.sin(v3_alpha_rad),
        np.sin(v3_delta_rad)
    ])
    assert v1_unit_vec.shape == (3, 360)
    # Take the cross product to get v2 (v2 = v3 x v1)
    v2_unit_vec = np.zeros((3, len(v3_delta_rad)))
    for i in range(len(v3_delta_rad)):
        v2_unit_vec[:, i] = np.cross(v3_unit_vec[:,i], v1_unit_vec[:,i])
    assert v1_unit_vec.shape == (3, 360)

    # Now make unit vector arrays including all vehicle roll angles
    # Make the rotation matrix about the v1 axis
    vroll = (2 * np.arange(nrolls, dtype=np.float64) / (nrolls - 1) - 1) * maxvroll
    cosvroll = np.cos(vroll * deg2rad)
    sinvroll = np.sin(vroll * deg2rad)

    # Find the v1 axis in cartesian space
    utensu = np.array([
        [ux * ux, uy * ux, uz * ux],
        [ux * uy, uy * uy, uz * uy],
        [ux * uz, uy * uz, uz * uz]
    ])
    ucpm = np.array([
        [  0, -uz,  uy],
        [ uz,   0, -ux],
        [-uy,  ux,   0]
    ])
    ident = np.array([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
    ])
    v1_uva = np.zeros((3, nrolls, npoints))
    v2_uva = np.zeros((3, nrolls, npoints))
    v3_uva = np.zeros((3, nrolls, npoints))
    for i in range(npoints):
        for j in range(nrolls):
            # a rotation matrix about the v1 axis...
            R = cosvroll[j] * ident + sinvroll[j] * ucpm + (1-cosvroll[j]) * utensu
            # rotate those unit vectors
            v1_uva[:, j, i] = np.dot(R, v1_unit_vec[:, i])
            assert np.allclose(np.dot(R, v1_unit_vec[:, i]), v1_unit_vec[:, i])
            v2_uva[:, j, i] = np.dot(R, v2_unit_vec[:, i])
            v3_uva[:, j, i] = np.dot(R, v3_unit_vec[:, i])


    # Now that we have the unit vectors at all vehicle roll angles and elongations,
    # we use the dot product method of finding the PA
    # Let's add a tiny displacement in the north and east directions and
    # make unit vectors out of them...
    ddelta_rad = 1. / (1000. * 60. * 60. * 180. / np.pi) # 1 milliarcsec in rad
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
    Nvec = unit_vec2 - unit_vec  # a single unit vector
    Nvec /= np.sqrt(np.sum(Nvec * Nvec))  # normalize it to make it a unit vector
    assert Nvec.shape == (3,)
    Nvec_a = np.zeros((3, nrolls, npoints))  # make it a 3 x nrolls x npoints
    # TODO refactor this to remove loops
    for i in range(npoints):
        for j in range(nrolls):
            # j, j transposed from idl bc of indexing
            Nvec_a[:, j, i] = Nvec  # just for making this the right shape to go with npoints and nrolls
    # Subtract off the unit vector pointing to the star to determine the direction of the E vector at the star
    Evec = unit_vec3 - unit_vec
    Evec /= np.sqrt(np.sum(Evec * Evec))  # make it a unit vector
    Evec_a = np.zeros((3, nrolls, npoints))  # make it a npoints x nrolls x 3 array
    for i in range(npoints):
        for j in range(nrolls):
            Evec_a[:, j, i] = Evec

    # Take the dot products of the V3 vectors with the N unit vectors

    Ntheta = np.arccos(np.sum(v3_uva * Nvec_a, axis=0))
    # DIY dot product, both are unit vectors so acos(sum(elemwise mult, 3)) gives angle in rad...
    # in the IDL, the first dim in this `shape` is the 3 elements of the vec
    # so we sum over axis 0
    assert Ntheta.shape == (nrolls, npoints)
    # Take the dot products of the V3 vectors with the E unit vectors
    Etheta = np.arccos(np.sum(v3_uva * Evec_a, axis=0))
    roll_rad = Ntheta  # this is the PA of the V3 axis, not a telescope roll !!!!!
    j = np.where(Etheta > np.pi / 2)
    if len(j[0]) > 0:
        roll_rad[j] = 2 * np.pi - roll_rad[j]
    # We're not quite done...we need to know the orientation...
    # Now translate the v2 vector into ecliptic coords to determine quadrant of PA
    # ~~~ this might have to do with when JWST has to flip to reach things?
    v2_delta_rad = np.arcsin(v2_uva[:,:,2])
    v2_sinalpha = v2_uva[:,:,1] / np.cos(v2_delta_rad)
    v2_cosalpha = v2_uva[:,:,0] / np.cos(v2_delta_rad)
    v2_alpha_rad = np.arctan2(v2_sinalpha, v2_cosalpha)
    j = np.where(v2_alpha_rad < 0)
    if len(j[0]) > 0:
        v2_alpha_rad[j] += 2 * np.pi
    v2_beta_rad = np.arcsin(np.sin(v2_delta_rad) * np.cos(obliq_rad) - np.cos(v2_delta_rad) * np.sin(obliq_rad) * np.sin(v2_alpha_rad))

    # Finally, we have the approximate PA of the V3 axis!
    # ---------------------

    # ---------------------
    # Now that we have the approximate PA, we can calculate the
    # approximate (V2,V3) coordinates of the target and PA reference points.
    # The target & PA reference points are specified relative to North.
    # So do this, we rotate the points by the the V3 PA, which is measured
    # relative to North. Then we shift them to the (V2,V3) location of the coronagraph center.
    #
    # BEGIN DETECTOR COORDS TRANSFORM SECTION

    # Calculate approximate celestial ([alpha, delta], i.e. [ra, dec]) coordinates of
    # companions relative to stars. This can be added to [ra,dec] of
    # stars to get absolute positions (approximately)

    # tiny offset in north and east so they can be rotated identically to companions for plotting
    dcoordsN_rad = (0.1 / 206264.806247) * np.array([np.sin(paN_rad), np.cos(paN_rad)])
    dcoordsE_rad = (0.1 / 206264.806247) * np.array([np.sin(paE_rad), np.cos(paE_rad)])
    # dcoords 1-3 are coords of companions
    dcoords1_rad = (separation_as1 / 206264.806247) * np.array([np.sin(pa1_rad), np.cos(pa1_rad)])
    dcoords2_rad = (separation_as2 / 206264.806247) * np.array([np.sin(pa2_rad), np.cos(pa2_rad)])
    dcoords3_rad = (separation_as3 / 206264.806247) * np.array([np.sin(pa3_rad), np.cos(pa3_rad)])

    # First, the rotation...
    targceloffset_rad = np.array([0, 0])  # assume telescope is pointed at star
    targteloffset_rad = np.zeros((nrolls, npoints, 2))
    targteloffset_rad[:, :, 0] = np.cos(roll_rad) * targceloffset_rad[0] - np.sin(roll_rad) * targceloffset_rad[1]
    targteloffset_rad[:, :, 1] = np.sin(roll_rad) * targceloffset_rad[0] + np.cos(roll_rad) * targceloffset_rad[1]
    refceloffsetN_rad = dcoordsN_rad  # PA reference point
    refceloffsetE_rad = dcoordsE_rad  # PA reference point
    refceloffset1_rad = dcoords1_rad  # PA reference point
    refceloffset2_rad = dcoords2_rad  # PA reference point
    refceloffset3_rad = dcoords3_rad  # PA reference point
    refteloffsetN_rad = np.zeros((nrolls, npoints, 2))
    refteloffsetE_rad = np.zeros((nrolls, npoints, 2))
    refteloffset1_rad = np.zeros((nrolls, npoints, 2))
    refteloffset2_rad = np.zeros((nrolls, npoints, 2))
    refteloffset3_rad = np.zeros((nrolls, npoints, 2))
    refteloffsetN_rad[:, :, 0] = np.cos(roll_rad) * refceloffsetN_rad[0] - np.sin(roll_rad) * refceloffsetN_rad[1]
    refteloffsetE_rad[:, :, 0] = np.cos(roll_rad) * refceloffsetE_rad[0] - np.sin(roll_rad) * refceloffsetE_rad[1]
    refteloffset1_rad[:, :, 0] = np.cos(roll_rad) * refceloffset1_rad[0] - np.sin(roll_rad) * refceloffset1_rad[1]
    refteloffset2_rad[:, :, 0] = np.cos(roll_rad) * refceloffset2_rad[0] - np.sin(roll_rad) * refceloffset2_rad[1]
    refteloffset3_rad[:, :, 0] = np.cos(roll_rad) * refceloffset3_rad[0] - np.sin(roll_rad) * refceloffset3_rad[1]
    refteloffsetN_rad[:, :, 1] = np.sin(roll_rad) * refceloffsetN_rad[0] + np.cos(roll_rad) * refceloffsetN_rad[1]
    refteloffsetE_rad[:, :, 1] = np.sin(roll_rad) * refceloffsetE_rad[0] + np.cos(roll_rad) * refceloffsetE_rad[1]
    refteloffset1_rad[:, :, 1] = np.sin(roll_rad) * refceloffset1_rad[0] + np.cos(roll_rad) * refceloffset1_rad[1]
    refteloffset2_rad[:, :, 1] = np.sin(roll_rad) * refceloffset2_rad[0] + np.cos(roll_rad) * refceloffset2_rad[1]
    refteloffset3_rad[:, :, 1] = np.sin(roll_rad) * refceloffset3_rad[0] + np.cos(roll_rad) * refceloffset3_rad[1]
    # Now, the shift...
    targtelcoords_rad = np.zeros((nrolls, npoints, 2))
    targtelcoords_rad[:, :, 0] = targteloffset_rad[:, :, 0] + cortelcoords_rad[0]
    targtelcoords_rad[:, :, 1] = targteloffset_rad[:, :, 1] + cortelcoords_rad[1]
    reftelcoordsN_rad = np.zeros((nrolls, npoints, 2))
    reftelcoordsE_rad = np.zeros((nrolls, npoints, 2))
    reftelcoords1_rad = np.zeros((nrolls, npoints, 2))
    reftelcoords2_rad = np.zeros((nrolls, npoints, 2))
    reftelcoords3_rad = np.zeros((nrolls, npoints, 2))
    reftelcoordsN_rad[:, :, 0] = refteloffsetN_rad[:, :, 0] + cortelcoords_rad[0]
    reftelcoordsE_rad[:, :, 0] = refteloffsetE_rad[:, :, 0] + cortelcoords_rad[0]
    reftelcoords1_rad[:, :, 0] = refteloffset1_rad[:, :, 0] + cortelcoords_rad[0]
    reftelcoords2_rad[:, :, 0] = refteloffset2_rad[:, :, 0] + cortelcoords_rad[0]
    reftelcoords3_rad[:, :, 0] = refteloffset3_rad[:, :, 0] + cortelcoords_rad[0]
    reftelcoordsN_rad[:, :, 1] = refteloffsetN_rad[:, :, 1] + cortelcoords_rad[1]
    reftelcoordsE_rad[:, :, 1] = refteloffsetE_rad[:, :, 1] + cortelcoords_rad[1]
    reftelcoords1_rad[:, :, 1] = refteloffset1_rad[:, :, 1] + cortelcoords_rad[1]
    reftelcoords2_rad[:, :, 1] = refteloffset2_rad[:, :, 1] + cortelcoords_rad[1]
    reftelcoords3_rad[:, :, 1] = refteloffset3_rad[:, :, 1] + cortelcoords_rad[1]

    # Transform from (v2,v3) telescope coordinates to science coordinates
    targdetcoords = np.zeros((nrolls, npoints, 2))
    targscicoords = np.zeros((nrolls, npoints, 2))
    for i in range(npoints):
        for j in range(nrolls):
            temptelcoords = np.array([targtelcoords_rad[j, i, 0], targtelcoords_rad[j, i, 1]]) * 206264.806247  # arcseconds
            tempidlcoords = aper.Tel2Idl(temptelcoords[0], temptelcoords[1])
            tempscicoords = aper.Idl2Sci(tempidlcoords[0], tempidlcoords[1])
            targscicoords[j, i] = np.asarray(tempscicoords)
            tempdetcoords = aper.Sci2Det(tempscicoords[0], tempscicoords[1])
            targdetcoords[j, i] = np.asarray(tempdetcoords)


    refdetcoordsN = np.zeros((nrolls, npoints, 2))
    refdetcoordsE = np.zeros((nrolls, npoints, 2))
    refdetcoords1 = np.zeros((nrolls, npoints, 2))
    refdetcoords2 = np.zeros((nrolls, npoints, 2))
    refdetcoords3 = np.zeros((nrolls, npoints, 2))
    refscicoordsN = np.zeros((nrolls, npoints, 2))
    refscicoordsE = np.zeros((nrolls, npoints, 2))
    refscicoords1 = np.zeros((nrolls, npoints, 2))
    refscicoords2 = np.zeros((nrolls, npoints, 2))
    refscicoords3 = np.zeros((nrolls, npoints, 2))
    for i in range(npoints):
        for j in range(nrolls):
            temptelcoordsN = np.array([reftelcoordsN_rad[j, i, 0], reftelcoordsN_rad[j, i, 1]]) * 206264.806247  # arcseconds
            temptelcoordsE = np.array([reftelcoordsE_rad[j, i, 0], reftelcoordsE_rad[j, i, 1]]) * 206264.806247  # arcseconds
            temptelcoords1 = np.array([reftelcoords1_rad[j,i,0], reftelcoords1_rad[j,i,1]]) * 206264.806247  # arcseconds
            temptelcoords2 = np.array([reftelcoords2_rad[j,i,0], reftelcoords2_rad[j,i,1]]) * 206264.806247  # arcseconds
            temptelcoords3 = np.array([reftelcoords3_rad[j,i,0], reftelcoords3_rad[j,i,1]]) * 206264.806247  # arcseconds
            tempidlcoordsN = aper.Tel2Idl(temptelcoordsN[0], temptelcoordsN[1])
            tempidlcoordsE = aper.Tel2Idl(temptelcoordsE[0], temptelcoordsE[1])
            tempidlcoords1 = aper.Tel2Idl(temptelcoords1[0], temptelcoords1[1])
            tempidlcoords2 = aper.Tel2Idl(temptelcoords2[0], temptelcoords2[1])
            tempidlcoords3 = aper.Tel2Idl(temptelcoords3[0], temptelcoords3[1])
            tempscicoordsN = aper.Idl2Sci(tempidlcoordsN[0], tempidlcoordsN[1])
            tempscicoordsE = aper.Idl2Sci(tempidlcoordsE[0], tempidlcoordsE[1])
            tempscicoords1 = aper.Idl2Sci(tempidlcoords1[0], tempidlcoords1[1])
            tempscicoords2 = aper.Idl2Sci(tempidlcoords2[0], tempidlcoords2[1])
            tempscicoords3 = aper.Idl2Sci(tempidlcoords3[0], tempidlcoords3[1])
            refscicoordsN[j, i] = tempscicoordsN
            refscicoordsE[j, i] = tempscicoordsE
            refscicoords1[j, i] = tempscicoords1
            refscicoords2[j, i] = tempscicoords2
            refscicoords3[j, i] = tempscicoords3
            tempdetcoordsN = aper.Sci2Det(tempscicoordsN[0], tempscicoordsN[1])
            tempdetcoordsE = aper.Sci2Det(tempscicoordsE[0], tempscicoordsE[1])
            tempdetcoords1 = aper.Sci2Det(tempscicoords1[0], tempscicoords1[1])
            tempdetcoords2 = aper.Sci2Det(tempscicoords2[0], tempscicoords2[1])
            tempdetcoords3 = aper.Sci2Det(tempscicoords3[0], tempscicoords3[1])
            refdetcoordsN[j,i] = tempdetcoordsN
            refdetcoordsE[j,i] = tempdetcoordsE
            refdetcoords1[j,i] = tempdetcoords1
            refdetcoords2[j,i] = tempdetcoords2
            refdetcoords3[j,i] = tempdetcoords3

    # Detector coordinates of star and companion
    s_x = targscicoords[:, :, 0] - corscicoords[0]
    s_y = targscicoords[:, :, 1] - corscicoords[1]
    c1_x = refscicoords1[:, :, 0] - corscicoords[0]
    c1_y = refscicoords1[:, :, 1] - corscicoords[1]
    c2_x = refscicoords2[:, :, 0] - corscicoords[0]
    c2_y = refscicoords2[:, :, 1] - corscicoords[1]
    c3_x = refscicoords3[:, :, 0] - corscicoords[0]
    c3_y = refscicoords3[:, :, 1] - corscicoords[1]
    n_x = refscicoordsN[:, :, 0] - corscicoords[0]
    n_y = refscicoordsN[:, :, 1] - corscicoords[1]
    e_x = refscicoordsE[:, :, 0] - corscicoords[0]
    e_y = refscicoordsE[:, :, 1] - corscicoords[1]

    # END OF DETECTOR POS SECTION

    # Determine when the target is observable
    # we're going to explore other roll angles now, so
    # these are going to be 2D arrays (roll angle, elongation)
    vroll_rad = np.zeros((nrolls, npoints))
    vroll_rad[:, :] = vroll[:, np.newaxis] * deg2rad
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
    sroll = sroll_rad / deg2rad
    spitch = spitch_rad / deg2rad

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

    return (x, observable.astype(np.uint8), elongation_rad, roll_rad, s_x, s_y, c1_x, c1_y,
            c2_x, c2_y, c3_x, c3_y,
            scisize, sciscale, np.deg2rad(sciyangle), n_x, n_y, e_x, e_y)

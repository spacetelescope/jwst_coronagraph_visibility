JWST Target Visibility Calculator
=================================

**This initial release is for user testing and (as with any software) bugs may remain! Report any issues at https://github.com/spacetelescope/jwst_visibility/issues/new or via email to the authors**

*Authors: Christopher Stark (cstark@stsci.edu), Joseph Long (jlong@stsci.edu)*

Introduction
------------

The allowed pointing of JWST leads to target visibility that depends on ecliptic latitude, and the range of roll angles allowed depends on solar elongation. The allowed PAs for a target can thus be a complicated function of time. As a result, it can be difficult to 1) understand the possible orientations of a given target on the detector, 2) determine the ideal roll angle offsets for multi-roll observations, and 3) determine a group of targets that are simultaneously visible. The JWST Target Visibility Calculator (TVC) was created to address these issues and assist with creating APT programs and diagnosing scheduling errors.

We stress that the TVC is designed to provide quick illustrations of the possible observable orientations for a given target. As such, the TVC rapidly approximates JWST’s pointing restrictions and **does not query the official JWST Proposal Constraint Generator (PCG)**. The TVC does not include detailed pointing restrictions like Earth and Moon avoidance, etc. Additionally, results may differ from official constraints by a degree or so. **Users should treat the results as close approximations.**

Additionally, detector geometry (e.g. conversion from companion PA to detector pixel positions) is provided by the SIAF (Science Image Aperture File). The SIAF is a standardized format for manipulating instrument apertures and coordinate conversions, maintained by STScI as part of operating JWST. As of this writing (April 2016), the code includes its own copies of the `PRDDEVSOC-D-012` version of the NIRCam and MIRI SIAFs.

Using the GUI
-------------

After installation, the ``jwst-visibility-gui`` command will be available. Launch the GUI by executing this command.

In the GUI you will see a control panel on the left and a plot panel on the right. Within the control panel there are:

  * a SIMBAD Target Resolver frame
  * RA and Dec fields to supply coordinates manually
  * a Companions frame
  * an Instrument/Mask Selector frame
  * buttons to toggle between detector position angles and V3 (telescope) position angles in the observability plot
  * an Update Plot button
  * a few buttons to load example calculations

To find a target, type the target name into the SIMBAD Target Resolver text box and click Search. If SIMBAD is unable to find a match, the result “Target not found” will be displayed. Upon successful matching, the target’s SIMBAD ID, RA, and declination will be displayed. Users may also enter RA and Dec values as decimal degrees in the boxes provided.

Once a target has been resolved, you may click Update Plot to calculate the target’s visibility. The calculation usually takes less than a minute (as estimated on a 2015 Macbook Pro).

For additional documentation explaining how to interpret the plots produced, see TODO URL.

Using the Python API
--------------------

The current Python API is a direct translation of the original IDL Coronagraph Visibility Tool for JWST by Christopher Stark. The function in which the bulk of the calculation happens is ``jwst_visibility.skyvec2ins.skyvec2ins`` (in ``jwst_visibility/skyvec2ins.py``).

Inputs
^^^^^^

  * ``ra`` -- right ascension in decimal degrees (0 <= ``ra`` <= 360)
  * ``dec`` -- declination in decimal degrees (-90 <= ``dec`` <= 90)
  * ``pa1``, ``pa2``, ``pa3`` -- position angles of the three possible companions (can be 0.0)
  * ``separation_as1``, ``separation_as2``, ``separation_as3`` -- separations of the three possible companions (can be 0.0)
  * ``instrname`` -- string (one of 'NIRCam', 'NIRISS', 'NIRSpec', 'MIRI')
  * ``apername`` -- string (an aperture name as defined in the JWST SIAF -- the Science Instrument Aperture File)
  * ``npoints`` -- number of points to sample in a complete orbit (default: 360)
  * ``nrolls`` -- number of increments in (-10 degrees, +10 degrees) to evaluate JWST roll angles (default: 20)

Outputs
^^^^^^^

n.b. 2D arrays are of shape (nrolls, npoints) in the Python (rows, cols) indexing convention, though some quantities are the same for all roll angles.

  * ``x`` -- 1D array, days
  * ``observable`` -- 2D mask array (0 = point unobservable, 1 = observable)
  * ``elongation_rad`` -- 1D array, elongation in radians
  * ``roll_rad`` -- 2D array, v3 PA in radians
  * ``s_x``, ``s_y`` -- 2D arrays, stellar location on detector in pix
  * ``c1_x``, ``c1_y``, ``c2_x``, ``c2_y``, ``c3_x``, ``c3_y`` -- 2D arrays, companion locations on detector in pix
  * ``scisize`` -- tuple of x, y dimensions of science aperture in pix
  * ``sciscale`` -- tuple of x, y scale of science aperture in arcsec / pix
  * ``n_x``, ``n_y`` -- 2D arrrays of north vector coordinates in science frame
  * ``e_x``, ``e_y`` -- 2D arrays of east vector coordinates in science frame

Should there be user demand, a more conventional Python API may be developed, but the current plan is to support the GUI as the primary interface.

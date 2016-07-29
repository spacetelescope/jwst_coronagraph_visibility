JWST Target Visibility Calculator
=================================

**This initial release is for user testing and (as with any software) bugs may remain! Report any issues at https://github.com/spacetelescope/jwst_visibility/issues/new or via email to the authors**

*Authors: Christopher Stark (cstark@stsci.edu), Joseph Long (jlong@stsci.edu)*

.. image:: screenshot.png
   :width: 60%
   :align: center
   :alt: Screenshot of the JWST Target Visibility Calculator showing target HR 8799 with three companions plotted.

Introduction
------------

The allowed pointing of JWST leads to target visibility that depends on ecliptic latitude, and the range of roll angles allowed depends on solar elongation. The allowed PAs for a target can thus be a complicated function of time. As a result, it can be difficult to 1) understand the possible orientations of a given target on the detector, 2) determine the ideal roll angle offsets for multi-roll observations, and 3) determine a group of targets that are simultaneously visible. The JWST Target Visibility Calculator (TVC) was created to address these issues and assist with creating APT programs and diagnosing scheduling errors.

We stress that the TVC is designed to provide quick illustrations of the possible observable orientations for a given target. As such, the TVC rapidly approximates JWST’s pointing restrictions and **does not query the official JWST Proposal Constraint Generator (PCG)**. The TVC does not include detailed pointing restrictions like Earth and Moon avoidance, etc. Additionally, results may differ from official constraints by a degree or so. **Users should treat the results as close approximations.**

Additionally, detector geometry (e.g. conversion from sky coordinates to the instrument's ``Idl`` frame) is provided by the SIAF (Science Image Aperture File). The SIAF is a standardized format for manipulating instrument apertures and coordinate conversions, maintained by STScI as part of operating JWST. As of this writing (April 2016), the code includes its own copies of the PRDDEVSOC-D-012 version of the NIRCam and MIRI SIAFs.

Installation
------------

We have done our best to support a variety of Python installations, but different users may have particular customizations to their environment that we have not anticipated. Your mileage may vary!

For the most predictable installation experience, we recommend installing Python with the conda package manager and using the `AstroConda <http://astroconda.readthedocs.io/en/latest/index.html>`_ package set supported by STScI. After following the `AstroConda installation instructions <http://astroconda.readthedocs.io/en/latest/installation.html>`_, the following command will install the JWST Target Visibility Calculator:

   pip install -e git+https://github.com/spacetelescope/jwst_visibility.git#egg=jwst-visibility

If you would rather not use AstroConda, read on.

What are the prerequisites for installing the JWST Target Visibility Calculator?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The JWST Target Visibility Calculator depends on:

 * Python 2.7, 3.4, or 3.5
 * NumPy version 1.9.0 (or greater)
 * Matplotlib version 1.4.2 (or greater)
 * requests version 2.8.1 (or greater)

Running the automated test suite additionally depends on pytest version 2.9.1 (or greater).

Can I install the JWST Target Visibility Calculator without using conda at all?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The pip command should be sufficient to install the prerequisites in most cases. If you experience difficulties, try installing NumPy alone first, then using the given installation command.

Can I install the JWST Target Visibility Calculator without all of AstroConda?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Yes! If you're already using conda, but do not want to install everything else included in AstroConda, simply install the prerequisite packages with ``conda install`` before running the ``pip`` command listed above.

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

The plot shows the solar elongation required to observe the target as a black line, with the observable portions highlighted in red. For each red portion, the plot shows the range of allowed position angles in blue ("Aperture PA" denotes the angle from North to the y axis of the science frame on the detector in the eastward direction). To plot the PA of the V3 axis instead, click on the "V3 PA" button on the control panel and then Update Plot.

You can inspect any individually calculated blue point by clicking on it. You can also zoom in on any region of the plot using the `standard matplotlib controls <http://matplotlib.org/users/navigation_toolbar.html>`_ along the toolbar at the bottom of the plotting area.

To display the orientation of the target in the science frame on the detector, we must choose an astronomical reference direction. To do so, enable one of the Companions by clicking on the check box in the left column. Specify the companion’s PA (in degrees E of N) and separation (in arcseconds). We note that a companion can be thought of as a binary star, an exoplanet, the location of a disk’s major or minor axis, or any sort of reference applicable to the astrophysical scene of interest. One can add up to 3 companions.

Before we update the plot, select the instrument and coronagraphic mask that you’d like to use to observe the target using the drop-down menus. Finally, click Update Plot again to refresh the plots. When clicking on the plotted points, the crosshairs select the nearest blue point. The corresponding companion points are marked on the science frame panel. The North and East axes are also shown on the science frame panel as solid red and yellow lines, respectively.

When clicking on the science frame panel, the nearest companion point is selected and highlighted. The corresponding companion points are also highlighted, and the corresponding PA is shown with crosshairs in the left panel. Below the science frame panel, the separation (in pixels) and angle on the detector (CCW relative to +y axis) are displayed for each companion point.

By default, the science frame panel is drawn to show the full extent of the coronagraph aperture selected. The aperture is outlined in dashed red. (For MIRI coronagraphs, the clear aperture is outlined in solid red within the dashed red overall aperture.) For close companions, the zoom tool from the toolbar will let you draw a box around the region of interest. When in zoom mode, clicking companion or PA points will not have any effect. Click the zoom toolbar button again to exit zoom mode and restore the original click behavior.

Using the Python API
--------------------

The current Python API is a direct translation of the original IDL Coronagraph Visibility Tool for JWST by Christopher Stark. The function in which the bulk of the calculation happens is ``jwst_visibility.skyvec2ins.skyvec2ins`` (in ``jwst_visibility/skyvec2ins.py``).

Parameters
^^^^^^^^^^

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
lambda_rad0 : float
    ecliptic longitude of quadrature with the sun, in radians,
    at the beginning of the year-long interval sampled by
    this function (indirectly, this specifies the start date).
npoints : int
    number of points to sample in the year-long interval
    to find observable dates (default: 360)
nrolls : int
    number of roll angles in the allowed roll angle range to
    sample at each date (default: 14)
maxvroll : float
    maximum number of degrees positive or negative roll around
    the boresight to allow (as designed: 7.0)

Note: ``lambda_rad0`` is the longitude of quadrature at
day 0 of the code, so it should be 90 deg W of the
solar longitude.

Returns
^^^^^^^

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

Should there be user demand, a more conventional Python API may be developed, but the current plan is to support the GUI as the primary interface.

Automated software testing
--------------------------

To ensure the correctness of the Python to IDL translation, the ``jwst_visibility/tests/targets/`` folder contains CSV files output by the IDL tool for the ``skyvec2ins`` procedure. The Python code contains automated tests that run the Python ``skyvec2ins`` function and compare the output.

To run the test suite, install ``pytest`` and run the command::

    $ py.test

from this folder.

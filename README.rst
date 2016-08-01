JWST Target Visibility Tool
=================================

**This initial release is for user testing and (as with any software) bugs may remain! Report any issues at https://github.com/spacetelescope/jwst_visibility/issues/new or via email to the authors**

*Authors: Christopher Stark (cstark@stsci.edu), Joseph Long (jlong@stsci.edu)*

.. image:: screenshot.png
   :width: 60%
   :align: center
   :alt: Screenshot of the JWST Target Visibility Tool showing target HR 8799 with three companions plotted.

The allowed pointing of JWST leads to target visibility that depends on ecliptic latitude, and the range of roll angles allowed depends on solar elongation. The allowed PAs for a target can thus be a complicated function of time. As a result, it can be difficult to 1) understand the possible orientations of a given target on the detector, 2) determine the ideal roll angle offsets for multi-roll observations, and 3) determine a group of targets that are simultaneously visible. The JWST Target Visibility Tool (TVT) was created to address these issues and assist with creating APT programs and diagnosing scheduling errors.

We stress that the TVT is designed to provide quick illustrations of the possible observable orientations for a given target. As such, the TVT rapidly approximates JWST’s pointing restrictions and **does not query the official JWST Proposal Constraint Generator (PCG)**. The TVT does not include detailed pointing restrictions like Earth and Moon avoidance, etc. Additionally, results may differ from official constraints by a degree or so. **Users should treat the results as close approximations.**

Additionally, detector geometry (e.g. conversion from sky coordinates to the instrument's ``Idl`` frame) is provided by the SIAF (Science Image Aperture File). The SIAF is a standardized format for manipulating instrument apertures and coordinate conversions, maintained by STScI as part of operating JWST. As of this writing (April 2016), the code includes its own copies of the PRDDEVSOC-D-012 version of the NIRCam and MIRI SIAFs.

For installation instructions and usage instructions, see the `documentation <https://github.com/spacetelescope/jwst_visibility/blob/master/docs/index.rst>`_ on GitHub or ``docs/index.rst`` in this repository.

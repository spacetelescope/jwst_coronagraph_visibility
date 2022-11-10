# jwst_coronagraph_visibility: James Webb Space Telescope Coronagraph Visibility Tool

[![Current Release](https://img.shields.io/github/v/release/spacetelescope/jwst_coronagraph_visibility.svg)](https://github.com/spacetelescope/jwst_coronagraph_visibility/releases/latest/)
[![License](https://img.shields.io/github/license/spacetelescope/jwst_coronagraph_visibility)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.6-blue.svg)](https://www.python.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4488421.svg)](https://doi.org/10.5281/zenodo.4488421)
[![STScI](https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat)](http://www.stsci.edu)
**SIAF version:** PRDOPSSOC-027

The **J**ames **W**ebb **S**pace **T**elescope **Coronagraph Visibility** Tool (CVT) is a GUI-based target visibility 
tool for assessing target visibilities and available position angles versus time relative to 
MIRI<sup id="a1">[[1]](#f1)</sup> and NIRCam<sup id="a2">[[2]](#f1)</sup> coronagraphic masks. The CVT is available as 
a standalone gui-based python tool (AstroConda package) or as a macOS app bundle. 

The allowed pointing of JWST leads to target visibility that depends on ecliptic latitude, and the range of roll angles 
allowed depends on solar elongation. The allowed position angles for a target can thus be a complicated function of 
time. As a result, it can be difficult to: (1) understand the possible orientations of a given target on the detector, 
especially in relation to any instrumental obscurations; (2) determine the ideal roll angle offsets for multi-roll 
observations; and (3) determine a group of targets that are simultaneously visible. The CVT was created to address 
these issues and assist in planning MIRI and NIRCam programs prior to entering targets and observations into the 
APT<sup id="a3">[[3]](#f1)</sup>. 

We stress that the CVT is designed to provide quick illustrations of the possible observable orientations for a given 
target. As such, the CVT rapidly approximates JWST’s pointing restrictions and does not query the official JWST Proposal 
Constraint Generator (PCG), nor include detailed pointing restrictions like Earth and Moon avoidance, etc. Therefore, 
CVT results should be treated as useful approximations that may differ from official APT constraints by a degree or so.

Documentation can be found online at [JWST Coronagraph Visibility Tool Help](https://jwst-docs.stsci.edu/jwst-other-tools/jwst-target-visibility-tools/jwst-coronagraphic-visibility-tool-help).

**Authors:** Christopher Stark, Joseph Long, J. Brendan Hagan, Mees Fix and Bryony Nickson

<p align="center">
  <img src="screenshot.png" alt="Screenshot of the JWST Coronagraph Visibility Tool showing target HR 8799 with three companions plotted."/>
</p>

<a name="user-install"></a>
## Installation for Users 

### Installing the Python Package

#### Installing with pip

CVT may be installed from the [Python Package Index](https://pypi.org/) in the usual manner for Python packages.
 
    $ pip install jwst_coronagraph_visibility 

### Installing the macOS application

If you are running macOS and would like a double-clickable application, click on the following link:<br> **[Download for macOS (86.4 MB)](https://github.com/spacetelescope/jwst_coronagraph_visibility/releases/download/0.4.4/jwst_coronagraph_visibility_tool_macos.zip)**. 

Simply extract the downloaded zip file to obtain the .app bundle, then double-click to run the JWST Coronagraph Visibility Tool.

## Installation for Contributors

For those wishing to contribute to the code base, you can install `jwst_coronagraph_visibility` by cloning and 
installing the repository. This is only recommended for those looking to help with development. In general, those 
wishing only to use the jwst_coronagraph_visibility tool should install the latest stable version from using Astroconda, 
as described in the [instructions above](#user-install).

### Prerequisites 

It is highly recommended that contributors have a working installation of [Miniconda](https://conda.io/miniconda.html) 
or [Anaconda](Anaconda) for Python 3.6. Package requirements for contributing to `jwst_coronagraph_visibility` will be 
provided by a `setup.py` script included in the repository. 

### Clone the repository:

Clone the `jwst_coronagraph_visibility` GitHub repository as follows:

    $ git clone https://github.com/brynickson/jwst_coronagraph_visibility.git
    $ cd jwst_coronagraph_visibility

### Environment Installation 

Following the download of the `jwst_coronagraph_visibility` repository, create and activate a new 
`jwst_coronagraph_visibility` environment:

    $ conda create -n jwst_coronagraph_visibility-3.7 python=3.7
    $ conda activate jwst_coronagraph_visibility-3.7
    
### Package installation

Next, you need to install the `jwst_coronagraph_visibility` package. This can be accomplished by running the `setup.py` 
script:

    (jwst_coronagraph_visibility-3.7)$ python setup.py develop
    
The package should now appear if you run `conda list jwst_coronagraph_visibility`. 


## Citation 

If you use the CVT for work/research presented in a publication (whether directly, or as a dependency to another 
package), please consider citing the Zenodo record using the DOI page above. Please find additional instructions in 
[CITATION](CITATION).


## Software Contributions

Contributors should use a ["forking workflow"](https://github.com/spacetelescope/style-guides/blob/master/guides/git-workflow.md#the-forking-workflow-) 
when making contributions to the project. 

## Code of Conduct 

Users and contributors to the `jwst_coronagraph_visibility` repository should adhere to the 
[Code of Conduct](CODE_OF_CONDUCT). Any issues or violations pertaining to the Code of Conduct should be brought to 
the attention of a `jwst_coronagraph_visibility` team member or to `conduct@stsci.edu`.

## Questions

For any questions about the `jwst_coronagraph_visibility` project or its software or documentation, please 
[open an Issue](https://github.com/spacetelescope/jwst_coronagraph_visibility/issues).

## Known Issues

  * The CVT does not (and will not) query the JWST Proposal Constraint Generator. The only constraint on the field of regard is the Sun and anti-Sun avoidance angle.
  * Target name resolution depends on the availability of the SIMBAD service. If the service cannot be reached, you will have to enter coordinates yourself.
  * The CVT does not currently provide a way to export the plotted points as text. Plots can be saved from the GUI using the save icon below the plot panel.
  * The CVT has only been tested on Mac and Linux. Issue reports from Windows users are welcome, and we will do our best to address them, but we are not testing the tool on Windows.

*See issue tracker at* https://github.com/spacetelescope/jwst_coronagraph_visibility/issues.

## Current Development Team
- Mees Fix [@mfixstsci](https://github.com/mfixstsci)
- Bryony Nickson [@brynickson](https://github.com/brynickson)
<br>

##  Acronyms 
<b id="f1">[1]</b> - Mid-Infrared Instrument (see [documentation](https://jwst-docs.stsci.edu/mid-infrared-instrument/miri-observing-modes/miri-coronagraphic-imaging)) [ ↩](#a1) <br>
<b id="f1">[2]</b> - Near-Infrared Instrument (see [documentation](https://jwst-docs.stsci.edu/near-infrared-camera/nircam-observing-modes/nircam-coronagraphic-imaging)) [ ↩](#a2) <br>
<b id="f1">[3]</b> - Astronomer's Proposal Tool (see [documentation](https://jwst-docs.stsci.edu/jwst-astronomers-proposal-tool-overview)) [ ↩](#a2)

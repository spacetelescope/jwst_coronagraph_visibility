# Releases

## 0.4.5 

### New Features 

#### Project & API Documentation

 - Added project citation information to `README` and `CITATION`, along with a Zenodo badge.
 - Updated `README` to include: links to JDOX documentation; installation instructions for both general users and 
 contributors; link to `CODE_OF_CONDUCT`; contact information for current developers and provided link to the appropriate
 git workflow.
 - Added release notes (`CHANGES`)
 - Added a documented release procedure (`RELEASES`) 
 
## 0.4.4 (2020-02-13)

### New Features

#### Repository

- Updated to use `pysiaf` version `0.7.1` which uses **PRD PRDOPSSOC-027**
- Updated unit tests to reflect changes from latest SIAF in **PRDOPSSOC-027**

### Bug Fixes

#### Mac OS app bundle

 - Fixed Mac bundle issue that was causing application not to open


## 0.4.2 (2020-01-16)

### Bug Fixes

#### Mac OS app bundle 

 - Fixed bug in Mac OS app bundle
 
 
## 0.4.0 (2020-01-16)
 
### New Features
 
#### Repository
 
 - Added use of version `0.6.3` of `pysiaf` which uses **PRDOPSSOC-M-026** as default for the JWST SIAF.
 
### Bug Fixes

#### Repository

 - Discontinued use of the `jwxml` package. Tool now uses the STScI supported `pysiaf` package for all information using the SIAF. 
 - Discontinued support for python 2.7*
 
 
## 0.3.0 (2017-07-24)
 
### Bug Fixes
 
#### Repository

 - Fixed bug in GUI so that NIRCam A long-wavelength bar mask is now oriented correctly in GUI (previously flipped left-to-right)
 - Fixed bug in GUI such that The SIMBAD search field is now cleared when user enters RA/Dec or chooses an example
 

## 0.2.0 (2017-03-12)

### New Features

#### Project & API Documentation

 - Updated development instructions in `README`.
 
#### Repository 

 - Incorporated revised definitions of the science instrument apertures in the SIAF (version **PRDOPSSOC-F-008**)
 - Added MIRI target acquisition locations on the detector plot 
 - Updated the MIRI TA spots to be translucent and renamed by APT numbers
 - Added a 'zoom to fit' button to the GUI
 - Added day of year to plot overlay in GUI
 
 ### Bug Fixes
 
#### Repository 

 - Changed visibility calculation plot such that it now starts on Jan 1 instead of Oct 1 
 
## 0.1.0 (2017-01-17)
 
### New Features 

#### Project and API Documentation

 - Added a `LICENSE` file has been to reflect availability under the 3-Clause BSD license.

#### Repository

 - Renamed tool to `jwst_coronagraph_visibility` with the command to launch the GUI taking on the name `jwst-coronagraph-visibility-gui`
 - Updated `jwxml` to **PRDOPSSOC-E-002** (version `0.2.0`), which is required for new aperture names used by this tool
 - Added NIRCam coronagraph ND squares and clips in the right hand (detector) plot in GUI


## 0.0.2 (2016-11-03)

### New Features

#### Repository 

- SIAF is now bundled with the `jwxml` package, and targets **PRDDEVSOC-D-012**. This includes minor refinements to the transformations from sky coordinates to the coronagraph aperture coordinates.
 - Added controls for the sampling of the roll angle and time of year to the GUI

### Bug Fixes

#### Repository

- Removed NIRCam Module B is as a selectable instrument

## 0.0.1 (2016-08-01)

### New Features

#### Repository 

- First tagged release of the tool for internal and external testing. Target visibility plots from this tool should always be checked against APT for consistency, as this tool does not account for all of the same constraints (and is not intended to).


 
 
 
 
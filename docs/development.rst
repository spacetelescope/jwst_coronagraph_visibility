*********************************
Developing and releasing the tool
*********************************

Ensure that pytest and Sphinx are installed.

Building the docs
=================

Documentation is reST processed by Sphinx. Generate the HTML copy of the documentation with::

   cd docs/
   make html

and open ``_build/html/index.html`` to see the result.

AstroConda
==========

*Instructions for updating the tool in AstroConda will be added after the tool is packaged for AstroConda.*

App bundles for macOS
=====================

To build the Mac ``.app`` bundle, run ``make app`` from the repository root with PyInstaller and pytest installed.

To build a maximally backwards compatible bundle using the ``banana`` build machine at STScI::

   ssh banana
   mkdir -p /tmp/$USER
   export HOME=/tmp/$USER
   export PATH="/tmp/$USER/Library/Python/2.7/bin:$PATH"
   export PYTHONPATH="/tmp/$USER/Library/Python/2.7/lib/python/site-packages:$PYTHONPATH"
   cd
   pip install --user -U numpy
   pip install --user -U matplotlib
   git clone https://github.com/spacetelescope/jwst_visibility.git
   cd jwst_visibility
   git pull
   pip install --user .
   make app
   cd dist
   zip -r jwst_visibility_calculator_osx.zip ./JWST\ Visibility\ Calculator.app/
   mv jwst_visibility_calculator_osx.zip /user/$USER/

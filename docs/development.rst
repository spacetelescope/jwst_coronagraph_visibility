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

Preparing to release
====================

1. Ensure that the test suite passes: ``pytest``.
2. Remove ``dev`` from the version number in ``setup.py``.
3. Commit and tag: ``git commit -am "Releasing 0.x.y"``, ``git tag 0.x.y``.
4. Push the tag: ``git push --tags origin`` or ``git push --tags upstream``.
5. Update the version number in ``setup.py`` and add ``dev`` to the end.
6. Commit: ``git commit -am "Back to development: 0.x.y"``.
7. Push changes to GitHub.

Releasing
=========

To release on PyPI
------------------

1. Follow the `instructions from PyPA <https://packaging.python.org/distributing/#uploading-your-project-to-pypi>`_ to set up ``twine`` and login credentials.
2. Check out the just-tagged commit: ``git checkout 0.x.y``
3. Create the sdist: ``python setup.py sdist``
4. Create a universal wheel: ``python setup.py bdist_wheel --universal``
5. Check distribution is valid: ``twine check dist/*``
6. Test upload with ``twine upload --repository-url https://test.pypi.org/legacy/ dist/*``
7. Upload to PyPI: ``twine upload dist/*``

To release through AstroConda
-----------------------------

1. Fork `astroconda/astroconda-contrib <https://github.com/astroconda/astroconda-contrib>`_ to your GitHub, if you haven't already.
2. Clone your fork locally (or pull upstream changes, if you already have a local clone).
3. Check out a feature branch like this: ``git checkout -b jwst_coronagraph_visibility-0.x.y``.
4. Update the version number in ``astroconda-contrib/jwst_coronagraph_visibility/meta.yaml`` (along with any version requirements on upstream dependencies).
5. Commit your changes.
6. Push to GitHub.
7. Open a pull request against `astroconda/astroconda-contrib <https://github.com/astroconda/astroconda-contrib>`_.

App bundles for macOS
=====================

To build the macOS ``.app`` bundle, run ``make app`` from the repository root with PyInstaller and pytest installed.

To build a maximally backwards compatible bundle using the ``banana`` build machine at STScI::

   ssh banana
   mkdir -p /tmp/$USER
   export HOME=/tmp/$USER
   export PATH="/tmp/$USER/Library/Python/2.7/bin:$PATH"
   export PYTHONPATH="/tmp/$USER/Library/Python/2.7/lib/python/site-packages:$PYTHONPATH"
   cd
   curl -O https://bootstrap.pypa.io/get-pip.py
   python get-pip.py --user
   pip install --user -U pyinstaller
   pip install --user -U numpy
   pip install --user -U matplotlib
   git clone https://github.com/spacetelescope/jwst_coronagraph_visibility.git
   cd jwst_coronagraph_visibility
   git pull
   pip install --user .
   make app
   cd dist
   zip -r jwst_coronagraph_visibility_tool_macos.zip ./JWST\ Coronagraph\ Visibility\ Tool.app/
   mv jwst_coronagraph_visibility_tool_macos.zip /user/$USER/

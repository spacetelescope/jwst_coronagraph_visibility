********************************************************
Using the Python API for the JWST Target Visibility Tool
********************************************************

The current Python API is a direct translation of the original IDL Coronagraph Visibility Tool for JWST by Christopher Stark. The function in which the bulk of the calculation happens is ``jwst_visibility.skyvec2ins.skyvec2ins`` (in ``jwst_visibility/skyvec2ins.py``).

.. autofunction:: jwst_visibility.skyvec2ins.skyvec2ins

Should there be user demand, a more conventional Python API may be developed, but the current plan is to support the GUI as the primary interface.

Automated software testing
==========================

To ensure the correctness of the Python to IDL translation, the ``jwst_visibility/tests/targets/`` folder contains CSV files output by the IDL tool for the ``skyvec2ins`` procedure. The Python code contains automated tests that run the Python ``skyvec2ins`` function and compare the output.

To run the test suite, install ``pytest`` and run the command::

    $ py.test

from the root folder of this repository.

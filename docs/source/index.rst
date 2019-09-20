.. mkstuff documentation master file, created by
   sphinx-quickstart on Mon Oct 24 16:46:22 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation of the `mkstuff` package
======================================

.. include:: install.rst

Content
-------

Overview

* `mkstuff`
   * `mkstuff.py`
      General library and helper functions
   * `mkplot.py`
      For plotting
   * `athena.py`
      To deal with `athena` (http://www.cosmostat.org/software/athena) input and output
   * `populationWL.py`
      Population weak lensing (lensing by foreground objects)

.. toctree::
   :numbered:
   :maxdepth: 3
   :caption: General

   examples
   mkstuff

.. toctree::
   :maxdepth: 3
   :caption: Executable scripts

   ascii2fits
   fits2ascii

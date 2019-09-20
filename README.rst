|Travis|_ |Coveralls|_

.. |Travis| image:: https://api.travis-ci.org/martinkilbinger/mkstuff.svg?branch=master
.. _Travis: https://travis-ci.org/martinkilbinger/mkstuff


.. |Coveralls| image:: https://coveralls.io/repos/github/martinkilbinger/mkstuff/badge.svg?branch=master
.. _Coveralls: https://coveralls.io/github/martinkilbinger/mkstuff?branch=master



The `mkstuff` package
=====================

Library for random scientific python modules and scripts.

:Author: Martin Kilbinger `(martin.kilbinger@cea.fr) <martin.kilbinger@cea.fr>`_

:Version: 1.2.0

:Date: April 2019

:Documentation: |link-to-docs|

.. |link-to-docs| raw:: html

  <a href="https://martinkilbinger.github.io/mkstuff"
  target="_blank">https://martinkilbinger.github.io/mkstuff</a>


.. include:: docs/source/install.rst

Content
-------

Modules, in `mkstuff/mkstuff`
        * `mkstuff.py`
          General library and helper functions
        * `mkplot.py`
          For plotting
        * `athena.py`
          To deal with `athena` (http://www.cosmostat.org/software/athena) input and output
        * `populationWL.py`
          Population weak lensing (lensing by foreground objects)

Executable scripts, in `mkstuff/bin`
        * `fits2ascii.py`
        * `ascii2fits.py`

|Travis|_ |Coveralls|_
#|Python35|_ |Python36|_ |Python37|_ |PyPi|_

.. |Travis| image:: https://api.travis-ci.org/martinkilbinger/mkstuff.svg?branch=master
.. _Travis: https://travis-ci.org/martinkilbinger/mkstuff

.. |Coveralls| image:: https://coveralls.io/github/martinkilbinger/mkstuff/badge.svg?branch=master&service=github
.. _Coveralls: https://coveralls.io/github/martinkilbinger/mkstuff

#.. |Python35| image:: https://img.shields.io/badge/python-3.5-blue.svg
#.. _Python35: https://badge.fury.io/py/python-pySAP

#.. |Python36| image:: https://img.shields.io/badge/python-3.6-blue.svg
#.. _Python36: https://badge.fury.io/py/python-pySAP

#.. |Python37| image:: https://img.shields.io/badge/python-3.7-blue.svg
#.. _Python37: https://badge.fury.io/py/python-pySAP

#.. |PyPi| image:: https://badge.fury.io/py/python-pySAP.svg
#.. _PyPi: https://badge.fury.io/py/python-pySAP



mkstuff
=======

Library for random scientific python modules and scripts.

:Author: Martin Kilbinger `(martin.kilbinger@cea.fr) <martin.kilbinger@cea.fr>`_

:Version: 1.2.0

:Date: April 2019

:Documentation: |link-to-docs|

.. |link-to-docs| raw:: html

  <a href="https://martinkilbinger.github.io/mkstuff"
  target="_blank">https://martinkilbinger.github.io/mkstuff</a>


Installation
------------

First, clone the ``mkstuff`` repository from GitHUB using the command

.. code-block:: sh

        git clone https://github.com/martinkilbinger/mkstuff.git

Next, run the ``setup.py`` script in the ``mkstuff`` directory

.. code-block:: sh

        cd mkstuff
        python setup.py install [OPTIONS]

OPTIONS can for example specify installation to a specific directory with ``--prefix DIR``. Type

.. code-block:: sh

        python setup.py install --help

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

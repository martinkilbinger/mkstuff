|Travis|_

.. |Travis| image:: https://api.travis-ci.org/martinkilbinger/mkstuff.svg?branch=master
.. _Travis: https://travis-ci.org/martinkilbinger/mkstuff
.. image:: https://coveralls.io/repos/github/martinkilbinger/mkstuff/badge.svg?branch=master
:target: https://coveralls.io/github/martinkilbinger/mkstuff?branch=master



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

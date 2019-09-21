.. mkstuff documentation master file, created by
   sphinx-quickstart on Mon Oct 24 16:46:22 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation of the `mkstuff` package
======================================

.. about

About
^^^^^

Library for random scientific python modules and scripts.

:Author: Martin Kilbinger `(martin.kilbinger@cea.fr) <martin.kilbinger@cea.fr>`_

:Version: 1.2.0

:Date: April 2019

:Documentation: |link-to-docs|

.. |link-to-docs| raw:: html

  <a href="https://martinkilbinger.github.io/mkstuff"
  target="_blank">https://martinkilbinger.github.io/mkstuff</a>

Installation
^^^^^^^^^^^^

First, clone the ``mkstuff`` repository from `github <https://github.com/>`_ using the command

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

.. overview.rst

Overview
^^^^^^^^

* `mkstuff`
   * `mkstuff.py`
      General library and helper functions
   * `mkplot.py`
      For plotting
   * `athena.py`
      To deal with `athena` (http://www.cosmostat.org/software/athena) input and output
   * `populationWL.py`
      Population weak lensing (lensing by foreground objects)
* `bin`
   Executable scripts

Detailed documentation
----------------------

mkstuff package
===============

Subpackages
-----------

.. toctree::

    mkstuff.example
    mkstuff.tests

Submodules
----------

mkstuff.athena module
---------------------

.. automodule:: mkstuff.athena
    :members:
    :undoc-members:
    :show-inheritance:

mkstuff.info module
-------------------

.. automodule:: mkstuff.info
    :members:
    :undoc-members:
    :show-inheritance:

mkstuff.mkplot module
---------------------

.. automodule:: mkstuff.mkplot
    :members:
    :undoc-members:
    :show-inheritance:

mkstuff.mkstuff module
----------------------

.. automodule:: mkstuff.mkstuff
    :members:
    :undoc-members:
    :show-inheritance:

mkstuff.populationWL module
---------------------------

.. automodule:: mkstuff.populationWL
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: mkstuff
    :members:
    :undoc-members:
    :show-inheritance:

.. executable scripts

.. toctree::
   :maxdepth: 4
   :caption: Executable scripts

   ascii2fits
   fits2ascii


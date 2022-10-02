Quickstart
==========

Installation
------------

First make sure that you have python3 and pip installed.
Then execute the following in your terminal:

.. code-block:: console

   $ pip install rad-tools

.. _update:

Update
------

.. code-block:: console

   $ pip install rad-tools --upgrade

TB2J plotter
-------------

Display help message:

.. code-block:: python

    tb2j_plotter.py -h

Plot all distances and all isotropic exchange parameters 
from *exchange.out* file:

.. code-block:: python

    tb2j_plotter.py -f exchange.out

Plot all distances and all isotropic exchange parameters 
from *exchange.out* file for the bonds with distance <= 5 Angstrom:

.. code-block:: python

    tb2j_plotter.py -f exchange.out -maxd 5

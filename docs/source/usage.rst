Quickstart
==========

Installation
------------

First make sure that you have python3 and pip installed.
Then execute the following in your terminal:

.. code-block:: console

   pip install rad-tools

.. _update:

Update
------

If you want to update the package to the last available version (|release|)
run the following in your terminal:

.. code-block:: console

   pip install rad-tools --upgrade

TB2J plotter
-------------

Whenever you forget how to use the script do not hesitate and display the help message!

.. code-block:: python

    tb2j_plotter.py -h

Simpliest usage example: plot all distances and all isotropic exchange parameters 
from *exchange.out* file:

.. code-block:: python

    tb2j_plotter.py -f exchange.out

Plot all distances and all isotropic exchange parameters 
from *exchange.out* file for the bonds with distance <= 5 Angstrom:

.. code-block:: python

    tb2j_plotter.py -f exchange.out -maxd 5

Phonopy plotter
---------------

Whenever you forget how to use the script do not hesitate and display the help message!

.. code-block:: python

    phonopy_plotter.py -h

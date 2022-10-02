Script Collection
=================
Various scripts which were useful during my PhD.

Installation
------------

First make sure that you have python3 and pip installed.
Then execute the following in your terminal:

.. code-block:: console

   pip install rad-tools

tb2j_plotter.py `docs <https://github.com/adrybakov/rad-tools/blob/master/doc/tb2j_plotter.rst>`_
============================================================================================
Script for visualisation of TB2J results.

Display Isotropic exchange and distances (one output file for each). 
Output files will have the following name structure: 
*output-name.display_data_type.png*

Currently filtering by R vectors, distances and template file 
is supported.

``--filename`` (or ``-f``) argument is required, the rest of them are optional.


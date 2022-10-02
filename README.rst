Script Collection
=================
Various scripts which were useful during my PhD.

.. image:: https://readthedocs.org/projects/rad-tools/badge/?version=latest
    :target: https://rad-tools.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Installation
------------

First make sure that you have python3 and pip installed.
Then execute the following in your terminal:

.. code-block:: console

   pip install rad-tools

tb2j_plotter.py `docs <https://rad-tools.adrybakov.com/en/latest/tb2j_plotter.html>`_
============================================================================================
Script for visualisation of TB2J results.

Display Isotropic exchange and distances (one output file for each). 


Supports filtering by R vectors, distances and template file. 
The result will be defined by logical conjugate of the specified conditions.

``--filename`` (or ``-f``) argument is required, the rest of them are optional.

Output files will have the following name structure: 
*output-name.display_data_type.png*

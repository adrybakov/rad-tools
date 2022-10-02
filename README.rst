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

Currently filtering by 
R vectors (see :ref:`R-vector <tb2j_plotter_R-vector>`), 
distances (see :ref:`max-distance <tb2j_plotter_max-distance>`,
:ref:`min-distance <tb2j_plotter_min-distance>` and
:ref:`distance <tb2j_plotter_distance>`), 
and template file (see :ref:`template <tb2j_plotter_template>`), is supported. 
The result will be defined by logical conjugate of the specified conditions.

``--filename`` (or ``-f``) argument is required, the rest of them are optional.

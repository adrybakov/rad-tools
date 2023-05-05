.. _user-guide:

********************
RAD-tools user guide
********************

The package covers post-processing scenarios for the 
results of |QE|_, |TB2J|_ and |Wannier90|_, as well as provides some custom scripts.

It is expected to be used in two ways:

* :doc:`As Python library <../api/index>`
* Via :doc:`scripts <scripts>` (i.e. usage from command line)

Scripts
=======

For any script use ``--help`` or ``-h`` option in console in order to display 
help message with the short summary of the arguments.

.. code-block:: console
   
   script-name.py -h

All the files used in the usage examples are available :examples:`here <>`.

.. toctree::
   :caption: for Quantum Espresso:
   :maxdepth: 1
   :hidden:

   rad-plot-dos

.. toctree::
   :caption: for TB2J:
   :maxdepth: 1


   rad-make-template
   tb2j-extractor
   tb2j-plotter

.. toctree::
   :caption: for Wannier90:
   :maxdepth: 1

   identify-wannier-centres

.. toctree::

   output-notes



   
There are no more guides. You are now guideless. Good luck.

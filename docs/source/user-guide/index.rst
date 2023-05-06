.. _user-guide:

********************
RAD-tools user guide
********************

The package covers post-processing scenarios for the 
results of |QE|_, |TB2J|_ and |Wannier90|_, as well as provides some custom scripts.

It is expected to be used in two ways:

* :doc:`As Python library <../api/index>`
* :ref:`Via scripts <scripts>` (i.e. usage from command line)

.. _scripts:

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
   rad-extract-tb2j
   rad-plot-tb2j

.. toctree::
   :caption: for Wannier90:
   :maxdepth: 1

   rad-identify-wannier-centres

.. toctree::

   output-notes



   
There are no more guides. You are now guideless. Good luck.

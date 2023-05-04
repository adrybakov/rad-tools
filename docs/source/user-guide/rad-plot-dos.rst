.. _rad-plot-dos:

**********************
``rad-plot-dos.py``
**********************

.. note::
    Name changed in 0.5.21 from ``rad-dos-plotter.py`` to ``rad-plot-dos.py``

Script for visualisation of projected density of states from 
`Quantum-ESPRESSO <https://www.quantum-espresso.org/>`_ 
(`QE input description <https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html>`_).

The script relies on the QE output format and may not work if the QE  output files 
are modified or renamed. Minimum input for the script is the folder 
with QE output files inside:

.. note::
    Currently script only supports collinear cases and non-collinear non-spinorbit case.

Usage example
=============
Minimal input looks like the following:

.. code-block:: bash

    rad-plot-dos.py -if collinear

where "collinear" is a path to the folder with output files from QE PDOS calculations.

If you want to choose particular energy window use an 
option :ref:`--window <rad-plot-dos_energy-window>`:

.. code-block:: bash

    rad-plot-dos.py -if collinear -w -10 5


Arguments
=========

.. _rad-plot-dos_input-path:

-ip, --input-path
-----------------
Relative or absolute path to the folder with dos files.

.. code-block:: text

    required


.. _rad-plot-dos_seedname:

-s, --seedname
--------------
Prefix for output files containing PDOS(E). 
As specified in the QE projwfc.x input file.

If it is not provided the script will try to 
detect it automatically in the "input-path" folder.

.. code-block:: text

    default : None

Renamed in 0.5.21: from "filpdos" to "seedname".


.. _rad-plot-dos_output-path:

-op, --output-path
------------------
Relative or absolute path to the folder for saving outputs.

.. code-block:: text

    default : current directory (".")


.. _rad-plot-dos_energy-window:

-ew, --energy-window
--------------------
Energy window for the plots.  
By default whole range present in the files is plotted.

.. code-block:: text

    default : None

Renamed in 0.5.21: from "window" to "energy-window".


.. _rad-plot-dos_dos-window:

-dw, --dos-window
-----------------
DOS window for the plots.  
By default whole range present in the files is plotted.

.. code-block:: text

    default : None

Added in 0.5.21.


.. _rad-plot-dos_efermi:

-ef, --efermi
-------------
Fermi energy. If specified zero will be shift to Fermi energy.

.. code-block:: text

    default : 0


.. _rad-plot-dos_separate:

-sep, --separate
----------------
Whenever to plot projected DOS for each atom  of the same type separately.

.. code-block:: text

    default : False


.. _rad-plot-dos_relative:

-r, --relative
--------------
Whenever to use relative style.

.. code-block:: text

    default : False


.. _rad-plot-dos_normalize:

-n, --normalize
---------------
Whenever to use normalize relative style.

.. code-block:: text

    default : False


.. _rad-plot-dos_verbose:

-v, --verbose
-------------
Verbose output, propagates to the called methods.

.. code-block:: text

    default : False


.. _rad-plot-dos_interactive:

-i, --interactive
-----------------
Interactive plotting.

.. code-block:: text

    default : False


.. _rad-plot-dos_save-pickle:

-sp, --save-pickle
------------------
Whenever to save figures as .pickle files.

.. code-block:: text

    default : False

Added in 0.5.21.


.. _rad-plot-dos_save-txt:

-st, --save-txt
---------------
Whenever to save some data as txt files.

.. code-block:: text

    default : False

Added in 0.5.21.

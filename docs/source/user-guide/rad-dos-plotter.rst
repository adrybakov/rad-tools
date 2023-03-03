.. _rad-dos-plotter:

**********************
``rad-dos-plotter.py``
**********************

Script for visualisation of projected density of states from 
`Quantum-ESPRESSO <https://www.quantum-espresso.org/>`_ 
(`QE input description <https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html>`_).

The script rely on the QE output format and may not work if the QE  output files 
are modified or renamed. Minimum input for the script is the folder 
with QE output files inside:

Arguments
=========

.. _rad-dos-plotter_input-path:

-ip, --input-path
-----------------
Relative or absulute path to the folder with dos files.

.. code-block:: text

    required


.. _rad-dos-plotter_filpdos:

-f, --filpdos
-------------
Prefix for output files containing PDOS(E). 
As specified in the QE projwfc.x input file.

If it is not provided the script will try to 
detect it automatically in the "input-path" folder.

.. code-block:: text

    default : None


.. _rad-dos-plotter_output-path:

-op, --output-path
------------------
Relative or absolute path to the folder for saving outputs.

.. code-block:: text

    default : current directory (".")


.. _rad-dos-plotter_window:

-w, --window
------------
Energy window for the plots.  
By default whole range present in the files is plotted.

.. code-block:: text

    default : None


.. _rad-dos-plotter_efermi:

-ef, --efermi
-------------
Fermi energy. If specified zero will be shift to Fermi energy.

.. code-block:: text

    default : 0


.. _rad-dos-plotter_separate:

-s, --separate
--------------
Whenever to plot projected DOS for each atom  of the same type separately.

.. code-block:: text

    default : False


.. _rad-dos-plotter_relative:

-r, --relative
--------------
Whenever to use relative style.

.. code-block:: text

    default : False


.. _rad-dos-plotter_normalize:

-n, --normalize
---------------
Whenever to use normalize relative style.

.. code-block:: text

    default : False

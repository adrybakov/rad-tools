.. _rad-dos-plotter:

**********************
``rad-dos-plotter.py``
**********************

Script for visualisation of projected density of states from 
`Quantum-ESPRESSO <https://www.quantum-espresso.org/>`_.

Detect the labels based on the name of the file and print its content.
Pdos for *up* states is plotted with positive sign, 
Pdos for *down* states is plotted with negative sign. 
Summed Pdos from the file (2 and 3 column) is plotted in grey background.

Arguments
=========

.. _rad-dos-plotter_filename:

-f, --filename
--------------
Relative or absolute path to the dos file,
including the name and extention of the file itself.

.. code-block:: text

    required

.. _rad-dos-plotter_output-dir:

-op, --output-dir
-----------------
Relative or absolute path to the folder for saving outputs.

.. code-block:: text

    default : current directory

See also: :ref:`example <scripts_output-notes>`.


.. _rad-dos-plotter_output-name:

-on, --output-name
------------------
Seedname for the output files.

.. code-block:: text

    default : pdos

See also: :ref:`example <scripts_output-notes>`.


.. _rad-dos-plotter_window:

-w, --window
------------
Energy window.

DOS will be plotted in the energy limits ``window[0]< E < window[1]``

.. code-block:: text

    default : None

.. _rad-dos-plotter_interactive:

-i, --interactive
------------------
Interactive mode flag.

If specified then :py:func:`plt.show()` function will be called 
instead of saving graph to the file, which will open standart 
pyplot interactive window.

.. code-block:: text

    default : False

.. _rad-dos-plotter_legend-location:

-ll, --legend-location
----------------------
Legend location, will be passed to the plt.legend().

.. code-block:: text

    default "best"

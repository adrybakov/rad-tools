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

    *required* : True

    *type* : str


.. _rad-dos-plotter_output-dir:

-op, --output-dir
-----------------
Relative or absolute path to the folder for saving outputs.

    *default* : current directory

    *type* : str

See also: :ref:`example <scripts_output-notes>`.


.. _rad-dos-plotter_output-name:

-on, --output-name
------------------
Seedname for the output files.

    *default* : pdos

    *type* : str

See also: :ref:`example <scripts_output-notes>`.


.. _rad-dos-plotter_window:

-w, --window
------------
Energy window.

DOS will be plotted in the energy limits ``window[0]< E < window[1]``

    *default* : None

    *type* : float

    *nargs* : 2


.. _rad-dos-plotter_interactive:

-i, --interactive
------------------
Interactive mode flag.

If specified then :py:func:`plt.show()` function will be called 
instead of saving graph to the file, which will open standart 
pyplot interactive window.

    *default* : False

    *action* : store_true

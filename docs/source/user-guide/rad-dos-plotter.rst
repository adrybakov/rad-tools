.. _rad-dos-plotter:

**********************
``rad-dos-plotter.py``
**********************

Script for visualisation of projected density of states.

Currently plot the content of separate file. Pdos for *up* states is
be plotted with positive sign, Pdos for *down* states is plotted with 
negative sign. Summed Pdos from the file (2 and 3 column) is
plotted in grey background.

Arguments
=========

.. _rad-dos-plotter_filename:

``--filename``, ``-f``
----------------------
Relative or absulute path to the dos file,
including the name and extention of the file itself.

    *required* : True

    *type* : str


.. _rad-dos-plotter_output-dir:

``--output-dir``, ``-op``
-------------------------
Relative or absolute path to the folder for saving outputs.

    *default* : current directory

    *tyoe* : str


.. _rad-dos-plotter_output-name:

``--output-name``, ``-on``
--------------------------
Seedname for the output files.

    *default* : pdos

    *type* : str


.. _rad-dos-plotter_window:

``--window``, ``-w``
--------------------
Energy window.

DOS will be plotted in the energy limits ``window[0]< E < window[1]``

    *default* : None

    *type* : float

    *nargs* : 2


.. _rad-dos-plotter_interactive:

``--interactive``, ``-i``
-------------------------
Interactive mode flag.

If specified then :py:func:`plt.show()` function will be called 
instead of saving graph to the file, which will open standart 
pyplot interactive window with the plot.

    *default* : False

    *action* : store_true
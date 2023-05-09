.. _rad-plot-dos:

**********************
``rad-plot-dos.py``
**********************

.. versionchanged:: 0.5.21 Renamed from ``rad-dos-plotter.py``

Script for visualization of projected density of states.
It supports the result of the calculation from:

* |QE|_ (|projwfc|_).
    All type of calculations are supported 
    (Collinear, spin-unpolarized; 
    collinear, spin-polarized; 
    non-collinear, non-spin-orbit; 
    non-collinear, spin-orbit)

By default the script sums PDOS over the atoms of the same type 
(and divide by the amount of atoms of the same type). 
Option :ref:`--separate <rad-plot-dos_separate>` plots PDOS individually for each atom.

K-resolved PDOS are handled by summing over kpoints (and dividing by the number of kpoints).

Output data
===========
Script tries to detect all seednames present 
in :ref:`--input-path <rad-plot-dos_input-path>` by default. 

:ref:`--seedname <rad-plot-dos_seedname>` option specifies 
particular seedname to work with.

For each seedname a separate folder "seedname-suffix" is created.

Suffix is a combination of any number of the following words:
    * "separate" - appears if :ref:`--separate <rad-plot-dos_separate>` option is used.
    * "relative" - appears if :ref:`--relative <rad-plot-dos_relative>` option is used.
    * "normalized" - appears if :ref:`--normalize <rad-plot-dos_normalize>` option is used.
    * "vstotal" - appears if :ref:`--background-total <rad-plot-dos_background-total>` option is used.

The structure of the :ref:`output folder <rad-plot-dos_output-path>` is the following:

.. code-block:: text

    output_path/
    ├── ....
    ├── seedname_1/
    ├── ...
    └── seedname_n/

Each seedname folder has the structure:

.. code-block:: text

    seedname/
    ├── pdos-vs-dos.png
    ├── atomic-contributions.png
    └── atom-resolved/
        ├── output_name_1
        └── output_name_2
    └── orbital-resolved/
        ├── output_name_1
        └── output_name_2

By default all possible graphs are plotted:

* pdos-vs-dos.png
    Total partial density of states (sum over all projectors) vs 
    total density of states (directly from the plane-wave basis).
    Affected by :ref:`--save-pickle <rad-plot-dos_save-pickle>`.
* atomic-contributions.png
    Contribution of each atom (summed over all projectors) 
    into the total partial density of states.
    Affected by :ref:`--save-pickle <rad-plot-dos_save-pickle>`.
    Affected by :ref:`--save-txt <rad-plot-dos_save-txt>`.
* atom-resolved/
    Contribution of each projectors group (i.e. :math:`s`, :math:`p`, :math:`d`, :math:`f`) 
    into the partial density of state of each atom.
    Affected by :ref:`--save-pickle <rad-plot-dos_save-pickle>`.
    Affected by :ref:`--save-txt <rad-plot-dos_save-txt>`.
* orbital-resolved/
    Contribution of each projector (i. e. :math:`p_z`, :math:`p_x`, :math:`p_y`) into the total 
    partial density of states of each group (i.e. :math:`p`).
    Affected by :ref:`--save-pickle <rad-plot-dos_save-pickle>`.
    Affected by :ref:`--save-txt <rad-plot-dos_save-txt>`.

By default only the pictures (.png) are created. Two additional formats of the output are:

* txt
    Content of the plots in txt format. First line is the header with projectors.
* pickle
    Python-specific format, which allowed to pick up the ``figure`` 
    from the python code and modify it:

    .. code-block:: python

        import pickle
        import matplotlib.pyplot as plt

        fig = pickle.load(open('filename.pickle', 'rb'))
        axes = fig.get_axes()

        for ax in axes:
            ax.set_xlabel("Custom x label")
            ax.set_ylabel("Custom y label")
            ax.set_title("Custom title")

        fig.savefig("filename.png", dpi=400, bbox_inches="tight")

    If ``fig.show()`` or ``plt.show()`` does not work the following fix may help
    (`credit <https://stackoverflow.com/a/54579616>`_):

    .. code-block:: python

        def show_figure(fig):
            dummy = plt.figure()
            new_manager = dummy.canvas.manager
            new_manager.canvas.figure = fig
            fig.set_canvas(new_manager.canvas)

        show_figure(fig)
        plt.show()

Usage example
=============
Minimal possible input is:

.. code-block:: bash

    rad-plot-dos.py 

It will try to detect PDOS output files in the current directory and plot them.

To choose energy window use an 
option :ref:`--energy-window <rad-plot-dos_energy-window>`:

.. code-block:: bash

    rad-plot-dos.py -ew -10 5

To choose :ref:`input <rad-plot-dos_input-path>` or 
:ref:`output <rad-plot-dos_output-path>` path use:

.. code-block:: bash

    rad-plot-dos.py -ip "input_path" -op "output_path" -ew -10 5

.. _rad-plot-dos_arguments:

Arguments
=========

.. _rad-plot-dos_input-path:

-ip, --input-path
-----------------
Relative or absolute path to the folder with PDOS files.

.. code-block:: text

    default : current directory (".")


.. _rad-plot-dos_seedname:

-s, --seedname
--------------
Prefix for input files with PDOS(E). 

In the case of Quantum Espresso-produced pdos it is the same
as specified in the QE projwfc.x input file (filpdos).

If it is not provided the script tries to 
detect it automatically in the 
:ref:`--input-path <rad-plot-dos_input-path>` folder.

.. code-block:: text

    default : None

.. versionchanged:: 0.5.21 from "filpdos" to "seedname".


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

By default the whole energy range present in the files is plotted.

.. code-block:: text

    default : None

Renamed in version 0.5.21: from "window" to "energy-window".


.. _rad-plot-dos_dos-window:

-dw, --dos-window
-----------------
DOS window for the plots. 

By default the whole energy range present in the 
:ref:`--energy-window <rad-plot-dos_energy-window>` is plotted.

.. code-block:: text

    default : None

.. versionadded:: 0.5.21


.. _rad-plot-dos_efermi:

-ef, --efermi
-------------
Fermi energy. 

Zero is shifted to Fermi energy.

.. code-block:: text

    default : 0


.. _rad-plot-dos_separate:

-sep, --separate
----------------
Whether to plot projected DOS for each atom of the same type separately.

.. code-block:: text

    default : False


.. _rad-plot-dos_relative:

-r, --relative
--------------
Whether to use relative style.

.. code-block:: text

    default : False


.. _rad-plot-dos_normalize:

-n, --normalize
---------------
Whether to normalized PDOS values to 1.

(with respect to LDOS of each plot or to total PDOS if
:ref:`--background-total <rad-plot-dos_background-total>`).

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
Whether to save figures as .pickle files.

.. code-block:: text

    default : False

.. versionadded:: 0.5.21


.. _rad-plot-dos_save-txt:

-st, --save-txt
---------------
Whether to save the data as txt files.

.. note::
    It does not affect "pdos-vs-dos.png", 
    because these data are accessible directly from PDOS input files.

.. code-block:: text

    default : False

.. versionadded:: 0.5.21


.. _rad-plot-dos_background-total:

-bt, --background-total
-----------------------
Whether to use total PDOS as the background for all plots.

Total partial density of states is used instead of corresponding 
local density of states in all background data 
(and in all normalization routines as well) .

.. code-block:: text

    default : False

.. versionadded:: 0.5.21

.. _rad-plot-dos:

***************
rad-plot-dos.py
***************

.. versionchanged:: 0.5.21 Renamed from ``rad-dos-plotter.py``

Script for visualization of projected density of states.
It supports the result of the calculation from:

* |QE|_ (|projwfc|_).
    All type of calculations are supported 
    (Collinear, spin-unpolarized; 
    collinear, spin-polarized; 
    non-collinear, non-spin-orbit; 
    non-collinear, spin-orbit)

When summation over the atoms of the same type is involved PDOS is just summed.
In contrary, when the summation over k-points is involved, 
PDOS is summed and divided by the number of k points.

K-resolved PDOS are handled by summing over kpoints (and dividing by the number of kpoints).

There are two ways to use this script:

* :ref:`rad-plot-dos_predefined-plots`
* :ref:`rad-plot-dos_custom-plots`

Before we discuss the difference between those two types, we describe 
parameters, which affects both of them

Common remarks
==============

Script tries to detect all seednames present 
in :ref:`rad-plot-dos_input-path` by default, 
unless :ref:`rad-plot-dos_seedname` option specifies 
particular seedname to work with.

Plot styles
-----------

There are two main plot style:

* Default
    PDOS of each projector appears on the separate plot within the same figure.


.. figure:: /../examples/rad-plot-dos/style-examples/csp/custom.png
    :align: center
    :target: ../../../../../_images/custom.png

    Default style

* :ref:`Relative <rad-plot-dos_relative>`
    PDOS of all projectors appear on the same plot, 
    PDOS of each projector starts at the end of the previous projector`s PDOS.


.. figure:: /../examples/rad-plot-dos/style-examples/csp-relative/custom.png
    :align: center
    :target: ../../../../../_images/custom1.png

    Relative style

Each of these styles could be modified by the following "substyles":

* :ref:`Normalized <rad-plot-dos_normalize>`
    PDOS of each projector is normalized with respect to the local DOS. 
    Local DOS could be the sum of all PDOS, as well as total PDOS
    (see :ref:`rad-plot-dos_background-total`).

.. figure:: /../examples/rad-plot-dos/style-examples/csp-normalized/custom.png
    :align: center
    :target: ../../../../../_images/custom2.png

    With defautlt style

.. figure:: /../examples/rad-plot-dos/style-examples/csp-relative-normalized/custom.png
    :align: center
    :target: ../../../../../_images/custom3.png

    With relative style

* :ref:`Total as a background <rad-plot-dos_background-total>`
    Total PDOS is used as the background values instead 
    of the sum of the PDOS from the plot (which is used by default).

.. figure:: /../examples/rad-plot-dos/style-examples/csp-vstotal/custom.png
    :align: center
    :target: ../../../../../_images/custom4.png

    With defautlt style

.. figure:: /../examples/rad-plot-dos/style-examples/csp-relative-vstotal/custom.png
    :align: center
    :target: ../../../../../_images/custom5.png

    With relative style

* Both of them

.. figure:: /../examples/rad-plot-dos/style-examples/csp-normalized-vstotal/custom.png
    :align: center
    :target: ../../../../../_images/custom6.png

    With defautlt style

.. figure:: /../examples/rad-plot-dos/style-examples/csp-relative-normalized-vstotal/custom.png
    :align: center
    :target: ../../../../../_images/custom7.png

    With relative style

In addition one could modify the colours used in the 
:ref:`relative <rad-plot-dos_relative>` or :ref:`custom <rad-plot-dos_custom-plots>` plots
with the :ref:`rad-plot-dos_colours` parameter.

:ref:`rad-plot-dos_efermi` allows to shift zero to the value of Fermi energy.

:ref:`rad-plot-dos_energy-window` and :ref:`rad-plot-dos_dos-window` allows to specify 
energy and states/eV windows.

Interactive plot
----------------

:ref:`rad-plot-dos_interactive` opens each plot in an interactive matplotlib window.
It allows one to modify the range and appearance of the plot (to drag the legend).

Output remarks
--------------

For each seedname a separate folder "seedname-suffix" is created.

Suffix is a combination of any number of the following words:

* "separate" - appears if :ref:`rad-plot-dos_separate` option is used.
* "relative" - appears if :ref:`rad-plot-dos_relative` option is used.
* "normalized" - appears if :ref:`rad-plot-dos_normalize` option is used.
* "vstotal" - appears if :ref:`rad-plot-dos_background-total` option is used.

.. note::
    :ref:`rad-plot-dos_separate` option contribute to the suffix in both cases, 
    but affects only the :ref:`rad-plot-dos_predefined-plots`.

The structure of the :ref:`output folder <rad-plot-dos_output-path>` is the following:

.. code-block:: text

    output_path/
    ├── ....
    ├── seedname_1-suffixes/
    ├── ...
    └── seedname_n-suffixes/

By default only the pictures (.png) are created. Two additional formats of the output are:

Output formats
--------------

* txt (:ref:`rad-plot-dos_save-txt`)
    Content of the plots in txt format. First line is the header with projectors. 
    It has the same name as the corresponding picture, but the extension is ".txt".
* pickle (:ref:`rad-plot-dos_save-pickle`)
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


.. _rad-plot-dos_custom-plots:

Custom plots
============

Custom plots allows the user to create plots with the hand-picked PDOS.

In order to get custom plot one have to provide :ref:`rad-plot-dos_custom` argument.

As a parameters this argument requires any number of strings, 
where each string specifies one PDOS for the plot. 
This string specifies the set of atoms and projectors, which are summed to produce PDOS.

.. note::
    In reality projector specify the set of projectors: :math:`p` but not :math:`p_x`

The following rules apply to the construction of the input string:

.. role:: color1
.. role:: color2
.. role:: color3
.. role:: color4

* :color1:`atom_type` is required
    Each string can correspond only to one atom type.
* Atom numbers (:color2:`n1` and :color2:`n2`) are optional.
    Each atom number is preceded by exactly one "#" symbol. 
    If no numbers are provided, then the sum is carried out over all atoms of the type
    :color1:`atom_type`.
* Projectors section is optional.
    Projector section is enclosed in parenthesis. 
    It is either absent or contains at least one :color3:`projector_type`. 
    If projectors are not specified, then the sum is carried out 
    over all projectors for each atom.
* Projectors are separated by commas.
    Each comma has to be followed by the projector.
* Projector numbers (:color4:`m1`, :color4:`m2`, :color4:`k2`, :color4:`k2`) are optional.
    Each projector number is preceded by exactly one "#" symbol.
    If no numbers are provided for :color3:`projector_type`, then the sum is carried out
    over all projectors of the type :color3:`projector_type` for each atom.
* Spaces are ignored.
    Feel free to add as many space as you wish. Keep in mind that input string serves 
    as a label in the plot as is.

The format of the string:

:color1:`atom_type`\#\ :color2:`n1`\#\ :color2:`n2`... 
(:color3:`projector_type1`\#\ :color4:`m1`\#\ :color4:`m2`, 
:color3:`projector_type2`:color4:`\#\ k1`\#\ :color4:`k2`, ...)

Here is an example of the set of PDOS file from |projwfc|_ output:

#. seedname.pdos_atm\#\ :color2:`1`\(\ :color1:`Ni`)_wfc\#\ :color4:`1`\(\ :color3:`s`)

#. seedname.pdos_atm\#\ :color2:`1`\(\ :color1:`Ni`)_wfc\#\ :color4:`2`\(\ :color3:`p`)

#. seedname.pdos_atm\#\ :color2:`1`\(\ :color1:`Ni`)_wfc\#\ :color4:`3`\(\ :color3:`d`)

#. seedname.pdos_atm\#\ :color2:`1`\(\ :color1:`Ni`)_wfc\#\ :color4:`4`\(\ :color3:`s`)

#. seedname.pdos_atm\#\ :color2:`2`\(\ :color1:`I`)_wfc\#\ :color4:`1`\(\ :color3:`s`)

#. seedname.pdos_atm\#\ :color2:`2`\(\ :color1:`I`)_wfc\#\ :color4:`2`\(\ :color3:`p`)

#. seedname.pdos_atm\#\ :color2:`3`\(\ :color1:`I`)_wfc\#\ :color4:`1`\(\ :color3:`s`)

#. seedname.pdos_atm\#\ :color2:`3`\(\ :color1:`I`)_wfc\#\ :color4:`2`\(\ :color3:`p`)

Where the colour code specify the correspondence of the input string parts to the
|projwfc|_ output files. 

Here are few examples of the input strings:

* ":color1:`Ni`" 
    Sums over all projectors of Ni: 1-4. 
    Equivalent to: ":color1:`Ni`\#\ :color2:`1"` or 
    ":color1:`Ni` \(\ :color3:`s`, :color3:`p`, :color3:`d`)" or 
    ":color1:`N   i`"
* ":color1:`Ni` \(\ :color3:`s`\#\ :color4:`1`, :color3:`d`)"
    Sums over one s and d projector of Ni: 1, 3
* ":color1:`I`"
    Sums over all projectors of all I atoms: 5-8. 
    Equivalent to ":color1:`I`\#\ :color2:`2`\#\ :color2:`3` \(\ :color3:`s`, :color3:`p`)"
* ":color1:`I`\#\ :color2:`3` \(\ :color3:`p`)"
    Sums over p projector of the second I atom: 8.
* ":color1:`I` \(\ :color3:`p`\#\ :color4:`2`)"
    Sums over p projector of all I atoms: 6, 8. 
    Equivalent to ":color1:`I` \(\ :color3:`p`)" or 
    ":color1:`I`\#\ :color2:`2`\#\ :color2:`3` \(\ :color3:`p`)"

Output file of the custom plot is located in the output folder with the name "custom.png"
(with corresponding txt or pickle output if any). 
If "custom.png" already exists in the output folder, 
then integer number is added to the end ("custom1.txt") 
in order to prevent accidental loss of the previous files. 
Integer is the smallest one, which provides unique name.


.. _rad-plot-dos_predefined-plots:

Predefined plots
================

The predefined plots are:

* pdos-vs-dos.png
    Total partial density of states (sum over all projectors) vs 
    total density of states (directly from the plane-wave basis).
    Affected by :ref:`rad-plot-dos_save-pickle`.
* atomic-contributions.png
    Contribution of each atom (summed over all projectors) 
    into the total partial density of states.
    Affected by :ref:`rad-plot-dos_save-pickle`.
    Affected by :ref:`rad-plot-dos_save-txt`.
* atom-resolved/
    Contribution of each projectors group (i.e. :math:`s`, :math:`p`, :math:`d`, :math:`f`) 
    into the partial density of state of each atom.
    Affected by :ref:`rad-plot-dos_save-pickle`.
    Affected by :ref:`rad-plot-dos_save-txt`.
* orbital-resolved/
    Contribution of each projector (i. e. :math:`p_z`, :math:`p_x`, :math:`p_y`) into the total 
    partial density of states of each group (i.e. :math:`p`).
    Affected by :ref:`rad-plot-dos_save-pickle`.
    Affected by :ref:`rad-plot-dos_save-txt`.

Option :ref:`rad-plot-dos_separate` plots PDOS individually for each atom.

Each seedname folder has the structure:

.. code-block:: text

    seedname/
    ├── pdos-vs-dos.png
    ├── atomic-contributions.png
    ├── atom-resolved/
    │   ├── output_name_1
    │   └── output_name_2
    └── orbital-resolved/
        ├── output_name_1
        └── output_name_2


Usage example
=============
Minimal possible input is:

.. code-block:: bash

    rad-plot-dos.py 

It will try to detect PDOS output files in the current directory and plot them.

To choose energy window use an 
option :ref:`rad-plot-dos_energy-window`:

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
:ref:`rad-plot-dos_input-path` folder.

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

By default the whole states/eV range present in the 
:ref:`rad-plot-dos_energy-window` is plotted.

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
:ref:`rad-plot-dos_background-total` is used).

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


.. _rad-plot-dos_custom:

--custom
--------
Custom PDOS plot. See :ref:`rad-plot-dos_custom-plots` for info.

.. code-block:: text

    default : None
    nargs : any

.. versionadded:: 0.7.5


.. _rad-plot-dos_colours:

-cls, --colours
---------------
Colours for the relative and custom plots.

Values are passed directly to the matplotlib as strings, 
therefore any valid value is allowed. Examples: "red" or "#FF0000".
When :ref:`rad-plot-dos_custom` is used the order of colours is the same as for 
the values of the :ref:`rad-plot-dos_custom`.

.. code-block:: text

    default : None
    nargs : any

.. versionadded:: 0.7.5

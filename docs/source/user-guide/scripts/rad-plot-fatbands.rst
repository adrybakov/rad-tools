.. _rad-plot-fatbands:

********************
rad-plot-fatbands.py
********************


Script for visualization of projected density of states.
It supports the result of the calculation from:

* |QE|_ (|projwfc|_).
    All type of calculations are supported 
    (Collinear, spin-unpolarized; 
    collinear, spin-polarized; 
    non-collinear, non-spin-orbit; 
    non-collinear, spin-orbit)

.. warning::

    Preliminary version, work in progress.


.. _rad-plot-fatbands_arguments:

Arguments
=========

.. _rad-plot-fatbands_input-folder:

-if, --input-folder
-------------------
Relative or absolute path to the folder with PDOS files.

.. code-block:: text

    default : current directory (".")


.. _rad-plot-fatbands_seedname:

-s, --seedname
--------------
Prefix for input files with PDOS(E). 

In the case of Quantum Espresso-produced pdos it is the same
as specified in the QE projwfc.x input file (filpdos).

If it is not provided the script tries to 
detect it automatically in the 
:ref:`rad-plot-fatbands_input-folder` folder.

.. code-block:: text

    default : None



.. _rad-plot-fatbands_output-name:

-on, --output-name
------------------
Relative or absolute name to the folder for saving outputs.

.. code-block:: text

    default : current directory (".")

See also :ref:`output-notes`.


.. _rad-plot-fatbands_energy-window:

-ew, --energy-window
--------------------
Energy window for the plots.  

By default the whole energy range present in the files is plotted.

.. code-block:: text

    default : None



.. _rad-plot-fatbands_k-window:

-kw, --k-window
---------------
K-point window for the plots.

By default whole range present in the files is plotted.

.. code-block:: text

    default : None
    nargs : 2


.. _rad-plot-fatbands_efermi:

-ef, --efermi
-------------
Fermi energy. 

Zero is shifted to Fermi energy.

.. code-block:: text

    default : 0


.. _rad-plot-fatbands_verbose:

-v, --verbose
-------------
Verbose output, propagates to the called methods.

.. code-block:: text

    default : False


.. _rad-plot-fatbands_interactive:

-i, --interactive
-----------------
Interactive plotting.

.. code-block:: text

    default : False

.. _rad-plot-fatbands_separate:

-sep, --separate
----------------
Whether to plot projected DOS for each atom of the same type separately.

.. code-block:: text

    default : False


.. _rad-plot-fatbands_save-pickle:

-sp, --save-pickle
------------------
Whether to save figures as .pickle files.

.. code-block:: text

    default : False


.. _rad-plot-fatbands_save-txt:

-st, --save-txt
---------------
Whether to save the data as txt files.

.. note::
    It does not affect "pdos-vs-dos.png", 
    because these data are accessible directly from PDOS input files.

.. code-block:: text

    default : False


.. _rad-plot-fatbands_custom:

--custom
--------
Custom PDOS plot. 

.. code-block:: text

    default : None
    nargs : any



.. _rad-plot-fatbands_colours:

-cls, --colours
---------------
Colours for the relative and custom plots.

Values are passed directly to the matplotlib as strings, 
therefore any valid value is allowed. Examples: "red" or "#FF0000".
When :ref:`rad-plot-fatbands_custom` is used the order of colours is the same as for 
the values of the :ref:`rad-plot-fatbands_custom`.

.. code-block:: text

    default : None
    nargs : any



.. _rad-plot-fatbands_labels:

-lbs, --labels
--------------
Labels for the custom plots.

Amount of labels have to be the same as the amount of custom strings, or one more.
If one more, then first one is interpreted as the label for the background 
(Use "None" to switch it off). If the amount of argument is one more  and the first one is None, 
then the label for the total PDOS is switched off and the total PDOS itself is not plotted.


.. code-block:: text

    default : None
    nargs : any



.. _rad-plot-fatbands_axes-labels-fontsize:

-alfs, --axes-labels-fontsize
-----------------------------
Fontsize of the labes of the axes.

.. code-block:: text

    default : 14
    type : int



.. _rad-plot-fatbands_legend-fontsize:

-lfs, --legend-fontsize
-----------------------
Fontsize of the legend.

.. code-block:: text

    default : 12
    type : int



.. _rad-plot-fatbands_title-fontsize:

-tfs, --title-fontsize
----------------------
Fontsize of the title.

.. code-block:: text

    default : 18
    type : int

.. _rad-plot-fatbands_k-points:

-kp, --k-points
---------------
List of high symmetry points.

.. code-block:: text

    default : None



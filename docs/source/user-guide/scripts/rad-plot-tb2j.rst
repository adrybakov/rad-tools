.. _rad-plot-tb2j:

****************
rad-plot-tb2j.py
****************

Script for visualizations of 
|TB2J|_  results.

.. versionchanged:: 0.6 Renamed from ``tb2j-plotter.py``

The script displays isotropic exchange, distances and DMI 
(one output file for each). 

Supports filtering by 
R vectors (see :ref:`--R-vector <rad-plot-tb2j_R-vector>`), 
distances (see :ref:`--max-distance <rad-plot-tb2j_max-distance>`,
:ref:`--min-distance <rad-plot-tb2j_min-distance>` and
:ref:`--distance <rad-plot-tb2j_distance>`), 
and template file (see :ref:`--template-file <rad-plot-tb2j_template-file>`). 
The result is defined by logical conjugate of the specified conditions.

:ref:`--input-filename <rad-plot-tb2j_input-filename>` 
(or :ref:`-if <rad-plot-tb2j_input-filename>`) argument is required, 
the rest of them are optional.


Output files have the following name structure: 
"output-name.display-data-type.png" 

.. _rad-plot-tb2j_example:

Usage example
=============

Example is based on the exchange.out file from 
:examples:`examples folder <rad-plot-tb2j>`. 

There is one required argument in the script 
(:ref:`--input-filename <rad-plot-tb2j_input-filename>`), therefore a minimum input 
for the script to run is:

.. code-block:: bash

    rad-plot-tb2j.py -if exchange.out

which produces three pictures "exchange.iso.png", 
"exchange.distance.png", "exchange.distance.png". each file name has a shared
default seedname ("exchange", use 
:ref:`--output-name <rad-plot-tb2j_output-name>` to change it), data type 
("iso", "dmi", "distance", see 
:ref:`--what-to-plot <rad-plot-tb2j_what-to-plot>`) and file extension ".png".

.. dropdown:: Output images

    .. figure:: /../examples/rad-plot-tb2j/exchange.iso.png
        :align: center

        exchange.iso.png

    .. figure:: /../examples/rad-plot-tb2j/exchange.dmi.png
        :align: center

        exchange.dmi.png

    .. figure:: /../examples/rad-plot-tb2j/exchange.distance.png
        :align: center

        exchange.distance.png

.. note::

    In the following text only "exchange.iso.png" file is produced with the help 
    of  :ref:`--what-to-plot <rad-plot-tb2j_what-to-plot>` argument.

Basic adjustments
-----------------

Since "the exchange.out" file contains a lot of exchange bonds the pictures with 
all of them are not really useful. Lets plot the isotropic exchange picture 
with some adjustments:

    * Filter the model by maximum distance (:ref:`-md <rad-plot-tb2j_max-distance>`).
    * Draw unit cell borders (:ref:`-dc <rad-plot-tb2j_draw-cells>`)
    * Scale the marks of atoms (:ref:`-sa <rad-plot-tb2j_scale-atoms>`)
    * Scale the data (:ref:`-sd <rad-plot-tb2j_scale-data>`)
    * Add the title to the plot (:ref:`-t <rad-plot-tb2j_title>`)

.. code-block:: bash

    rad-plot-tb2j.py -if exchange.out -wtp iso -maxd 5 -dc -sa 1.2 -sd 1.2 -t "First neighbour exchange" -on exchange_filtered

.. figure:: /../examples/rad-plot-tb2j/exchange_filtered.iso.png
    :align: center

    exchange_filtered.iso.png

Filtering
---------

For filtering the exchange model there are a few options available:

    * :ref:`--max_distance <rad-plot-tb2j_max-distance>`
    * :ref:`--min_distance <rad-plot-tb2j_min-distance>`
    * :ref:`--distance <rad-plot-tb2j_distance>`
    * :ref:`--R-vector <rad-plot-tb2j_R-vector>`
    * :ref:`--template <rad-plot-tb2j_template-file>`

Here is an example of how to filter exchange model in order to show 
first exchange neighbour, using different options:
:ref:`--max_distance <rad-plot-tb2j_max-distance>`
:ref:`--R-vector <rad-plot-tb2j_R-vector>` or
:ref:`--template <rad-plot-tb2j_template-file>` arguments:

.. code-block:: bash

    rad-plot-tb2j.py -if exchange.out -wtp iso -maxd 5 -dc -sa 1.2 u-sd 1.2 -t "First neighbour exchange" -on exchange_filtered
    rad-plot-tb2j.py -if exchange.out -wtp iso -tf template.txt -dc -sa 1.2 -sd 1.2 -t "First neighbour exchange" -on exchange_template
    rad-plot-tb2j.py -if exchange.out -wtp iso -R 1 0 0 1 1 0  0 1 0 -1 1 0 -1 0 0 -1 -1 0 0 -1 0 1 -1 0 -dc -sa 1.2 -sd 1.2 -t "First neighbour exchange" -on exchange_R

where template file is the following:

.. literalinclude:: /../examples/rad-plot-tb2j/template.txt
    :language: text

The images should be the same:

.. dropdown:: Output images

    .. figure:: /../examples/rad-plot-tb2j/exchange_filtered.iso.png
        :align: center

        exchange_filtered.iso.png

    .. figure:: /../examples/rad-plot-tb2j/exchange_R.iso.png
        :align: center

        exchange_R.iso.png

    .. figure:: /../examples/rad-plot-tb2j/exchange_template.iso.png
        :align: center

        exchange_template.iso.png

Modifying the model
-------------------

By default ``rad-plot-tb2j.py`` displays the bonds as it is in the model.
:ref:`-fs <rad-plot-tb2j_force-symmetry>` argument helps
to reproduce particular exchange model:

.. code-block:: bash

    rad-plot-tb2j.py -if exchange.out -tf template.txt -fs -dc -sa 1.2 -sd 1.2 -t "Forced symmetry exchange" -on exchange_forced_symmetry

.. figure:: /../examples/rad-plot-tb2j/exchange_forced_symmetry.iso.png
    :align: center

    exchange_forced_symmetry.iso.png

.. dropdown:: DMI and distances

    .. figure:: /../examples/rad-plot-tb2j/exchange_forced_symmetry.dmi.png
        :align: center

        exchange_forced_symmetry.dmi.png

    .. figure:: /../examples/rad-plot-tb2j/exchange_forced_symmetry.distance.png
        :align: center

        exchange_forced_symmetry.distance.png

Only one exchange parameter is present in the template file, 
therefore the model is filtered with respect to the template 
and then the value of the exchange for each bond is set to 
the medium value of all bonds from the same exchange group.
The direction of the DMI vectors is kept, but the magnitude 
of the DMI vector is scaled to the medium value.

.. note::

    When :ref:`--force-symmetry <rad-plot-tb2j_force-symmetry>` argument is provided
    :ref:`--template-file <rad-plot-tb2j_template-file>` is required.

.. _rad-plot-tb2j_arguments:

Arguments
=========

.. _rad-plot-tb2j_input-filename:

-if, --input-filename
---------------------
Relative or absolute path to the "exchange.out" file, 
including name and extension of the file.

.. code-block:: text

    required


.. _rad-plot-tb2j_output-path:

-op, --output-path
------------------
Relative or absolute path to the folder for saving outputs.

If the folder does not exist then it is created from the specified path.
The creation is applied recursively to the path, starting from the right
until the existing folder is reached.

.. code-block:: text

    default : current directory

See also: :ref:`example <output-notes>`.


.. _rad-plot-tb2j_output-name:

-on, --output-name
------------------
Seedname for the output files.

Output files have the following name structure:
"output-name.display-data-type.png"

.. code-block:: text

    default : exchange

See also: :ref:`example <output-notes>`.


.. _rad-plot-tb2j_what-to-plot:

-wtp, --what-to-plot
--------------------
Type of data for display.

Specifying the data to be displayed in the graphs. 
Everything is displayed by default, each value in a separate picture. 
Currently available for display: Isotropic exchange parameter, distance, \|DMI\|.

.. code-block:: text

    default : all


.. _rad-plot-tb2j_draw-cells:

-dc, --draw-cells
-----------------
Whether to draw the cells.

If specified then the shapes of all cells 
presented in the model (after filtering) are drawn. (0, 0, 0) is red.

.. code-block:: text

    default : False


.. _rad-plot-tb2j_R-vector:

-R, --R-vector
--------------
R vectors for filtering the model.

In TB2J outputs the bond is defined by atom 1 (from) and atom 2 (to). 
Atom 1 is always located in (0, 0, 0) unit cell, while atom 2 is located in 
R = (i, j, k) unit cell. This parameter tells the script to keep only the 
bonds for which atom 2 is located in one of specified R supercells. 
Supercells are specified by a set of integers separated by spaces. 
They are grouped by three starting from the left and forms a set 
of R vectors. If the last group contains 1 or 2 integers they are ignored.

.. code-block:: text

    default : None


.. _rad-plot-tb2j_max-distance:

-maxd, --max-distance
---------------------
(<=) Maximum distance.

All the bonds with the distance between atom 1 and atom 2 
greater than maximum distance are excluded from the model.

.. code-block:: text

    default : None


.. _rad-plot-tb2j_min-distance:

-mind, --min-distance
---------------------
(>=) Minimum distance.

All the bonds with the distance between atom 1 and atom 2 
lower than minimum distance are excluded from the model.

.. code-block:: text

    default : None


.. _rad-plot-tb2j_distance:

-d, --distance
--------------
(=) Exact distance.

Only the bonds with the exact distance remains in the model.

.. code-block:: text

    default : None

.. hint::
    There is no point in specifying maximum or minimum distance when 
    this parameter is provided.


.. _rad-plot-tb2j_template-file:

-tf, --template-file
--------------------
Relative or absolute path to the template file, 
including the name and extension of the file.

.. code-block:: text

    default : None

See also: :ref:`template <rad-make-template>`


.. _rad-plot-tb2j_scale-atoms:

-sa, --scale-atoms
------------------
Scale for the size of atom marks.

Use it if you want to make atom marks bigger (>1) or smaller (<1). 
Has to be positive.

.. code-block:: text

    default : 1


.. _rad-plot-tb2j_scale-data:

-sd, --scale-data
-----------------
Scale for the size of data text.

Use it if you want to make data text marks bigger (>1) or smaller (<1). 
Has to be positive.

.. code-block:: text

    default : 1


.. _rad-plot-tb2j_title:

-t, --title
-----------
Title for the plots.

Title is displayed in the picture.

.. code-block:: text

    default : None


.. _rad-plot-tb2j_force-symmetry:

-fs, --force-symmetry
---------------------
Force the exchange model to have the symmetry of the template.

.. code-block:: text

    default : False


.. _rad-plot-tb2j_verbose:

-v, --verbose
-------------
Verbose output, propagates to the called methods.

.. code-block:: text

    default : False
 
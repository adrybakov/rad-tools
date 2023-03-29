.. _tb2j-plotter:

*******************
``tb2j-plotter.py``
*******************

Script for visualisation of 
`TB2J <https://tb2j.readthedocs.io/en/latest/>`_ results.

The script displays isotropic exchange, distances and DMI 
(one output file for each). 

Supports filtering by 
R vectors (see :ref:`--R-vector <tb2j-plotter_R-vector>`), 
distances (see :ref:`--max-distance <tb2j-plotter_max-distance>`,
:ref:`--min-distance <tb2j-plotter_min-distance>` and
:ref:`--distance <tb2j-plotter_distance>`), 
and template file (see :ref:`--template-file <tb2j-plotter_template-file>`). 
The result is defined by logical conjugate of the specified conditions.

``--input-filename`` (or ``-if``) argument is required, the rest of them are optional.


Output files will have the following name structure: 
*output-name.display-data-type.png* 

.. _tb2j-plotter_example:

Usage example
=============

Example is based on the exchange.out file from 
:examples:`examples folder <tb2j-plotter>`. 

There is one required argument in the script 
(:ref:`--input-filename <tb2j-plotter_input-filename>`), so there is a minimum input 
for the script to run:

.. code-block:: bash

    tb2j-plotter.py -if exchange.out

which will produce three pictures *exchange.iso.png*, 
*exchange.distance.png*, *exchange.distance.png*. each file name have a common 
default seedname (*exchange*, use 
:ref:`--output-name <tb2j-plotter_output-name>` to change it), data type 
(*iso*, *dmi*, *distance*, see 
:ref:`--what-to-plot <tb2j-plotter_what-to-plot>`) and file extension.

.. dropdown:: Output images

    .. figure:: /../examples/tb2j-plotter/exchange.iso.png
        :align: center

        exchange.iso.png

    .. figure:: /../examples/tb2j-plotter/exchange.dmi.png
        :align: center

        exchange.dmi.png

    .. figure:: /../examples/tb2j-plotter/exchange.distance.png
        :align: center

        exchange.distance.png

.. note::

    In the following text only exchange.iso.png file will be produced with the help 
    of  :ref:`--what-to-plot <tb2j-plotter_what-to-plot>` argument.

Basic adjustments
-----------------

Since the exchange.out file contains a lot of exchange bonds the pictures with 
all of them are not really useful. Lets plot the isotropic exchange picture 
with some adjustments:

    * Filter the model by maximum distance (:ref:`-md <tb2j-plotter_max-distance>`).
    * Draw unit cell borders (:ref:`-dc <tb2j-plotter_draw-cells>`)
    * Scale the marks of atoms (:ref:`-sa <tb2j-plotter_scale-atoms>`)
    * Scale the data (:ref:`-sd <tb2j-plotter_scale-data>`)
    * Add the title to the plot (:ref:`-t <tb2j-plotter_title>`)

.. code-block:: bash

    tb2j-plotter.py -if exchange.out -wtp iso -maxd 5 -dc -sa 1.2 -sd 1.2 -t "First neighbour exchange" -on exchange_filtered

.. figure:: /../examples/tb2j-plotter/exchange_filtered.iso.png
    :align: center

    exchange_filtered.iso.png

Filtering
---------

For filtering the exchange model there are a few options available:

    * :ref:`--max_distance <tb2j-plotter_max-distance>`
    * :ref:`--min_distance <tb2j-plotter_min-distance>`
    * :ref:`--distance <tb2j-plotter_distance>`
    * :ref:`--R-vector <tb2j-plotter_R-vector>`
    * :ref:`--template <tb2j-plotter_template-file>`

Here is an example of how to filter exchange model in order to show 
first exchange neighbour, using different options:
:ref:`--max_distance <tb2j-plotter_max-distance>`
:ref:`--R-vector <tb2j-plotter_R-vector>` or
:ref:`--template <tb2j-plotter_template-file>` arguments:

.. code-block:: bash

    tb2j-plotter.py -if exchange.out -wtp iso -maxd 5 -dc -sa 1.2 u-sd 1.2 -t "First neighbour exchange" -on exchange_filtered
    tb2j-plotter.py -if exchange.out -wtp iso -tf template.txt -dc -sa 1.2 -sd 1.2 -t "First neighbour exchange" -on exchange_template
    tb2j-plotter.py -if exchange.out -wtp iso -R 1 0 0 1 1 0  0 1 0 -1 1 0 -1 0 0 -1 -1 0 0 -1 0 1 -1 0 -dc -sa 1.2 -sd 1.2 -t "First neighbour exchange" -on exchange_R

where template file is the following:

.. literalinclude:: /../examples/tb2j-plotter/template.txt
    :language: text

The images are the same:

.. dropdown:: Output images

    .. figure:: /../examples/tb2j-plotter/exchange_filtered.iso.png
        :align: center

        exchange_filtered.iso.png

    .. figure:: /../examples/tb2j-plotter/exchange_R.iso.png
        :align: center

        exchange_R.iso.png

    .. figure:: /../examples/tb2j-plotter/exchange_template.iso.png
        :align: center

        exchange_template.iso.png

Modifying the model
-------------------

By default ``tb2j-plotter.py`` display the bonds as it is in the model.
One could use :ref:`-fs <tb2j-plotter_force-symmetry>` argument in order 
to reproduce particular exchange model:

.. code-block:: bash

    tb2j-plotter.py -if exchange.out -tf template.txt -fs -dc -sa 1.2 -sd 1.2 -t "Forced symmetry exchange" -on exchange_forced_symmetry

.. figure:: /../examples/tb2j-plotter/exchange_forced_symmetry.iso.png
    :align: center

    exchange_forced_symmetry.iso.png

.. dropdown:: DMI and distances

    .. figure:: /../examples/tb2j-plotter/exchange_forced_symmetry.dmi.png
        :align: center

        exchange_forced_symmetry.dmi.png

    .. figure:: /../examples/tb2j-plotter/exchange_forced_symmetry.distance.png
        :align: center

        exchange_forced_symmetry.distance.png

Only one exchange parameter present in the template file, therefore the model 
is filtered with respect to the template and then the value of the exchange 
for each bond is set to the medium value of all bonds from the same exchange group.
DMI interaction is modified as well in a way that direction of the DMI vector is 
kept, but the magnitude of the DMI vector is scaled to the medium value.

.. note::

    When :ref:`--force-symmetry <tb2j-plotter_force-symmetry>` argument is provided
    :ref:`--template-file <tb2j-plotter_template-file>` is required.



Arguments
=========

.. _tb2j-plotter_input-filename:

-if, --input-filename
---------------------
Relative or absolute path to the TB2J exchange output file, 
including the name and extension of the file.

.. code-block:: text

    required


.. _tb2j-plotter_output-path:

-op, --output-path
------------------
Relative or absolute path to the folder for saving outputs.

If the folder does not exist then it is created from the specified path.
The creation is applied recursively to the path, starting from the right
until the existing folder is reached.

.. code-block:: text

    default : current directory

See also: :ref:`example <scripts_output-notes>`.


.. _tb2j-plotter_output-name:

-on, --output-name
------------------
Seedname for the output files.

Output files will have the following name structure:
*output-name.display-data-type.png*

.. code-block:: text

    default : exchange

See also: :ref:`example <scripts_output-notes>`.


.. _tb2j-plotter_what-to-plot:

-wtp, --what-to-plot
--------------------
Type of data for display.

Specifying the data which will be displayed in the graphs. 
Everything is displayed by default, each value in a separate picture. 
Currently available for display: Isotropic exchange parameter, distance, \|DMI\|.

.. code-block:: text

    default : all


.. _tb2j-plotter_draw-cells:

-dc, --draw-cells
-----------------
Whenever to draw the cells.

If specified then the shape of all cells 
presented in the model (after filtering) is drawn. (0, 0, 0) will be in red.

.. code-block:: text

    default : False


.. _tb2j-plotter_R-vector:

-R, --R-vector
--------------
R vectors for filtering the model.

In TB2J outputs the bond is defined by atom 1 (from) and atom 2 (to). 
Atom 1 is always located in (0, 0, 0) supercell, while atom 2 is located in 
R = (i, j, k) supercell. This parameter tells the script to keep only the 
bonds for which atom 2 is located in one of specified R supercells. 
In order to specify supercells provide a set of integers separated 
by spaces. They are grouped by three starting from the left to form a set 
of R vectors. If the last group will contain 1 or 2 integers they will be 
ignored.

.. code-block:: text

    default : None


.. _tb2j-plotter_max-distance:

-maxd, --max-distance
---------------------
(<=) Maximum distance.

All the bonds with the distance between atom 1 and atom 2 
greater than maximum distance are excluded from the model.

.. code-block:: text

    default : None


.. _tb2j-plotter_min-distance:

-mind, --min-distance
---------------------
(>=) Minimum distance.

All the bonds with the distance between atom 1 and atom 2 
lower than minimum distance are excluded from the model.

.. code-block:: text

    default : None


.. _tb2j-plotter_distance:

-d, --distance
--------------
(=) Exact distance.

Only the bonds with the exact distance remains in the model.

.. code-block:: text

    default : None

.. hint::
    There is no point in specifying maximum or minimum distance when 
    this parameter is provided.


.. _tb2j-plotter_template-file:

-tf, --template-file
--------------------
Relative or absolute path to the template file, 
including the name and extension of the file.

.. code-block:: text

    default : None

See also: :ref:`template <rad-make-template>`


.. _tb2j-plotter_scale-atoms:

-sa, --scale-atoms
------------------
Scale for the size of atom marks.

Use it if you want to display atom marks bigger or smaller. 
Have to be positive.

.. code-block:: text

    default : 1


.. _tb2j-plotter_scale-data:

-sd, --scale-data
-----------------
Scale for the size of data text.

Use it if you want to display data text marks bigger or smaller. 
Have to be positive.

.. code-block:: text

    default : 1


.. _tb2j-plotter_title:

-t, --title
-----------
Title for the plots

Title will be displayed in the picture.

.. code-block:: text

    default : None


.. _tb2j-plotter_force-symmetry:

-fs, --force-symmetry
---------------------
Force the exchange model to have the symmetry of the template.

.. code-block:: text

    default : False


.. _tb2j-plotter_verbose:

-v, -verbose
------------
Verbose output, propagates to the called methods.

.. code-block:: text

    default : False
 
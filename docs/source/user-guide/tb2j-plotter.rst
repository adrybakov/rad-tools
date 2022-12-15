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

``--filename`` (or ``-f``) argument is required, the rest of them are optional.


Output files will have the following name structure: 
*output-name.display-data-type.png* 

.. _tb2j-plotter_example:

Usage example
=============

Example is based on the exchange.out file from 
:examples:`examples folder <tb2j-plotter>`. 

There is one required argument in the script 
(:ref:`--filename <tb2j-plotter_filename>`), so there is a minimum input 
for the script to run:

.. code-block:: bash

    tb2j-plotter.py -f exchange.out

which will produce three pictures *exchange.iso.png*, 
*exchange.distance.png*, *exchange.distance.png*. each file name have a common 
default seedname (*exchange*, use 
:ref:`--output-name <tb2j-plotter_output-name>` to change it), data type 
(*iso*, *dmi*, *distance*, see 
:ref:`--what-to-plot <tb2j-plotter_what-to-plot>`) and file extention.

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

    In the following only exchange.iso.png file will be produced with the help 
    of  :ref:`--what-to-plot <tb2j-plotter_what-to-plot>` argument.

Since the exchange.out file contains a lot of exchange bonds the pictures with 
all of them are not really usefull. Lets plot the isotropic exchange picture 
with some adjustments:

    * Filter the model by maximum distance (:ref:`-md <tb2j-plotter_max-distance>`).
    * Draw unit cell borders (:ref:`-dc <tb2j-plotter_draw-cells>`)
    * Scale the marks of atoms (:ref:`-sa <tb2j-plotter_scale-atoms>`)
    * Scale the data (:ref:`-sd <tb2j-plotter_scale-data>`)
    * Add the title to the plot (:ref:`-t <tb2j-plotter_title>`)

.. code-block:: bash

    tb2j-plotter.py -f exchange.out -wtp iso -maxd 5 -dc -sa 1.2 -sd 1.3 -t "First neighbor exchange" -on exchange_filtered

.. figure:: /../examples/tb2j-plotter/exchange_filtered.iso.png
    :align: center

    exchange_filtered.iso.png

For filtering the exchange model there is a few options available:

    * :ref:`--max_distance <tb2j-plotter_max-distance>`
    * :ref:`--min_distance <tb2j-plotter_min-distance>`
    * :ref:`--distance <tb2j-plotter_distance>`
    * :ref:`--R-vector <tb2j-plotter_R-vector>`
    * :ref:`--template <tb2j-plotter_template-file>`

By default ``tb2j-plotter.py`` does not display symmetrically equivalent
(with respect to the translation symmetry) bonds. It could raise the problems 
with analysis of DMI or anisotropic echange, in order to display all bonds 
one can use :ref:`-db <tb2j-plotter_double-bonds>` argument:

.. code-block:: bash

    tb2j-plotter.py -f exchange.out -wtp iso -maxd 5 -dc -sa 1.2  -t "First neighbor exchange" -db -on exchange_double_bonds

.. figure:: /../examples/tb2j-plotter/exchange_double_bonds.iso.png
    :align: center

    exchange_double_bonds.iso.png

At the end there is an example of how to filter exchange model in order to show 
frist exchange neighbor, but using 
:ref:`--R-vector <tb2j-plotter_R-vector>` or
:ref:`--template <tb2j-plotter_template-file>` arguments:

.. code-block:: bash

    tb2j-plotter.py -f exchange.out -wtp iso -tf template.txt -dc -sa 1.2  -t "First neighbor exchange" -db -on exchange_template
    tb2j-plotter.py -f exchange.out -wtp iso -R 1 0 0 1 1 0  0 1 0 -1 1 0 -1 0 0 -1 -1 0 0 -1 0 1 -1 0 -dc -sa 1.2  -t "First neighbor exchange" -db -on exchange_R

where template file is the following:

.. literalinclude:: /../examples/tb2j-plotter/template.txt
    :language: text

The images should be the same:

.. dropdown:: Output images

    .. figure:: /../examples/tb2j-plotter/exchange_R.iso.png
        :align: center

        exchange_R.iso.png

    .. figure:: /../examples/tb2j-plotter/exchange_template.iso.png
        :align: center

        exchange_template.iso.png



Arguments
=========

.. _tb2j-plotter_filename:

-f, --filename
--------------
Relative or absulute path to the TB2J exchange output file, 
including the name and extention of the file.

    *required* : True

    *type* : str


.. _tb2j-plotter_output-dir:

-op, --output-dir
-----------------
Relative or absolute path to the folder for saving outputs.

If the folder does not exist then it is created from the specified path.
The creation is applied recursevly to the path, starting from the right
until the existing folder is reached.

    *default* : current directory
        
    *type* : str


.. _tb2j-plotter_output-name:

-on, --output-name
------------------
Seedname for the output files.

Output files will have the following name structure:
*output-name.display-data-type.png*

    *default* : exchange
        
    *type* : str

See also: :ref:`example <tb2j-plotter_example>`


.. _tb2j-plotter_what-to-plot:

-wtp, --what-to-plot
--------------------
Type of data for display.

Specifying the data which will be displayed in the graphs. 
Everything is displayed by default, each value in a separate picture. 
Currently available for display: Isotropic exchange parameter, distance, \|DMI\|.

    *default* : all

    *type* : str

    *choices* : all, iso, distance, dmi


.. _tb2j-plotter_draw-cells:

-dc, --draw-cells
-----------------
Whenever to draw the cells.

If specified then the shape of all cells 
presented in the model (after filtering) is drawn. (0, 0, 0) will be in red.

    *default* : False

    *action* : store_true


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

    *default* : None

    *type* : int

    *nargs* : *


.. _tb2j-plotter_max-distance:

-maxd, --max-distance
---------------------
(<=) Maximum distance.

All the bonds with the distance beetwen atom 1 and atom 2 
greater than maximum distance are excluded from the model.

    *default* : None

    *type* : float


.. _tb2j-plotter_min-distance:

-mind, --min-distance
---------------------
(>=) Minimum distance.

All the bonds with the distance beetwen atom 1 and atom 2 
lower than minimum distance are excluded from the model.

    *default* : None

    *type* : float


.. _tb2j-plotter_distance:

-d, --distance
--------------
(=) Exact distance.

Only the bonds with the exact distance remains in the model.

    *default* : None

    *type* : float

.. hint::
    There is no point in specifying maximum or minimum distance when 
    this parameter is provided.


.. _tb2j-plotter_template-file:

-tf, --template-file
--------------------
Relative or absolute path to the template file, 
including the name and extention of the file.

    *default* : None

    *type* : str

See also: :ref:`template <rad-make-template>`


.. _tb2j-plotter_double-bonds:

-db, --double-bonds
-------------------
Whenever to keep both bonds.

In TB2J file there are two bonds for the pair of atom 1 and atom 2: 
from 1 to 2 and from 2 to 1 (when R = (0, 0, 0)). Isotropic and 
anisotropic exchange and distance usially are exactly the same. 
DMI vector have the same module and opposite directions. 
If this parameter is specifyied then both bonds are displayed. 
Otherwise bonds are combined in one by taking the average beetween
exchange parameters. 

    *default* : False

    *action* : store_true

.. caution::
    If this parameter is not specified then it is highly probable that
    DMI will be equal to zero even if it is not zero in TB2J file. 
    Moreover, it is necessary to check anisotropy matrices as well.


.. _tb2j-plotter_scale-atoms:

-sa, --scale-atoms
------------------
Scale for the size of atom marks.

Use it if you want to display atom marks bigger or smaller. 
Have to be positive.

    *default* : 1

    *type* : float


.. _tb2j-plotter_scale-data:

-sd, --scale-data
-----------------
Scale for the size of data text.

Use it if you want to display data text marks bigger or smaller. 
Have to be positive.

    *default* : 1

    *type* : float


.. _tb2j-plotter_title:

-t, --title
-----------
Title for the plots

Title will be displayed in the picture.

    *default* : None

    *type* : str
 
.. _tb2j_plotter:

*******************
``tb2j_plotter.py``
*******************

Script for visualisation of 
`TB2J <https://tb2j.readthedocs.io/en/latest/>`_ results.

Display Isotropic exchange and distances (one output file for each). 

Supports filtering by 
R vectors (see :ref:`R-vector <tb2j_plotter_R-vector>`), 
distances (see :ref:`max-distance <tb2j_plotter_max-distance>`,
:ref:`min-distance <tb2j_plotter_min-distance>` and
:ref:`distance <tb2j_plotter_distance>`), 
and template file (see :ref:`template <tb2j_plotter_template>`). 
The result is defined by logical conjugate of the specified conditions.

``--filename`` (or ``-f``) argument is required, the rest of them are optional.


Output files will have the following name structure: 
*output-name.display_data_type.png*

.. _tb2j_plotter_example:

Usage example
=============

Imagine you are executing the ``tb2j_plotter.py`` sccript from the 
folder *example* and your file structure looks like the following

.. code-block:: text

    example/
    ├── exchange.out
    └── output/
        
Now lets run the script

.. code-block:: console

    tb2j_plotter.py -f exchange.out 

After the execution your *example* folder will look similar to this
    
.. code-block:: text

    example/
    ├── exchange.out
    ├── exchange.iso.png
    ├── exchange.distance.png
    └── output/

Script produced two output files *exchange.iso.png*
and *exchange.distance.png*. Common seedname *exchange* comes by default 
(see :ref:`output-name <tb2j_plotter_output-name>`). *iso* and *distance* 
indicate the plotted data 
(see :ref:`what-to-plot <tb2j_plotter_what-to-plot>`). 

.. important::
    That output files are not located in *output* folder since the 
    current folder is used for output by default
    (see :ref:`output-dir <tb2j_plotter_output-dir>`). 
    
Lets save the output in the *output* folder:

.. code-block:: console

    tb2j_plotter.py -f exchange.out -op output

Now *example* folder should look like this

.. code-block:: text

    example/
    ├── exchange.out
    ├── exchange.iso.png
    ├── exchange.distance.png
    └── output/
        ├── exchange.iso.png
        └── exchange.distance.png

Output files have the same names, but they are saved in the *output* 
folder as your specifyed by ``-op`` argument.

It is not necessary to specify a path to the existing folder, 
for example try to execute

.. code-block:: console

    tb2j_plotter.py -f exchange.out -op output/bar/foo

The sript will create folder *bar* inside of the folder *output* and folder 
*foo* inside of the folder *bar*. The structure of the *example* folder now 
should look like that:

.. code-block:: text

    example/
    ├── exchange.out
    ├── exchange.iso.png
    ├── exchange.distance.png
    └── output/
        ├── exchange.iso.png
        |── exchange.distance.png
        └── bar/
            └── foo/
                ├── exchange.iso.png
                └── exchange.distance.png


Arguments
=========

.. _tb2j_plotter_filename:

``--filename``, ``-f``
----------------------
Relative or absulute path to the TB2J exchange output file, 
including the name and extention of the file.

    *required* : True

    *type* : str


.. _tb2j_plotter_mode:

``--mode``, ``-m``
------------------
Mode of plotting.

Two modes are supported: structure with the view from above 
and the plots with *value* over distance between bond and 
the center of the molecule.

    *default* : 2d

    *type* : str

    *choices* : all, 2d, molecule
    
.. hint::
    If you are plotting in molecule mode it is recommended to specify 
    ``--substrate_atoms`` argument.


.. _tb2j_plotter_substrate_atoms:

``--substrate_atoms``, ``-suba``
--------------------------------
Atoms from the substrate

Marks of atoms from the substracte (Same as in TB2J). 
You can specify only names. For example instead of "Cr12" one can provide 
"Cr" and then all Cr atoms will be considered as a substrate ones. 

    *default* : :py:class:`magnetic_atoms <.rad_tools.tb2j_tools.file_logic.ExchangeModel`

    *type* : str

    *nargs* : *


.. _tb2j_plotter_output-dir:

``--output-dir``, ``-op``
-------------------------
Relative or absolute path to the folder for saving outputs.

If the folder does not exist then it is created from the specified path.
The creation is applied recursevly to the path, starting from the right
until the existing folder is reached.

    *default* : current directory
        
    *type* : str


.. _tb2j_plotter_output-name:

``--output-name``, ``-on``
--------------------------
Seedname for the output files.

Output files will have the following name structure:
*output-name.display_data_type.png*

    *default* : exchange
        
    *type* : str

See also: :ref:`example <tb2j_plotter_example>`


.. _tb2j_plotter_what-to-plot:

``--what-to-plot``, ``-wtp``
----------------------------
Type of data for display.

Specifying the data for display at the graph. 
Everything is displayed by default, each value in a separate picture. 
Currently available for display: Isotropic exchange parameter, distance.

    *default* : all

    *type* : str

    *choices* : all, iso, distance


``--draw-cells``, ``-dc``
-------------------------
Whenever to draw the supercell`s shape.

If specified then the shape of all supercells 
presented in the model (after filtering) is drawn.

    *default* : False

    *action* : store_true


.. _tb2j_plotter_R-vector:

``--R-vector``, ``-R``
----------------------
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


.. _tb2j_plotter_max-distance:

``--max-distance``, ``-maxd``
-----------------------------
(<=) Maximum distance.

All the bonds with the distance beetwen atom 1 and atom 2 
greater then maximum distance are excluded from the model.

    *default* : None

    *type* : float


.. _tb2j_plotter_min-distance:

``--min-distance``, ``-mind``
-----------------------------
(>=) Minimum distance.

All the bonds with the distance beetwen atom 1 and atom 2 
lower then minimum distance are excluded from the model.

    *default* : None

    *type* : float


.. _tb2j_plotter_distance:

``--distance``, ``-d``
----------------------
(=) Exact distance.

Only the bonds with the exact distance remains in the model.

    *default* : None

    *type* : float

.. hint::
    There is no point in specifying maximum or minimum distance when 
    this parameter is specified.


.. _tb2j_plotter_template:

``--template``, ``-tf``
-----------------------
Relative or absolute path to the template file, 
including the name and extention of the file.

    *default* : None

    *type* : str

See also: :ref:`template <rad_make_template>`


.. _tb2j_plotter_double-bonds:

``--double-bonds``, ``-db``
---------------------------
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


.. _tb2j_plotter_scale_atoms:

``--scale-atoms``, ``-sa``
--------------------------
Scale for the size of atom marks.

Use it if you want to display atom marks bigger or smaller. 
Have to be positive.

    *default* : 1

    *type* : float


.. _tb2j_plotter_scale_data:

``--scale-data``, ``-sd``
-------------------------
Scale for the size of data text.

Use it if you want to display data text marks bigger or smaller. 
Have to be positive.

    *default* : 1

    *type* : float


.. _tb2j_plotter_title:

``--title``, ``t``
------------------
Title for the plots

Title will be displayed in the picture.

    *default* : None

    *type* : str
 

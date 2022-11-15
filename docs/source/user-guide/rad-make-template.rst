.. _rad-make-template:

************************
``rad-make-template.py``
************************

Script for the creation of template`s draft.

This script can provide blank template file or template file based on the TB2J
*exchange.out* file (see :ref:`--filename <rad-make-template_filename>`). 
Several filtering options are supported for the case of TB2J-based template file 
(:ref:`--R-vector <rad-make-template_R-vector>`, 
:ref:`--max-distance <rad-make-template_max-distance>`,
:ref:`--min-distance <rad-make-template_min-distance>`,
:ref:`--distance <rad-make-template_distance>`).

.. important::

    When template file is made on the base of TB2J file it is still necessary 
    to group bonds into neighbors and add names by hand afterwards.


Template specification
======================

Template file helps to choose particular bonds for filtering of exchange model 
or for extracting model with grouped parameters 
(i.e. :math:`J_1`, :math:`J_2`, ...).

Here is the draft of the template file which is provided by the script:

.. code-block:: text

    1   ====================
    2   Neighbors template:
    3   i j R_a R_b R_c
    4   --------------------
    5   J1 $J_1$
    6   atom1 atom2 0 0 0
    7   atom1 atom2 1 0 0
    8   atom1 atom1 -1 0 2
    9   --------------------
    10  J2
    11  atom2 atom1 9 5 -3
    12  atom1 atom2 1 4 0
    13  atom2 atom2 1 0 2
    14  ====================    

Line 1: Header of the file. 20 or more "=" symbols.

Line 2: Flag of the neighbors section, have to be in the file.

Line 3: Format of the bond specification line.

Line 4: Neighbor separator. Separates different neighbors
(:math:`J_1`, :math:`J_2`, ...) in the neighbors template file. 
20 or more "-" symbols.

Line 5: Name of the neighbor and the 
`LaTeX <https://www.latex-project.org/>`_ version of that name. 
Name have to be specified. Latex name is optional. 
Name and LaTeX name are separated by one or more spaces, 
as a consequence no spaces are allowed for both of them.

Line 6: First bond, which corresponds to the first neighbor (:math:`J_1`).
Format of the bond specification: 

.. code-block:: text

    atom1_mark atom2_mark R

Where R is a real-space vector of the unit cell in which the second atom is 
located (atom1 always located in R = (0, 0, 0)).

Line 7, 8: Specification of the second and third bond from the first neighbor.

Line 9: Neighbor separator.

Line 10: Name of the second neighbor.

.. note::
    There is no LaTeX name specifies for the second neighbor.

Lines 11-13: Specifications of the first, second and third bond, which are 
associated with the second neighbor.

Lines 14: Footer of the file. 20 or more "=" symbols.


Arguments
=========

.. _rad-make-template_output-dir:

-op, --output-dir
-----------------
Relative or absolute path to the folder for saving outputs.

    *default* : current directory
        
    *type* : str


.. _rad-make-template_output-name:

-on, --output-name
------------------
Template filename.

    *default* : template.txt

    *type* : str


.. _rad-make-template_filename:

-f, --filename
--------------
Relative or absulute path to the TB2J exchange output file, 
including the name and extention of the file.

    *default* : None 

    *type* : str


.. _rad-make-template_R-vector:

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


.. _rad-make-template_max-distance:

-maxd, --max-distance
---------------------
(<=) Maximum distance.

All the bonds with the distance beetwen atom 1 and atom 2 
greater than maximum distance are excluded from the model.

    *default* : None

    *type* : float


.. _rad-make-template_min-distance:

-mind, --min-distance
---------------------
(>=) Minimum distance.

All the bonds with the distance beetwen atom 1 and atom 2 
lower than minimum distance are excluded from the model.

    *default* : None

    *type* : float


.. _rad-make-template_distance:

-d, --distance
--------------
(=) Exact distance.

Only the bonds with the exact distance remains in the model.

    *default* : None

    *type* : float

.. hint::
    There is no point in specifying maximum or minimum distance when 
    this parameter is provided.

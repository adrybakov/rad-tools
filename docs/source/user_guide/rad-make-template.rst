.. _rad_make_template:

************************
``rad-make-template.py``
************************

Script for the creation of template`s template.


Template specification
======================

Template helps to choose particular bonds for filtering exchange model or 
for extracting model with grouped parameters (:math:`J_1`, :math:`J_2`,...).

Here is the template of the template which is provided by the scripts:

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
(:math:`J_1`, :math:`J_2`,...) in the neighbors template. 
20 or more "-" symbols.

Line 5: Name of the neighbor and the `LaTeX <https://www.latex-project.org/>`_ version of that name. Name have to be 
specified. Latex name is optional.

Line 6: First bond, which corresponds to the first neighbor (:math:`J_1`).
Format of the bond specification: 

.. code-block:: text

    Atom_1_mark Atom_2_mark R

Where R is a real-space vector of the unit cell in which the second atom is 
located (Atom 1 always located in R = (0, 0, 0)).

Line 7, 8: Specification of the second and third bond from the first neighbor.

Line 9: Neighbor separator.

Line 10: Name of the second neighbor.

.. note::
    There is no Latex name specifies for the second neighbor.

Lines 11-13: Specifications of the first, second and third bond, which is 
associated to the second neighbor.

Lines 14: Footer of the file. 20 or more "=" symbols.


Arguments
=========

.. _rad_make_template_output-dir:

``--output-dir``, ``-op``
-------------------------
Relative or absolute path to the folder for saving outputs.

    *default* : current directory
        
    *type* : str


.. _rad_make_template_output-name:

``--output-name``, ``-on``
--------------------------
Template file name, default "template.txt"

    *default* : template.txt

    *type* : str

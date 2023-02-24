.. _rad-make-template:

************************
``rad-make-template.py``
************************

Script for the creation of template`s draft.

This script can provide blank template file or template file based on the TB2J
*exchange.out* file (see :ref:`--tb2j-filename <rad-make-template_tb2j-filename>`). 
Several filtering options are supported for the case of TB2J-based template file 
(:ref:`--R-vector <rad-make-template_R-vector>`, 
:ref:`--max-distance <rad-make-template_max-distance>`,
:ref:`--min-distance <rad-make-template_min-distance>`,
:ref:`--distance <rad-make-template_distance>`).

.. important::

    When template file is made on the base of TB2J file it is still necessary 
    to group bonds into neighbors and add names by hand afterwards.


Usage example
=============
Example is based on the exchange.out file from 
:examples:`examples folder <rad-make-template>`. 

The simpliest case consist of template draft`s creation 
(see :ref:`template specification <template-draft>`). 
The following command will produce it:

.. code-block:: bash

    rad-make-template.py

For more advance user-case the file exchange.out from 
:examples:`examples folder <rad-make-template>` is used. 

Run the code:

.. code-block:: bash

    rad-make-template.py -tf exchange.out -on full_template.txt

It will produce the following file with the full template from the file:

.. dropdown:: full_template.txt

   .. literalinclude:: /../examples/rad-make-template/full_template.txt
    :language: text

This template is very long since original TB2J file includes a lot of 
interaction pairs, lets filter some of them and keep only the interactions 
with the distance <= 8 Angstrom.

.. code-block:: bash

    rad-make-template.py -tf exchange.out -on filtered_template.txt -maxd 8

.. dropdown:: filtered_template.txt

   .. literalinclude:: /../examples/rad-make-template/filtered_template.txt
    :language: text

Now for futher usage one only have to introduce grouping with respect to 
some exchange model to produce the final template file:

.. dropdown:: filtered_template_grouped.txt

   .. literalinclude:: /../examples/rad-make-template/filtered_template_grouped.txt
    :language: text

Check the :ref:`rad-make-template_arguments` section for more sorting options.

.. _template-draft:

Template specification
======================
Template file helps to choose particular bonds for filtering of exchange model 
or for extracting model with grouped parameters 
(i.e. :math:`J_1`, :math:`J_2`, ...).

Here is the draft of the template file which is provided by the script:

.. literalinclude:: /../examples/rad-make-template/template_demo.txt
    :language: text
    :linenos:   

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
    There is no LaTeX name specified for the second neighbor.

Lines 11-13: Specifications of the first, second and third bond, which are 
associated with the second neighbor.

Lines 14: Footer of the file. 20 or more "=" symbols.

.. _rad-make-template_arguments:

Arguments
=========

.. _rad-make-template_output-name:

-on, --output-name
------------------
Relative or absolute path to the template output file.

.. code-block:: text

    default : template

See also: :ref:`example <scripts_output-notes>`.

.. _rad-make-template_tb2j-filename:

-tf, --tb2j-filename
--------------------
Relative or absulute path to the TB2J exchange output file, 
including the name and extention of the file.

.. code-block:: text

    default : None 

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

.. code-block:: text

    default : None

.. _rad-make-template_max-distance:

-maxd, --max-distance
---------------------
(<=) Maximum distance.

All the bonds with the distance beetwen atom 1 and atom 2 
greater than maximum distance are excluded from the model.

.. code-block:: text

    default : None

.. _rad-make-template_min-distance:

-mind, --min-distance
---------------------
(>=) Minimum distance.

All the bonds with the distance beetwen atom 1 and atom 2 
lower than minimum distance are excluded from the model.

.. code-block:: text

    default : None

.. _rad-make-template_distance:

-d, --distance
--------------
(=) Exact distance.

Only the bonds with the exact distance remains in the model.

.. code-block:: text

    default : None

.. hint::
    There is no point in specifying maximum or minimum distance when 
    this parameter is provided.

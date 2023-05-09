.. _rad-make-template:

************************
``rad-make-template.py``
************************

Script for the creation of template`s draft.

This script can provide blank template file or template file based on the TB2J
"exchange.out" file (see :ref:`--input-filename <rad-make-template_input-filename>`). 
Several filtering options are supported for the case of TB2J-based template file 
(:ref:`--R-vector <rad-make-template_R-vector>`, 
:ref:`--max-distance <rad-make-template_max-distance>`,
:ref:`--min-distance <rad-make-template_min-distance>`,
:ref:`--distance <rad-make-template_distance>`).

.. important::

    When template file is made on the base of TB2J file it is still necessary 
    to group bonds into neighbors and add names manually.


Usage example
=============
Example is based on the exchange.out file from the
:examples:`examples folder <rad-make-template>`. 

Minimal usage scenario creates template draft`s with the command
(see :ref:`template specification <template-draft>`):

.. code-block:: bash

    rad-make-template.py

For more advance user-case the file exchange.out from 
:examples:`examples folder <rad-make-template>` is used. 

The code:

.. code-block:: bash

    rad-make-template.py -if exchange.out -on full_template.txt

produces the file with the full template:

.. dropdown:: full_template.txt

   .. literalinclude:: /../examples/rad-make-template/full_template.txt
    :language: text

This template is very long since original TB2J file includes a lot of 
interaction pairs, lets filter some of them and keep only the interactions 
with the distance <= 8 Angstroms.

.. code-block:: bash

    rad-make-template.py -if exchange.out -on filtered_template.txt -maxd 8

.. dropdown:: filtered_template.txt

   .. literalinclude:: /../examples/rad-make-template/filtered_template.txt
    :language: text

For further usage of the template it is necessary to introduce 
grouping with respect to some exchange model:

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

Line 1: Date and time of file creation.

Line 2: Blank line.

Line 3: Header of the file (20 or more "=" symbols).

Line 4: Flag of the neighbors section, have to be in the file.

Line 5: Format of the bond specification line.

Line 6: Neighbour separator (20 or more "-" symbols). 
Separates different neighbors
:math:`(J_1, J_2, \dots)` in the template file.

Line 7: Name of the neighbour and the |latex|_ version of that name. 
Name have to be specified. LaTeX name is optional. 
Name and LaTeX name are separated by one or more spaces, 
as a consequence no spaces are allowed for both of them.

Line 8: First bond, which corresponds to the first neighbour (:math:`J_1`).
Format of the bond specification: 

.. code-block:: text

    atom1_mark atom2_mark R

Where R is a real-space vector of the unit cell in which the second atom is 
located (atom1 always located in R = (0, 0, 0)).

Line 9, 10: Specification of the second and third bond from the first neighbour.

Line 11: Neighbour separator.

Line 12: Name of the second neighbour.

.. note::
    There is no LaTeX name specified for the second neighbour.

Lines 14-15: Specifications of the first, second and third bond, which are 
associated with the second neighbour.

Lines 16: Footer of the file (20 or more "=" symbols).

.. _rad-make-template_arguments:

Arguments
=========

.. _rad-make-template_output-name:

-on, --output-name
------------------
Name for the template output file.

.. code-block:: text

    default : template

See also: :ref:`example <output-notes>`.

.. _rad-make-template_input-filename:

-if, --input-filename
---------------------
Relative or absolute path to the 'exchange.out' file, 
including name and extension of the file.

.. code-block:: text

    default : None 

.. versionchanged:: 0.5.12 Renamed from "tb2j_filename"

.. _rad-make-template_R-vector:

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

.. _rad-make-template_max-distance:

-maxd, --max-distance
---------------------
(<=) Maximum distance.

All the bonds with the distance between atom 1 and atom 2 
greater than maximum distance are excluded from the model.

.. code-block:: text

    default : None

.. _rad-make-template_min-distance:

-mind, --min-distance
---------------------
(>=) Minimum distance.

All the bonds with the distance between atom 1 and atom 2 
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


.. _rad-make-template_verbose:

-v, --verbose
-------------
Verbose output, propagates to the called methods.

.. code-block:: text

    default : False

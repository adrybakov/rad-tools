.. _guide_io_radtools:

.. currentmodule:: radtools

*******************
RAD-tools interface
*******************

Yes, it is interface to itself.

Function :py:func:`read_template` reads exchange Hamiltonian template file 
and constructs :py:class:`ExchangeTemplate` from it.

Here is full template file specification:


.. _template-draft:

Template specification
======================
Template file helps to choose particular bonds for filtering of exchange Hamiltonian 
or for extracting model with grouped parameters 
(i.e. :math:`J_1`, :math:`J_2`, ...).

Here is the draft of the template file which is provided by the script:

.. literalinclude:: ../../../../examples/rad-make-template/template_demo.txt
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
.. _rad-plot-tb2j-magnons:

************************
rad-plot-tb2j-magnons.py
************************

Script for plotting magnon dispersion from |TB2J|_ results.

.. versionadded:: 0.7.12 

Scrip plots magnon dispersion spectra following the algorithm described in 
:ref:`library_magnon-dispersion-method`.

It requires |TB2J|_ output file "exchange.out" and 
the information about the magnetic ground state.

Example is based on the files from 
:examples:`examples folder <rad-plot-tb2j-magnons>`.

Ground state input
==================

Two types of ground state is supported: 

* arbitrary directions of spin in the unit cell (or supercell)
* single-Q incommensurate structure.

To define direction of spins in the supercell use :ref:`rad-plot-tb2j-magnons_spin` parameter:

.. code-block:: bash

    rad-plot-tb2j-magnons.py -if exchange.out -s Cr1 0 0 1 Cr2 0 0 1

In the input file "exchange.out" six atoms are present: Br1, Cr1, S1, Br2, Cr2, S2.
Only Cr1 and Cr2 have exchange interaction between them. Therefore, it is neccesary to specify
spin vectors for Cr1 and Cr2. You can specify spin vectors for all atoms, but it is not
necessary.

Spin spiral is defined by two vectors: 

* :math:`\vec{Q}` (:ref:`rad-plot-tb2j-magnons_spiral-vector`)

It is relative to the model reciprocal cell.

* Global rotation axis :math:`\vec{n}` (:ref:`rad-plot-tb2j-magnons_rotation-axis`)

It is given in absolute coordinate in a real space. Only the direction of the vector matters.

Template file
=============

Exchange template file (see :ref:`template-draft`) can be used to forme the model or
to filter the spin Hamiltonian.

Filtering of the model
======================

For filtering the spin Hamiltonian there are a few options available:

* :ref:`--max_distance <rad-plot-tb2j_max-distance>`
* :ref:`--min_distance <rad-plot-tb2j_min-distance>`
* :ref:`--R-vector <rad-plot-tb2j_R-vector>`
* :ref:`--template <rad-plot-tb2j_template-file>`

.. _rad-plot-tb2j-magnons_arguments:

Arguments
=========

.. _rad-plot-tb2j-magnons_input-filename:

-if, --input-filename
---------------------
Relative or absolute path to the "exchange.out" file,
including the name and extension of the file itself.

.. code-block:: text

    required
    type : str

.. _rad-plot-tb2j-magnons_template-file:

-tf, --template-file
--------------------
Relative or absolute path to the template file, 
including the name and extension of the file.

.. code-block:: text

    required
    type : str

.. _rad-plot-tb2j-magnons_output-name:

-on, --output-name
------------------
Seedname for the output files.

If this parameter is not specified, the result are printed in 
standard output stream. 

.. code-block:: text

    default : None
    type : str

See also: :ref:`example <output-notes>`.

.. _rad-plot-tb2j-magnons_output-path:

-op, --output-path
------------------
Relative or absolute path to the folder for saving outputs.

If the folder does not exist then it is created from the specified path.
The creation is applied recursively to the path, starting from the right
until the existing folder is reached.

.. code-block:: text

    default : current directory

See also: :ref:`example <output-notes>`.

.. _rad-plot-tb2j-magnons_spin:

-s, --spin
----------
Spin of the atoms in the model.

For each atom, which has at least one bond connected to it is necessary to specify
spin vector. The spin vector is specified in the form of atom`s name followed by
three numbers, separated by spaces. 
The numbers represent the x, y, and z components of the spin vector.

.. code-block:: text

    default : None

.. _rad-plot-tb2j-magnons_spiral-vector:

-Q, --spiral-vector
-------------------
Spin spiral vector. Relative to the reciprocal cell.

.. code-block:: text

    default : None
    type : float

.. _rad-plot-tb2j-magnons_rotation-axis:

-ra, --rotation-axis
--------------------
Direction of global rotation axis. In absolute coordinates in real space.

.. code-block:: text

    default : None
    type : float
    nargs : 3
    
.. _rad-plot-tb2j-magnons_path:

-p, --path
----------
Path in reciprocal space for the magnon dispersion.

.. code-block:: text

    default : None
    type : str

.. _rad-plot-tb2j-magnons_form-model:

-fm, --form-model
---------------------
Whether to form the model based on the template.

.. code-block:: text

    default : False
    type : bool

.. _rad-plot-tb2j-magnons_R-vector:

-R, --R-vector
--------------
R vectors for filtering the spin Hamiltonian.

In TB2J outputs the bond is defined by atom 1 (from) and atom 2 (to). 
Atom 1 is always located in (0, 0, 0) unit cell, while atom 2 is located in 
R = (i, j, k) unit cell. This parameter tells the script to keep only the 
bonds for which atom 2 is located in one of specified R supercells. 
Supercells are specified by a set of integers separated by spaces. 
They are grouped by three starting from the left and forms a set 
of R vectors. If the last group contains 1 or 2 integers they are ignored.

.. code-block:: text

    default : None

.. _rad-plot-tb2j-magnons_max-distance:

-maxd, --max-distance
---------------------
(<=) Maximum distance.

All the bonds with the distance between atom 1 and atom 2 
greater than maximum distance are excluded from the model.

.. code-block:: text

    default : None

.. _rad-plot-tb2j-magnons_min-distance:

-mind, --min-distance
---------------------
(>=) Minimum distance.

All the bonds with the distance between atom 1 and atom 2 
lower than minimum distance are excluded from the model.

.. code-block:: text

    default : None

.. _rad-plot-tb2j-magnons_save-txt:

-st, --save-txt
---------------
Whether to save data to .txt file. Two files appears: 
"output-name.txt" and "output-name_info.txt". First one contains raw data of the graph,
second one contains information about the parameters.

.. code-block:: text

    default : False

.. _rad-plot-tb2j-magnons_interactive:

-i, --interactive
-----------------
Whether to show interactive plot.

.. code-block:: text

    default : False

.. _rad-plot-tb2j-magnons_verbose:

-v, --verbose
--------------
Verbose output, propagates to the called methods.

.. code-block:: text

    default : False

.. _rad-plot-tb2j-magnons_bravais-type:

-bt, --bravais-type
--------------------

Bravais lattice type. 
If not provided, the type is identified automatically.

It does not force the Bravais lattice type on the model,
but tries to reach the desired type by reducing the 
numerical accuracy in the :py:func:`lepage` algorithm.

.. code-block:: text

    default : None
    type : str
    choices : CUB, FCC, BCC, TET, BCT, ORC, ORCF, ORCI, ORCC, HEX, RHL, MCL, MCLC, TRI

.. _rad-plot-tb2j-magnons_join-output:

-jo, --join-output
------------------
Whether to join the output files into a single file.

.. code-block:: text

    default : False

.. _rad-plot-tb2j-magnons_nodmi:

-nodmi
------
Whether to ignore DMI in the spinham.

.. code-block:: text

    default : False

.. _rad-plot-tb2j-magnons_no-anisotropic:

-noa, --no-anisotropic
----------------------
Whether to ignore anisotropic symmetric exchange in the spinham.

.. code-block:: text

    default : False

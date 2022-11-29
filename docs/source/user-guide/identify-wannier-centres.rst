.. _identify-wannier-centres:

*******************************
``identify-wannier-centres.py``
*******************************

Identify wannier centres with respect to the atom or to 
the point between the atom`s pair.

Use :ref:`--span <identify-wannier-centres_span>` option to increase the 
distance span for the search.

* If the centre is associated with some atom
    then ``-> atom`` will appear in the output file.

* If the centre is associated with some centre point between two atoms
    then ``-> atom1-atom2`` will appear in the output file.

* If the centre is unidentified 
    then ``-> None`` 
    will appear in the output file and additional information 
    will be present in the console.

* If the centre is equally close to the atom and to the point between two atoms 
    then ``-> atom or atom1-atom2`` will appear in the output file 
    (tolerance :math:`10^{-5}`).

Usage example
=============

The example_centres.xyz file looks like this:

.. literalinclude:: /../examples/identify-wannier-centres/example_centres.xyz

Lets run the code:

.. code-block:: console

    identify-wannier-centres.py example_centres.xyz

This command creates an output file example_centre.xyz_identified 
in the directory of the input file with the following content:

.. literalinclude:: /../examples/identify-wannier-centres/example_centres.xyz_identified

and produce the following output in the console:

.. code-block::
    :emphasize-lines: 1,6

    Centre [0.90034182 1.20456455 0.02140611] unindentified, 
    try to increase --span
    span limit = 0.1
    minimum distance to the atom = 0.10018733 (Br1)
    minimum distance to the bond`s centre = 1.20808898 (Br1-Cr1)

    Centre [2.70102542 3.61369366 5.54938336] unindentified, 
    try to increase --span
    span limit = 0.1
    minimum distance to the atom = 0.10017450 (Br2)
    minimum distance to the bond`s centre = 1.20812019 (Br2-Cr2)
    
which means that two centres are not identified. 
The script provides the distance to the closest atom 
and to the closest centre of the bond between some pair of the atoms
for each unidentified centre. 

As one can see first unidentified centre is quite close to the Br1 and 
the second one to the Br2, let us extend the span a little bit in order to 
correctly identify all the atoms:

.. code-block:: console

    identify-wannier-centres.py example_centres.xyz --span 0.11 --output-name example_centres.xyz_bigger_span

This command produces example_centres.xyz_bigger_span file:

.. literalinclude:: /../examples/identify-wannier-centres/example_centres.xyz_bigger_span

where all the centres are correctly Identified.

Arguments
=========

.. _identify-wannier-centres_filename:

filename
--------
Rellative or absolute path to the _centres.xyz file

Identified Wannier centres will be store in the "filename_identified" file.

    *type : str*

.. _identify-wannier-centres_span:

-s, --span
----------
Distance tolerance between centre and atom. (in Angstrom)

If some atoms remains unidentified try to increase the span

    *default* : 0.1

    *type* : float

.. _identify-wannier-centres_output-dir:

-op, --output-dir
-----------------
Relative or absolute path to the folder for saving outputs.

    *default* : the directory of the input file
        
    *type* : str

.. _identify-wannier-centres_output-name:

-on, --output-name
------------------
Seedname for the output files.

    *default* : Name of the input file + "_identified"

    *type* : str

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

The example_centres.xyz file looks like this: ::

        28
     Wannier centres, written by Wannier90 on16Sep2022 at 12:23:35 
    X          2.70102543       1.20456455       1.74829100
    X          2.70102543       1.20456455       1.73136237
    X          2.70102542       1.20456454       1.71844477
    X          2.70102540       1.20456454       1.74850456
    X          2.70102542       1.20456455       1.74171453
    X          0.90034180       3.61369365       3.82242305
    X          0.90034181       3.61369365       3.83936190
    X          0.90034181       3.61369365       3.85226490
    X          0.90034182       3.61369365       3.82220340
    X          0.90034181       3.61369365       3.82899496
    X          2.70102543       3.61369365       2.21721708
    X          2.70102542       3.61369367       2.29147179
    X          2.70102540       3.61369365       2.17808141
    X          0.90034181       1.20456453       3.35355250
    X          0.90034181       1.20456457       3.27930790
    X          0.90034180       1.20456458       3.39269417
    X          0.90034181       1.20456455      -0.05695451
    X          0.90034182       1.20456455       0.02140611
    X          0.90034181       1.20456456      -0.04741017
    X          2.70102542       3.61369365       5.62773652
    X          2.70102542       3.61369366       5.54938336
    X          2.70102542       3.61369365       5.61819448
    Br         0.90034179       1.20456453      -0.07878122
    Br         2.70102536       3.61369358       5.64955786
    Cr         2.70102536       1.20456453       1.73263498
    Cr         0.90034179       3.61369358       3.83807392
    S          2.70102536       3.61369358       2.20990738
    S          0.90034179       1.20456453       3.36085407

Lets run the code:

.. code-block:: console

    identify-wannier-centres.py example_centres.xyz

This command will create an output file example_centre.xyz_identified 
in the directory of the input file with the following content: ::

    28
     Wannier centres, written by Wannier90 on16Sep2022 at 12:23:35 
    X          2.70102543       1.20456455       1.74829100   ->   Cr1 
    X          2.70102543       1.20456455       1.73136237   ->   Cr1 
    X          2.70102542       1.20456454       1.71844477   ->   Cr1 
    X          2.70102540       1.20456454       1.74850456   ->   Cr1 
    X          2.70102542       1.20456455       1.74171453   ->   Cr1 
    X          0.90034180       3.61369365       3.82242305   ->   Cr2 
    X          0.90034181       3.61369365       3.83936190   ->   Cr2 
    X          0.90034181       3.61369365       3.85226490   ->   Cr2 
    X          0.90034182       3.61369365       3.82220340   ->   Cr2 
    X          0.90034181       3.61369365       3.82899496   ->   Cr2 
    X          2.70102543       3.61369365       2.21721708   ->   S1  
    X          2.70102542       3.61369367       2.29147179   ->   S1  
    X          2.70102540       3.61369365       2.17808141   ->   S1  
    X          0.90034181       1.20456453       3.35355250   ->   S2  
    X          0.90034181       1.20456457       3.27930790   ->   S2  
    X          0.90034180       1.20456458       3.39269417   ->   S2  
    X          0.90034181       1.20456455      -0.05695451   ->   Br1 
    X          0.90034182       1.20456455       0.02140611   ->   None
    X          0.90034181       1.20456456      -0.04741017   ->   Br1 
    X          2.70102542       3.61369365       5.62773652   ->   Br2 
    X          2.70102542       3.61369366       5.54938336   ->   None
    X          2.70102542       3.61369365       5.61819448   ->   Br2 
    Br         0.90034179       1.20456453      -0.07878122   ->   Br1 
    Br         2.70102536       3.61369358       5.64955786   ->   Br2 
    Cr         2.70102536       1.20456453       1.73263498   ->   Cr1 
    Cr         0.90034179       3.61369358       3.83807392   ->   Cr2 
    S          2.70102536       3.61369358       2.20990738   ->   S1  
    S          0.90034179       1.20456453       3.36085407   ->   S2 

and produce the following output in the console:

.. code-block::
    :emphasize-lines: 1,6

    Centre [0.90034182 1.20456455 0.02140611] unindentified, try to increase --span
    span limit = 0.1
    centre`s min span = 0.10018733 (with Br1 atom)
    centre`s min span = 1.50401109 (with centre point between Br1-Br1 atoms)

    Centre [2.70102542 3.61369366 5.54938336] unindentified, try to increase --span
    span limit = 0.1
    centre`s min span = 0.10017450 (with Br2 atom)
    centre`s min span = 3.37808251 (with centre point between Br2-Cr1 atoms)

which means that two centres was not identified. The script provides the span 
the closest atom and to the closest centre point between some pair of the atoms
for each unidentified centre.

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

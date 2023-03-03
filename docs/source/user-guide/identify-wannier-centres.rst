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
Example is based on the files from 
:examples:`examples folder <identify-wannier-centres>`. 

The example_centres.xyz file looks like this:

.. literalinclude:: /../examples/identify-wannier-centres/example_centres.xyz

Lets run the code:

.. code-block:: bash

    identify-wannier-centres.py example_centres.xyz

This command creates an output file example_centre.xyz_identified 
in the directory of the input file with the following content:

.. literalinclude:: /../examples/identify-wannier-centres/example_centres.xyz_identified

and produce the following output in the console:

.. literalinclude:: /../examples/identify-wannier-centres/console_output.txt
   :language: text
    
which means that two centres are not identified. 
The script provides in the console the distance to the closest atom 
and to the closest centre of the bond between some pair of the atoms
for each unidentified centre. In the output file some information about 
inedentified centres are provided as well.

As one can see first unidentified centre is quite close to the Br1 and 
the second one to the Br2, let us extend the span a little bit in order to 
correctly identify all the atoms:

.. code-block:: bash

    identify-wannier-centres.py example_centres.xyz --span 0.11 --output-name example_centres.xyz_bigger_span

This command produces example_centres.xyz_bigger_span file:

.. literalinclude:: /../examples/identify-wannier-centres/example_centres.xyz_bigger_span

where all the centres are correctly Identified.

Arguments
=========

.. _identify-wannier-centres_input-filename:

input-filename
--------------
Rellative or absolute path to the _centres.xyz file

Identified Wannier centres will be store in the "input-filename"_identified file.

.. code-block:: text

    required

.. _identify-wannier-centres_span:

-s, --span
----------
Distance tolerance between centre and atom. (in Angstrom)

If some atoms remains unidentified try to increase the span

.. code-block:: text

    default : 0.1

.. _identify-wannier-centres_output-dir:

-op, --output-dir
-----------------
Relative or absolute path to the folder for saving outputs.

.. code-block:: text

    default : the directory of the input file

See also: :ref:`example <scripts_output-notes>`.

.. _identify-wannier-centres_output-name:

-on, --output-name
------------------
Seedname for the output files.

.. code-block:: text

    default : Name of the input file + "_identified"

See also: :ref:`example <scripts_output-notes>`.

.. _identify-wannier-centres_no-colour:

-nc, --no-colour
----------------
Turn off coloured output.

.. code-block:: text

    default : False

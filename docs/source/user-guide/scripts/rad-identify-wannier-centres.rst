.. _rad-identify-wannier-centres:

****************************
rad-identify-wannier-centres
****************************

Identifies wannier centres with respect to the atom or to
the point between the atom`s pair.

.. versionchanged:: 0.6 Renamed from ``identify-wannier-centres.py``

.. versionchanged:: 0.9.0 Renamed from ``rad-identify-wannier-centres.py`` to ``rad-identify-wannier-centres``

Use :ref:`--span <rad-identify-wannier-centres_span>` option to increase the
distance span for the search.

* If the centre is associated with some atom
    then ``-> atom`` appears in the output file.

* If the centre is unidentified
    then ``-> None``
    appears in the output file and additional information
    is present in the console and in the output file.

Usage example
=============
Example is based on the files from
:examples:`examples folder <rad-identify-wannier-centres>`.

The "example_centres.xyz" file looks like this:

.. literalinclude:: /../examples/rad-identify-wannier-centres/example_centres.xyz

Lets run the code:

.. code-block:: bash

    rad-identify-wannier-centres -if example_centres.xyz

This command creates an output file "example_centre.xyz_identified"
in the directory of the input file with the following content:

.. literalinclude:: /../examples/rad-identify-wannier-centres/example_centres.xyz_identified

which means that two centres are not identified.
The script provides the distance to the closest atom
for each unidentified centre. In the output file some information about
unidentified centres is provided as well.

As one can see first unidentified centre is quite close to the
:math:`\text{Br}_1` and the second one to the :math:`\text{Br}_2`,
let us extend the span a little bit in order to correctly identify
all the atoms:

.. code-block:: bash

    rad-identify-wannier-centres -if example_centres.xyz --span 0.11 --output-name example_centres.xyz_bigger_span

This command produces example_centres.xyz_bigger_span file:

.. literalinclude:: /../examples/rad-identify-wannier-centres/example_centres.xyz_bigger_span

where all centres are correctly identified.

.. hint::
    Set span to 0 to see the distances to the closest atom in the output file for each
    wannier center.

.. _rad-identify-wannier-centres_arguments:

Arguments
=========

.. _rad-identify-wannier-centres_input-filename:

-if, --input-filename
---------------------
Relative or absolute path to the "_centres.xyz" file

Identified Wannier centres are stored in the "filename_identified" file.

.. code-block:: text

    required
    type: str

.. versionchanged:: 0.8.0 Renamed from ``filename``

.. _rad-identify-wannier-centres_span:

-s, --span
----------
Distance tolerance between centre and atom. (in Angstroms)

If some centres remain unidentified try to increase the span.

.. code-block:: text

    default: 0.1
    type: float


.. _rad-identify-wannier-centres_output-name:

-on, --output-name
------------------
Seedname for the output files.

See also: :ref:`example <output-notes>`

.. code-block:: text

    default: ""
    type: str

.. _rad-extract-tb2j:

****************
rad-extract-tb2j
****************

Script for extracting of template-based Hamiltonian from
|TB2J|_ results.

.. versionchanged: 0.6 Renamed from ``tb2j-extractor.py``

.. versionchanged: 0.9.0 Renamed from ``rad-extract-tb2j.py`` to ``rad-extract-tb2j``

If :ref:`--output_name <rad-extract-tb2j_output-name>` is not provided, then the result is
printed to the console, otherwise it is written to the text file.

Exchange types
==============

There are several options to control which exchange types are included in
the summary:

* :ref:`--no-anisotropic <rad-extract-tb2j_no-anisotropic>`
* :ref:`--nodmi <rad-extract-tb2j_nodmi>`
* :ref:`--no-matrix <rad-extract-tb2j_no-matrix>`

Isotropic exchange is always included in the output. By default everything is
included.

Usage example
=============

Example is based on the files from
:examples:`examples folder <rad-extract-tb2j>`.

There are two modes in which exchange summary can be printed:

* With the model of the template
* Full Hamiltonian

In the first case the model of the template template is constructed from the Hamiltonian
and exchange output is grouped based on the names provided in the template:

.. code-block:: bash

    rad-extract-tb2j -fm -if exchange.out -tf template.txt -on summary_formed_model.txt

.. dropdown:: summary with forced symmetry

   .. literalinclude:: /../examples/rad-extract-tb2j/summary_formed_model.txt
    :language: text

In the second case exchange summary is printed for every bond in the Hamiltonian.
Hamiltonian is filtered with the :ref:`--max-distance <rad-extract-tb2j_max-distance>`,
:ref:`--min-distance <rad-extract-tb2j_min-distance>` and
:ref:`--template-file <rad-extract-tb2j_template-file>`.

.. code-block:: bash

    rad-extract-tb2j -if exchange.out -tf template.txt -on summary.txt

.. dropdown:: summary without forced symmetry

   .. literalinclude:: /../examples/rad-extract-tb2j/summary.txt
    :language: text

.. _rad-extract-tb2j_arguments:

Arguments
=========

.. _rad-extract-tb2j_input-filename:

-if, --input-filename
---------------------
Relative or absolute path to the "exchange.out" file,
including the name and extension of the file itself.

.. code-block:: text

    required
    type: str


.. _rad-extract-tb2j_template-file:

-tf, --template-file
--------------------
Relative or absolute path to the template file,
including the name and extension of the file.

.. code-block:: text

    optional
    type: str


.. _rad-extract-tb2j_output-name:

-on, --output-name
------------------
Name of the output files.

If this parameter is not specified, the result is printed in the
standard output stream.

.. code-block:: text

    optional
    type: str


.. _rad-extract-tb2j_decimals:

-d, --decimals
--------------
Decimals after the comma for the exchange values.

.. code-block:: text

    default: 4
    type: int

.. versionchanged:: 0.5.17 Renamed from ``-acc``/``--accuracy``.

.. _rad-extract-tb2j_form-model:

-fm, --form-model
-----------------
Whether to form the model from the template.

.. code-block:: text

    default: False
    type: bool

.. versionchanged:: 0.8.0 Renamed from ``-fs``/``--force-symmetry``.

.. _rad-extract-tb2j_no-anisotropic:

-noa, --no-anisotropic
----------------------
Whether to output anisotropic exchange.

.. code-block:: text

    default: False
    type: bool

.. versionchanged:: 0.8.0 Renamed from ``-a``/``--anisotropic``.

.. _rad-extract-tb2j_no-matrix:

-nom, --no-matrix
-----------------
Whether to output whole matrix exchange.

.. code-block:: text

    default: False
    type: bool

.. versionchanged:: 0.8.0 Renamed from ``-m``/``--matrix``.

.. _rad-extract-tb2j_nodmi:

-nodmi
------
Whether to output DMI exchange.

.. code-block:: text

    default: False
    type: bool

.. versionchanged:: 0.8.0 Renamed from ``-dmi``.

.. _rad-extract-tb2j_verbose:

-v, --verbose
-------------
Verbose output, propagates to the called methods.

.. code-block:: text

    default: False
    type: bool


.. _rad-extract-tb2j_max-distance:

-maxd, --max-distance
---------------------
(<=) Maximum distance.

All the bonds with the distance between atom 1 and atom 2
greater than maximum distance are excluded from the model.

.. code-block:: text

    optional
    type: float

.. versionadded:: 0.8.0

.. _rad-extract-tb2j_min-distance:

-mind, --min-distance
---------------------
(>=) Minimum distance.

All the bonds with the distance between atom 1 and atom 2
lower than minimum distance are excluded from the Hamiltonian.

.. code-block:: text

    optional
    type: float

.. versionadded:: 0.8.0

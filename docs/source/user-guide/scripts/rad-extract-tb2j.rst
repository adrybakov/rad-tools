.. _rad-extract-tb2j:

*******************
rad-extract-tb2j.py
*******************

Script for extracting of template-based model from 
|TB2J|_ results.

.. versionchanged: 0.6 Renamed from ``tb2j-extractor.py``

If :ref:`--output_name <rad-extract-tb2j_output-name>` is not provided the result is 
passed to the console, otherwise it is written to the file with first 3 lines 
giving information about the source of exchange, date and time.

Exchange types
==============

There are several options to control which exchange types are included in 
the summary:

* :ref:`--isotropic <rad-extract-tb2j_isotropic>`
* :ref:`--anisotropic <rad-extract-tb2j_anisotropic>`
* :ref:`--dmi <rad-extract-tb2j_dmi>`
* :ref:`--matrix <rad-extract-tb2j_matrix>`
* :ref:`--all <rad-extract-tb2j_all>`

Format of the output block for each exchange type is provided in the 
:ref:`arguments <rad-extract-tb2j_arguments>` section.

Usage example
=============

Example is based on the files from 
:examples:`examples folder <rad-extract-tb2j>`. 

There are two modes in which exchange summary can be printed: 

* with the symmetry of the template and 
* with the model filtered based on the template.

In the first case the symmetry of the template is forced on the model and 
exchange output is grouped based on the names provided in the template:

.. code-block:: bash

    rad-extract-tb2j.py -all -fs -if exchange.out -tf template.txt -on summary_forced_symmetry

.. dropdown:: summary with forced symmetry

   .. literalinclude:: /../examples/rad-extract-tb2j/summary_forced_symmetry.txt
    :language: text

In the second case exchange summary is printed for every bond in the 
template file and no additional symmetry constrains are assumed on the model:

.. code-block:: bash

    rad-extract-tb2j.py -all -if exchange.out -tf template.txt -on summary

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
    type : str

.. _rad-extract-tb2j_template-file:

-tf, --template-file
--------------------
Relative or absolute path to the template file, 
including the name and extension of the file.

.. code-block:: text

    required
    type : str


See also: :ref:`template <rad-make-template>`


.. _rad-extract-tb2j_output-path:

-op, --output-path
------------------
Relative or absolute path to the folder for saving outputs.

If the folder does not exist then it is created from the specified path.
The creation is applied recursively to the path, starting from the right
until the existing folder is reached.

.. code-block:: text

    default : current directory (".")
    type : str

See also: :ref:`example <output-notes>`.


.. _rad-extract-tb2j_output-name:

-on, --output-name
------------------
Seedname for the output files.

If this parameter is not specified, the result are printed in 
standard output stream. 

.. code-block:: text

    default : None
    type : str

See also: :ref:`example <output-notes>`.


.. _rad-extract-tb2j_decimals:

-d, --decimals
--------------
Decimals after the comma for the exchange values.

.. code-block:: text

    default : 4
    type : int

.. versionchanged:: 0.5.17 Renamed from "-acc"/"--accuracy".

.. _rad-extract-tb2j_force-symmetry:

-fs, --force-symmetry
---------------------
Whether to force the symmetry of the template on the model.

.. code-block:: text

    default : False
    type : bool


.. _rad-extract-tb2j_isotropic:

-i, --isotropic
---------------
Whether to output isotropic exchange.

.. code-block:: text

    default : False
    type : bool

Section format:

.. code-block:: text

        Isotropic: J


.. _rad-extract-tb2j_anisotropic:

-a, --anisotropic
-----------------
Whether to output anisotropic exchange.

.. code-block:: text

    default : False
    type : bool

Section format:

.. code-block:: text

        Anisotropic: 
            Jxx Jxy Jxz
            Jxy Jyy Jyz
            Jxz Jyz Jzz


.. _rad-extract-tb2j_matrix:

-m, --matrix
------------
Whether to output the whole matrix of exchange.

.. code-block:: text

    default : False
    type : bool

Section format:

.. code-block:: text

        Matrix: 
            Jxx Jxy Jxz
            Jyx Jyy Jyz
            Jzx Jzy Jzz


.. _rad-extract-tb2j_dmi:

-dmi
----
Whether to output DMI exchange.

.. code-block:: text

    default : False
    type : bool

Section format in the case of forced symmetry:

.. code-block:: text

        |DMI|: |DMI|
        |DMI/J|: |DMI/J|
        DMI: DMI_x DMI_y DMI_z (Atom1 Atom2 Ra Rb Rc)
        ...

Otherwise:

.. code-block:: text

        |DMI|: |DMI|
        |DMI/J|: |DMI/J|
        DMI: DMI_x DMI_y DMI_z


.. _rad-extract-tb2j_all:

-all
----
Whether to output all types of exchange.

.. code-block:: text

    default : False
    type : bool


.. _rad-extract-tb2j_verbose:

-v, --verbose
-------------
Verbose output, propagates to the called methods.

.. code-block:: text

    default : False

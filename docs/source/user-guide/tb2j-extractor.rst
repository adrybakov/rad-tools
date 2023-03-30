.. _tb2j-extractor:

*********************
``tb2j-extractor.py``
*********************

Script for extracting of template-based model from 
`TB2J <https://tb2j.readthedocs.io/en/latest/>`_ results.

If :ref:`--output_name <tb2j-extractor_output-name>` is not provided the result is 
passed to the console, otherwise it is written to the file with first 3 lines 
giving information about the source of exchange and datetime.

Exchange types
==============

There are several options to control which exchange types are included in 
the summary:

    * :ref:`--isotropic <tb2j-extractor_isotropic>`
    * :ref:`--anisotropic <tb2j-extractor_anisotropic>`
    * :ref:`--dmi <tb2j-extractor_dmi>`
    * :ref:`--matrix <tb2j-extractor_matrix>`
    * :ref:`--all <tb2j-extractor_all>`

They work for both modes with the correspondent blocks having the following form
Format of the output block for each exchange type is provided in the 
arguments section.

Usage example
=============

Example is based on the files from 
:examples:`examples folder <tb2j-extractor>`. 

There are two modes in which exchange summary can be printed: 
with the symmetry of the template and 
with the model filtered based on the template.

In the first case the symmetry of the template is forced on the model and 
exchange output is grouped based on the names provided in the template:

.. code-block:: bash

    tb2j-extractor.py -all -fs -if exchange.out -tf template.txt -on summary_forced_symmetry

.. dropdown:: summary with forced symmetry

   .. literalinclude:: /../examples/tb2j-extractor/summary_forced_symmetry.txt
    :language: text

In the second case exchange summary is printed for every bond in the 
template file and no additional symmetry constrains are assumed on the model:

.. code-block:: bash

    tb2j-extractor.py -all -if exchange.out -tf template.txt -on summary

.. dropdown:: summary without forced symmetry

   .. literalinclude:: /../examples/tb2j-extractor/summary.txt
    :language: text

Arguments
=========

.. _tb2j-extractor_input-filename:

-if, --input-filename
---------------------
Relative or absolute path to the *exchange.out* file,
including the name and extension of the file itself.

.. code-block:: text

    required
    type : str

.. _tb2j-extractor_template-file:

-tf, --template-file
--------------------
Relative or absolute path to the template file, 
including the name and extension of the file.

.. code-block:: text

    required
    type : str


See also: :ref:`template <rad-make-template>`


.. _tb2j-extractor_output-path:

-op, --output-path
------------------
Relative or absolute path to the folder for saving outputs.

If the folder does not exist then it is created from the specified path.
The creation is applied recursively to the path, starting from the right
until the existing folder is reached.

.. code-block:: text

    default : current directory (".")
    type : str

See also: :ref:`example <scripts_output-notes>`.


.. _tb2j-extractor_output-name:

-on, --output-name
------------------
Seedname for the output files.

Output files will have the following name structure: *output-name*
If this parameter is not specified then result will be printed in 
standard output stream. If none is specify, output is passed to the console.

.. code-block:: text

    default : None
    type : str

See also: :ref:`example <scripts_output-notes>`.


.. _tb2j-extractor_decimals:

-d, --decimals
--------------
Decimals after the comma for the exchange values.

.. code-block:: text

    default : 4
    type : int

.. note::
    Changed in the version 0.5.17 from "-acc"/"--accuracy".

.. _tb2j-extractor_force-symmetry:

-fs, --force-symmetry
---------------------
Whenever to force the symmetry of the template on the model.

.. code-block:: text

    default : False
    type : bool


.. _tb2j-extractor_isotropic:

-i, --isotropic
---------------
Whenever to output isotropic exchange.

.. code-block:: text

    default : False
    type : bool

Section format:

.. code-block:: text

        Isotropic: J


.. _tb2j-extractor_anisotropic:

-a, --anisotropic
-----------------
Whenever to output anisotropic exchange.

.. code-block:: text

    default : False
    type : bool

Section format:

.. code-block:: text

        Anisotropic: 
            Jxx Jxy Jxz
            Jxy Jyy Jyz
            Jxz Jyz Jzz


.. _tb2j-extractor_matrix:

-m, --matrix
------------
Whenever to output whole matrix exchange.

.. code-block:: text

    default : False
    type : bool

Section format:

.. code-block:: text

        Matrix: 
            Jxx Jxy Jxz
            Jyx Jyy Jyz
            Jzx Jzy Jzz


.. _tb2j-extractor_dmi:

-dmi
----
Whenever to output DMI exchange.

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


.. _tb2j-extractor_all:

-all
----
Whenever to all types of exchange.

.. code-block:: text

    default : False
    type : bool


.. _tb2j-extractor_verbose:

-v, -verbose
------------
Verbose output, propagates to the called methods.

.. code-block:: text

    default : False

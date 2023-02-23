.. _tb2j-extractor:

*********************
``tb2j-extractor.py``
*********************

Script for extracting of template-based model from 
`TB2J <https://tb2j.readthedocs.io/en/latest/>`_ results.


.. _tb2j-extractor_verbose-ref:

Verbose output
==============
If :ref:`--verbose <tb2j-extractor_verbose>` and 
:ref:`--output-name <tb2j-extractor_output-name>` arguments are provided to 
the script then output will be written in .txt and .py files.
Content of the .py file is organised in a following manner:

.. code-block:: python

    import numpy as np 
    iso = {
        name:
        {
            (i, j, k): J_iso,
            ...
        },
        ...
    }

    aniso = {
        name:
        {
            (i, j, k): J_aniso,
            ...
        },
        ...
    }

    dmi = {
        name:
        {
            (i, j, k): DMI,
            ...
        },
        ...
    }

    matrix = {
        name:
        {
            (i, j, k): J_matrix,
            ...
        },
        ...
    }

where 

.. code-block:: python

    type(J_iso) = float

    J_aniso = np.array(
        [[J_aniso_xx, J_aniso_xy, J_aniso_xz],
         [J_aniso_xy, J_aniso_yy, J_aniso_yz],
         [J_aniso_xz, J_aniso_yz, J_aniso_zz]])

    DMI = np.array([D_x, D_y, D_z])

    J_matrix = np.array(
        [[J_xx, J_xy, J_xz],
         [J_yx, J_yy, J_yz],
         [J_zx, J_zy, J_zz]])

Content of .txt file is organised as: ::

    name
      atom1 atom2 (i, j, k)
        Isotropic: J_iso
        Anisotropic:
            J_aniso_xx J_aniso_xy J_aniso_xz
            J_aniso_xy J_aniso_yy J_aniso_yz
            J_aniso_xz J_aniso_yz J_aniso_zz
        DMI: D_x D_y D_z
        |DMI|: |DMI|
        |DMI/J|: |DMI|/|J|
        Matrix:
            J_xx J_xy J_xz
            J_yx J_yy J_yz
            J_zx J_zy J_zz
      
      ...
    
    ...

Arguments
=========

.. _tb2j-extractor_filename:

-f, --filename
--------------
Relative or absulute path to the *exchange.out* file,
including the name and extention of the file itself.

.. code-block:: text

    required

.. _tb2j-extractor_template-file:

-tf, --template-file
--------------------
Relative or absolute path to the template file, 
including the name and extention of the file.

.. code-block:: text

    required


See also: :ref:`template <rad-make-template>`


.. _tb2j-extractor_output-dir:

-op, --output-dir
-----------------
Relative or absolute path to the folder for saving outputs.

If the folder does not exist then it is created from the specified path.
The creation is applied recursevly to the path, starting from the right
until the existing folder is reached.

.. code-block:: text

    default : current directory

See also: :ref:`example <scripts_output-notes>`.


.. _tb2j-extractor_output-name:

-on, --output-name
------------------
Seedname for the output files.

Output files will have the following name structure: *output-name*
If this parameter is not specified then result will be printed in 
standart output stream.

.. code-block:: text

    default : None

See also: :ref:`example <scripts_output-notes>`.


.. _tb2j-extractor_dmi:

-dmi
----
Whenever to print each dmi vector for each exchange group separately.

.. code-block:: text

    default : False


.. _tb2j-extractor_verbose:

-v, --verbose
-------------
Whenever to print each neighbor from a 
:ref:`template file <tb2j-extractor_template-file>` in a verbose way.

.. code-block:: text

    default : False

.. _tb2j-extractor_accuracy:

-acc, --accuracy
----------------
Accuracy for the exchange values.

.. code-block:: text

    default : 4

.. _tb2j-extractor_force-symmetry:

-fs, --force-symmetry
---------------------
Whenever to force the symmetry of the template on the model.

.. code-block:: text

    default : True



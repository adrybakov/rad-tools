.. _tb2j_refractor:

*********************
``tb2j_refractor.py``
*********************

Script for extracting of template-based model from 
`TB2J <https://tb2j.readthedocs.io/en/latest/>`_ results.

Arguments
=========

.. _tb2j_refractor_filename:

``--filename``, ``-f``
----------------------
Relative or absulute path to the *exchange.out* file,
including the name and extention of the file itself.

    *required* : True

    *type* : str


.. _tb2j_refractor_template:

``--template``, ``-t``
----------------------
Relative or absolute path to the template file, 
including the name and extention of the file.

    *required* : True

    *type* : str

See also: :ref:`template <rad-make-template>`


.. _tb2j_refractor_output-dir:

``--output-dir``, ``-op``
-------------------------
Relative or absolute path to the folder for saving outputs.

If the folder does not exist then it is created from the specified path.
The creation is applied recursevly to the path, starting from the right
until the existing folder is reached.

    *default* : current directory
        
    *type* : str


.. _tb2j_refractor_output-name:

``--output-name``, ``-on``
--------------------------
Seedname for the output files.

Output files will have the following name structure: *output-name*

    *default* : exchange_refr

    *type* : str

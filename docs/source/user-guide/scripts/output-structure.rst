.. _output-notes:

****************
Output structure
****************

.. note::
    In version 0.8.0 the ``--output_path`` argument was removed.

.. hint::
    Starting form the version 0.8.0 all output is saved at the same folder where the
    input file is located by default. The ``--output-name`` argument is used to change
    the default behavior.

For all scripts there is a common logic for the ``--output-name`` argument:

* Script produces only one file.

The value of ``--output-name`` argument is used as a name of the output file. If the 
``--output-name`` is a path to the folder (i.e. "foo/bar/") then the output file has 
the default value, but saved in a provided folder. If the ``--output-name`` is a path
to the file (i.e. "foo/bar/file.txt") then the output file has the provided name and
saved in a provided folder.

* Script produces multiple files.

Same logic as for the one file, but the provided name is used as a seedname for the
output files. 

* Script produces a folder.

The provided name is used as a name of the folder. If the 
``--output-name`` is a path to the folder (i.e. "foo/bar/") then the output folder has 
the default value, but saved in a provided folder. If the ``--output-name`` is a path
to the file (i.e. "foo/bar/file") then the output folder has the provided name and
saved in a provided folder.

* Script produces multiple folders.

Logic is the same as for the one folder, but the provided name is used as a seedname for the
output folders.


.. hint::
    In the last case "foo/bar/" and "foo/bar" are two different values. The first one
    creates output folder with the default name in the folder "foo/bar/": "foo/bar/output-folder/", 
    the second one creates output folder with the name "bar" in the folder "foo": "foo/bar/".


The ``script.py`` script is executed from the 
folder "example" and the file structure is:

.. code-block:: text

    example/
    ├── input-file
    └── output/

.. code-block:: bash

   script.py -if input-file 

After the execution the "example" folder looks similar to:
    
.. code-block:: text

    example/
    ├── input-file
    ├── output-name-1
    ├── output-name-2
    └── output/

Script produced two output files "output-name-1"
and "output-name-2". Shared seedname "output-name" is different for each 
script and comes by default.

.. important::
    The output files are not located in "output" folder since the 
    current folder is used for output by default.
    
Next command saves the output in the "output" folder:

.. code-block:: bash

    script.py -if input-file -on output/

After it's execution "example" folder should have the structure:

.. code-block:: text

    example/
    ├── input-file
    ├── output-name-1
    ├── output-name-2
    └── output/
        ├── output-name-1
        └── output-name-2

Output files have the same names, but they are saved in the "output" 
folder as specified by ``--output-name`` (``-on``) argument.

It is not necessary to specify a path to the existing folder, 
for example after the command:

.. code-block:: bash

    script.py -if input-file -on output/bar/foo/

the script creates folder "bar" inside of the folder "output" and folder 
"foo" inside of the folder "bar". The structure of the "example" folder 
should look like:

.. code-block:: text
    
    example/
    ├── input-file
    ├── output-name-1
    ├── output-name-2
    └── output/
        ├── output-name-1
        ├── output-name-2
        └── bar/
            └── foo/
                ├── output-name-1
                └── output-name-2

In order to change the shared output name one may run:

.. code-block:: bash

    script.py -if input-file -on output/custom-output-name

The structure of the "example" folder now should be the following: 

.. code-block:: text
    
    example/
    ├── input-file
    ├── output-name-1
    ├── output-name-2
    └── output/
        ├── output-name-1
        ├── output-name-2
        ├── custom-output-name-1
        ├── custom-output-name-2
        └── bar/
            └── foo/
                ├── output-name-1
                └── output-name-2


.. _output-notes:

****************
Output structure
****************

For all scripts there is a common logic for the ``--output-path`` and 
``--output-name`` arguments. In some scripts only ``--output-name`` 
is available and then the logic is straightforward: there is only one 
output file and it is saved as specified in ``--output-name`` argument.

Meanwhile, in some scripts several output files are expected and additional 
argument ``--output-path`` appear, which specifies the directory where all
output files are saved and ``--output-name`` specifies common name for 
all output files. Here is an example for this case:


The ``script.py`` script is executed from the 
folder "example" and the file structure is:

.. code-block:: text

    example/
    ├── input-file
    └── output/

.. code-block:: bash

   script.py -f input-file 

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

    script.py -f input-file -op output

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
folder as specified by ``--output-dir`` argument.

It is not necessary to specify a path to the existing folder, 
for example after the command:

.. code-block:: bash

    script.py -f input-file -op output/bar/foo

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

    script.py -f input-file -op output -on custom-output-name

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


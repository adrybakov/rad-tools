.. _output-notes:

************
Output notes
************

For all scripts there is a common logic for the ``--output-path`` and 
``--output-name`` arguments. In some scripts only ``--output-name`` 
is available and then the logic is straightforward: there is only one 
output file and it will be saved as specified in ``--output-name`` argument.

Meanwhile, in some scripts several output files are expected and additional 
argument ``--output-path`` appear, which will specify the directory where all
output files are saved and ``--output-name`` specify common name for all output
file. Here is an example for this case:


Imagine you are executing the ``script.py`` script from the 
folder *example* and your file structure looks like the following

.. code-block:: text

    example/
    ├── input-file
    └── output/
        
Now lets run the script

.. code-block:: bash

   script.py -f input-file 

After the execution your *example* folder will look similar to this
    
.. code-block:: text

    example/
    ├── input-file
    ├── output-name-1
    ├── output-name-2
    └── output/

Script produced two output files *output-name-1*
and *output-name-2*. Common seedname *output-name* comes by default.

.. important::
    That output files are not located in *output* folder since the 
    current folder is used for output by default.
    
Lets save the output in the *output* folder:

.. code-block:: bash

    script.py -f input-file -op output

Now *example* folder should look like this

.. code-block:: text

    example/
    ├── input-file
    ├── output-name-1
    ├── output-name-2
    └── output/
        ├── output-name-1
        └── output-name-2

Output files have the same names, but they are saved in the *output* 
folder as your specified by ``--output-dir`` argument.

It is not necessary to specify a path to the existing folder, 
for example try to execute

.. code-block:: bash

    script.py -f input-file -op output/bar/foo

The script will create folder *bar* inside of the folder *output* and folder 
*foo* inside of the folder *bar*. The structure of the *example* folder now 
should look like that:

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

In order to change the common output name run the following:

.. code-block:: bash

    script.py -f input-file -op output -on custom-output-name

The structure of the *example* folder now should look like that: 

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


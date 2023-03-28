*******
Scripts
*******

For any script use --help or -h option in console in order to display 
help message with the short summary of the arguments.

.. code-block:: console

   script-name.py -h

All the files used in the usage examples are available :examples:`here <>`.


Scripts for analysing `TB2J <https://tb2j.readthedocs.io/en/latest/>`_ 
outputs.

.. toctree::
   :maxdepth: 1

   tb2j-plotter
   tb2j-extractor
   rad-make-template

Other scripts.

.. toctree::
   :maxdepth: 1

   rad-dos-plotter
   identify-wannier-centres

.. _scripts_output-notes:

Notes on the output
-------------------

For all scripts there is a common logic for the ``--output_dir`` and 
``--output_name`` arguments. In some scripts only ``--output_name`` 
is available and then the logic is straightforward: there is only one 
output file and it will be saved as specified in ``--output_name`` argument.

Meanwhile, in some scripts several output files are expected and additional 
argument ``--output_dir`` appear, which will specify the directory where all
output files are saved and ``--output_name`` specify common name for all output
file. Here is an example for this case:


Imagine you are executing the ``script.py`` script from the 
folder *example* and your file structure looks like the following

.. code-block:: text

    example/
    ├── input_file
    └── output/
        
Now lets run the script

.. code-block:: bash

   script.py -f input_file 

After the execution your *example* folder will look similar to this
    
.. code-block:: text

    example/
    ├── input_file
    ├── output_name_1
    ├── output_name_2
    └── output/

Script produced two output files *output_name_1*
and *output_name_2*. Common seedname *output_name* comes by default.

.. important::
    That output files are not located in *output* folder since the 
    current folder is used for output by default.
    
Lets save the output in the *output* folder:

.. code-block:: bash

    script.py -f input_file -op output

Now *example* folder should look like this

.. code-block:: text

    example/
    ├── input_file
    ├── output_name_1
    ├── output_name_2
    └── output/
        ├── output_name_1
        └── output_name_2

Output files have the same names, but they are saved in the *output* 
folder as your specified by ``--output-dir`` argument.

It is not necessary to specify a path to the existing folder, 
for example try to execute

.. code-block:: bash

    script.py -f input_file -op output/bar/foo

The script will create folder *bar* inside of the folder *output* and folder 
*foo* inside of the folder *bar*. The structure of the *example* folder now 
should look like that:

.. code-block:: text
    
    example/
    ├── input_file
    ├── output_name_1
    ├── output_name_2
    └── output/
        ├── output_name_1
        ├── output_name_2
        └── bar/
            └── foo/
                ├── output_name_1
                └── output_name_2

In order to change the common output name run the following:

.. code-block:: bash

    script.py -f input_file -op output -on custom_output_name

The structure of the *example* folder now should look like that: 

.. code-block:: text
    
    example/
    ├── input_file
    ├── output_name_1
    ├── output_name_2
    └── output/
        ├── output_name_1
        ├── output_name_2
        ├── custom_output_name_1
        ├── custom_output_name_2
        └── bar/
            └── foo/
                ├── output_name_1
                └── output_name_2


``tb2j_plotter.py``
===================
Script for visualisation of TB2J exchange.out file.
---------------------------------------------------

Allows the user to visualise the exchange.out file. 
Currently sorting by R vector, distance and template file is supported.

Parameters
----------

``--filename``, ``-f``
~~~~~~~~~~~~~~~~~~~~~~
   Relative or absulute path to the *exchange.out* file, 
   including the name and extention of the file itself.

   *required* : True

   *type* : str

``--output-dir``, ``-op``
~~~~~~~~~~~~~~~~~~~~~~~~~
Relative or absolute path to the directory for saving output.

If the directory does not exist the it will be created from the path.
The creation will be applied recursevly until the existent directory 
is reacheed.

*default* : current directory
        
*type* : str

``--what-to-plot``, ``-wtp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Type of data for display.

Specifying the data for display at the graph. 
Isotropic exchange parameters are displayed by default. 
Currently available for display: Isotropic exchange parameter, distance.

*default* : iso 

*type* : str

*choices* : iso, distance

``--draw-cells``, ``-dc``
~~~~~~~~~~~~~~~~~~~~~~~~~
Whenever to draw the supercell`s shape.

If specified then the shape of all supercells 
presented in the model after filtering will be drawn.

*default* : False

*action* : store_true

``--R-vector``, ``-R``
~~~~~~~~~~~~~~~~~~~~~~
R vectors for filtering the model.

*default* : None

*type* : int

*nargs* : *





    


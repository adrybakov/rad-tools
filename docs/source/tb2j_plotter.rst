``tb2j_plotter.py``
===================
Script for visualisation of TB2J results.

Display Isotropic exchange or distances, each value in a separate picture
by default. Currently filtering by R vectors, distances and template file 
is supported.

Parameters
----------

``--filename``, ``-f``

    Relative or absulute path to the TB2J exchange output file, 
    including the name and extention of the file.

        *required* : True

        *type* : str

``--output-dir``, ``-op``

    Relative or absolute path to the directory for saving outputs.

    If the directory does not exist then it is created from the specified path.
    The creation is applied recursevly to the path, starting from the right
    until the existent directory is reached.

        *default* : current directory
        
        *type* : str

``--output-name``, ``-on``

    Seedname for the output files.

    Output file will have the following name structure:
    *seedname*.*display_data_type*.png

        *default* : exchange
        
        *type* : str

``--what-to-plot``, ``-wtp``

    Type of data for display.

    Specifying the data for display at the graph. 
    Everything is displayed by default, each value in a separate picture. 
    Currently available for display: Isotropic exchange parameter, distance.

        *default* : all

        *type* : str

        *choices* : all, iso, distance

``--draw-cells``, ``-dc``

    Whenever to draw the supercell`s shape.

    If specified then the shape of all supercells 
    presented in the model after filtering is drawn.

        *default* : False

        *action* : store_true

``--R-vector``, ``-R``

    R vectors for filtering the model.

    In TB2J outputs the bond is defined by atom 1 (from) and atom 2 (to). 
    Atom 1 always located in (0, 0, 0) supercell, while atom 2 is located in 
    R = (i, j, k) supercell. This parameter tells to keep only the bonds
    for which atom 2 is located in one of specified R supercells. 
    In order to specify supercells provide a set of integers 
    separated by spaces. They will be grouped by three to form a set of R vectors 
    starting from the left. If the last group will contain 1 or 2 integers 
    they will be ignored.

        *default* : None

        *type* : int

        *nargs* : *

``--max-distance``, ``-maxd``

    (<=) Maximum distance.

    All the bonds with the distance beetwen atom 1 and atom 2 
    greater then maximum distance will be excluded from the model.

        *default* : None

        *type* : float

``--min-distance``, ``-mind``

    (>=) Minimum distance.

    All the bonds with the distance beetwen atom 1 and atom 2 
    lower then minimum distance will be excluded from the model.

        *default* : None

        *type* : float

``--distance``, ``-d``

    (=) Exact distance.

    Only the bonds with the distance exact as this one remains in the model.
    Note: there is no point in specifying maximum or minimum distance when 
    this parameter is specified.

``--template``, ``-t``

    Relative or absolute path to the template file, 
    including the name and extention of the file.

    #TODO

    *default* : None

    *type* : str

``--double-bonds``, ``-db``

    Whenever to keep both bonds.

    In TB2J file there is two bonds for the pair of atom 1 and atom 2: 
    from 1 to 2 and from 2 to 1 (when R = (0, 0, 0)). Isotropic and 
    anisotropic exchange and distance usially are exactly the same. 
    DMI vector have the same module and opposite directions. 
    If this parameter is specifyied then both bonds are displayed. 
    Otherwise bonds are combined in one by taking the average beetween
    exchange parameters (Note that it forces DMI to be equal to zero).

        *default* : False

        *action* : store_true
 

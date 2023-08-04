.. _guide_crystal_kpoints:

.. currentmodule:: radtools

********
K points
********

For the full reference see :ref:`api_kpoints`

:py:class:`Kpoints` is a relatively small class, which represent the set of k-points
of a particular path of a given lattice.

It is used to generate full kpoints path both for calculation and for plotting.


Creation
========

Usually it is created from some :py:class`.Lattice` (or :py:class`.Crystal`):

.. doctest::

    >>> from radtools import Kpoints, lattice_example
    >>> lattice = lattice_example("CUB")
    >>> kp = lattice.kpoints
    >>> kp.hs_names
    ['G', 'M', 'R', 'X']

However, it could be created explicitly:

.. doctest::

    >>> from radtools import Kpoints
    >>> b1, b2, b3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    >>> names = ["G", "X"]
    >>> coordinates = [[0,0,0], [0.5,0,0]]
    >>> labels = [R"$\Gamma$", "X"]
    >>> kp = Kpoints(b1, b2, b3, names=names, coordinates=coordinates, labels=labels)
    >>> kp.hs_names
    ['G', 'X']

Settings
========

After the creation you can add high symmetry kpoints, set the amount of kpoints between them and path: 

* :py:attr:`.Kpoints.n`

The amount of kpoints to be generated between each pair of high symmetry points in the path.

.. doctest::

    >>> # Default value is 100
    >>> kp.n
    100
    >>> kp.n = 10
    >>> kp.n
    10

* :py:attr:`.Kpoints.path`

The path itself. We use a specific format in the package: "G-K-X|R-S".
"-" separates high symmetry points in each subpath, "|" separates sections of the path.  
In the example n points are generated between "G" and "K", between "K" ans "X", 
between "R" and "S", but not between "X" and "R".
By default path is constructed from the list of points.

.. doctest::

    >>> from radtools import Kpoints
    >>> b1, b2, b3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    >>> names = ["G", "K", "X", "R"]
    >>> coordinates = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0], [0.5, 0.5, 0.5]]
    >>> labels = ["$\Gamma$", "K", "X", "R"]
    >>> kp = Kpoints(b1, b2, b3, names=names, coordinates=coordinates, labels=labels)
    >>> kp.path
    [['G', 'K', 'X', 'R']]
    >>> # It cause an Error, because high symmetry point "S" is not defined
    >>> kp.path = "G-K-X|R-S"
    Traceback (most recent call last):
    ...
    ValueError: Point 'S' is not defined. Defined points are:
      G : [0 0 0]
      K : [0.5 0.5 0. ]
      X : [0.5 0.  0. ]
      R : [0.5 0.5 0.5]
    >>> kp.path = "G-K-X|R-G"
    >>> kp.path
    [['G', 'K', 'X'], ['R', 'G']]
    >>> kp.add_hs_point(name="S", coordinates=[0.5, 0.5, 0.5], label="S")
    >>> kp.path = "G-K-X|R-S"
    >>> kp.path
    [['G', 'K', 'X'], ['R', 'S']]
    >>> kp.path_string
    'G-K-X|R-S'

.. note::

    Internally RAD-tools stores the path as a list of subpaths, where each subpath
    is a list of high symmetry points. This format is also correct for assigning the path attribute.

Usage
=====

Once the setting of the Kpoints are done, it can be used for calculation or plotting.

Calculation
-----------

There is one property suitable for calculation: :py:attr:`Kpoints.points`. which is an array 
of all generated kpoints. For each pair of high symmetry points it generates :py:attr:`Kpoints.n`
between them. The first and the last points are always the high symmetry points of this section of the path.

.. doctest::

    >>> from radtools import Kpoints
    >>> b1, b2, b3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    >>> names = ["G", "K", "X"]
    >>> coordinates = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0]]
    >>> labels = ["$\Gamma$", "K", "X"]
    >>> kp = Kpoints(b1, b2, b3, names=names, coordinates=coordinates, labels=labels, n=4)
    >>> kp.points()
    array([[0. , 0. , 0. ],
           [0.1, 0.1, 0. ],
           [0.2, 0.2, 0. ],
           [0.3, 0.3, 0. ],
           [0.4, 0.4, 0. ],
           [0.5, 0.5, 0. ],
           [0.5, 0.5, 0. ],
           [0.5, 0.4, 0. ],
           [0.5, 0.3, 0. ],
           [0.5, 0.2, 0. ],
           [0.5, 0.1, 0. ],
           [0.5, 0. , 0. ]])

.. note::
    For each section the last point is repeated twice, because it is the first point 
    of the next section of the path.    

    .. code-block:: python

        array([[0. , 0. , 0. ], # <--- Gamma
               [0.1, 0.1, 0. ],
               [0.2, 0.2, 0. ],
               [0.3, 0.3, 0. ],
               [0.4, 0.4, 0. ],
               [0.5, 0.5, 0. ], # <--- K
               [0.5, 0.5, 0. ], # <--- K
               [0.5, 0.4, 0. ],
               [0.5, 0.3, 0. ],
               [0.5, 0.2, 0. ],
               [0.5, 0.1, 0. ],
               [0.5, 0. , 0. ]]) # <--- X

Plotting
========

For plotting there are three properties. Two of them are for the high symmetry points
and describe the labels and position of ticks on the x-axis:

.. doctest::

    >>> kp.labels
    ['$\\Gamma$', 'K', 'X']
    >>> import numpy as np
    >>> np.around(kp.coordinates(), decimals=4)    
    array([0.    , 0.7071, 1.2071])

The third property gives the coordinates of the :py:attr:`.Kpoints.points` for the plot:

.. doctest::

    >>> for point in kp.flatten_points():
    ...     print(round(point, 4))
    ... 
    0.0
    0.1414
    0.2828
    0.4243
    0.5657
    0.7071
    0.7071
    0.8071
    0.9071
    1.0071
    1.1071
    1.2071

.. note::
    Those coordinates are directly corresponds to the k-points from the previous subsection. 

    .. code-block:: python
 
        0.0    # <--- Gamma
        0.1414
        0.2828
        0.4243
        0.5657
        0.7071 # <--- K
        0.7071 # <--- K
        0.8071
        0.9071
        1.0071
        1.1071
        1.2071 # <--- X

.. hint::

    Repeated :py:attr:`.Kpoints.points` or :py:attr:`.Kpoints.flatten_points` 
    can be used to restore the position of high symmetry points in the path.
    



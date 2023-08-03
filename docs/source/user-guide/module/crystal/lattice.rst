.. _guide_crystal_lattice:

.. currentmodule:: radtools

*******
Lattice
*******

For the full reference see :ref:`api_lattice`

This guide describe attributes and methods of the general lattice class.
Each Bravais lattice class is a child of :py:class:`.Lattice`, 
thus all it`s method and attribute are available. For the description of each 
Bravais lattice class see:

.. toctree::
    :maxdepth: 2
    
    bravais-lattices/index

Creation
========

Lattice can be created in three different ways:

* From ``cell`` matrix:

.. doctest::

    >>> from radtools import Lattice
    >>> cell = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    >>> lattice = Lattice(cell)
    >>> lattice.cell
    array([[1, 0, 0],
           [0, 1, 0],
           [0, 0, 1]])

* From three lattice vectors :math:`\vec{a}_1`, :math:`\vec{a}_2`, :math:`\vec{a}_3`:

.. doctest::
    
        >>> from radtools import Lattice
        >>> a1 = [1, 0, 0]
        >>> a2 = [0, 1, 0]
        >>> a3 = [0, 0, 1]
        >>> lattice = Lattice(a1, a2, a3)
        >>> lattice.cell
        array([[1, 0, 0],
               [0, 1, 0],
               [0, 0, 1]])

* From lattice parameters :math:`a`, :math:`b`, :math:`c`, :math:`\alpha`, :math:`\beta`, :math:`\gamma`:

.. doctest::
    
        >>> from radtools import Lattice
        >>> lattice = Lattice(1, 1, 1, 90, 90, 90)
        >>> import numpy as np
        >>> np.round(lattice.cell, decimals=1)
        array([[1., 0., 0.],
               [0., 1., 0.],
               [0., 0., 1.]])

Identification of the lattice
=============================

Bravais lattice type can be calculated via :py:meth:`.Lattice.identify` method:

.. doctest::
    
        >>> from radtools import Lattice
        >>> lattice = Lattice(1, 1, 1, 90, 90, 90)
        >>> lattice.identify()
        'CUB'

Identification procedure is realised in the :py:func:`.lepage` function. 
For the algorithm description and reference see :ref:`rad-tools_lepage`.

:py:attr:`.Lattice.identify` method returns a string with the lattice type, 
but does not change the lattice itself. Therefore, if you want to get an instance of 
:py:class:`.CUB` lattice, you should create it manually:

.. doctest::
    
        >>> from radtools import Lattice
        >>> lattice = Lattice(1, 1, 1, 90, 90, 90)
        >>> lattice.identify()
        'CUB'
        >>> type(lattice)
        <class 'radtools.crystal.lattice.Lattice'>
        >>> from radtools import bravais_lattice_from_cell
        >>> lattice = bravais_lattice_from_cell(lattice.cell)
        >>> type(lattice)
        <class 'radtools.crystal.bravais_lattice.CUB'>


.. note::
    
    Lattice identification is not a trivial task. 
    The algorithm is based on the assumption that the lattice`s unit cell is primitive.


Reference attributes
====================

If lattice is an instance of one of the Bravais lattice classes,
then you can use the following attributes for the information about the lattice:

.. doctest::

    >>> from radtools import bravais_lattice_from_cell
    >>> lattice = bravais_lattice_from_cell([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    >>> lattice.pearson_symbol
    'cP'
    >>> lattice.crystal_family
    'c'
    >>> lattice.centring_type
    'P'

Lattice parameters
==================

All lattice parameters can be accessed as attributes:

* Real space

.. doctest::

    >>> from radtools import bravais_lattice_from_cell
    >>> lattice = bravais_lattice_from_cell([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    >>> lattice.a1  
    array([1, 0, 0])
    >>> lattice.a2
    array([0, 1, 0])
    >>> lattice.a3
    array([0, 0, 1])
    >>> lattice.cell
    array([[1, 0, 0],
           [0, 1, 0],
           [0, 0, 1]])
    >>> lattice.a
    1.0
    >>> lattice.b
    1.0
    >>> lattice.c
    1.0
    >>> lattice.alpha
    90.0
    >>> lattice.beta
    90.0
    >>> lattice.gamma
    90.0
    >>> lattice.unit_cell_volume
    1
    >>> lattice.parameters
    (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)

* Reciprocal space

.. doctest::

    >>> from radtools import bravais_lattice_from_cell
    >>> from math import pi
    >>> lattice = bravais_lattice_from_cell([[2*pi, 0, 0], [0, 2*pi, 0], [0, 0, 2*pi]])
    >>> lattice.b1
    array([1., 0., 0.])
    >>> lattice.b2
    array([0., 1., 0.])
    >>> lattice.b3
    array([0., 0., 1.])
    >>> lattice.reciprocal_cell
    array([[1., 0., 0.],
           [0., 1., 0.],
           [0., 0., 1.]])
    >>> lattice.k_a
    1.0
    >>> lattice.k_b
    1.0
    >>> lattice.k_c
    1.0
    >>> lattice.k_alpha
    90.0
    >>> lattice.k_beta
    90.0
    >>> lattice.k_gamma
    90.0

Variation of the lattice
========================
Some of the Bravais lattice classes have several variations.

To check the variation of the lattice use :py:attr:`.Lattice.variation` attribute:

.. doctest::

    >>> from radtools import lattice_example, Lattice
    >>> lattice = lattice_example("BCT")
    >>> lattice.variation
    'BCT1'
    >>> lattice = Lattice(1, 1, 1, 90, 90, 90)
    >>> lattice.variation
    'Lattice'

Plotting of the lattice
=======================

Lattice primitive, conventional unit cells as  well as the Wigner-Seitz cell 
and Brillouin zone can be plotted using :py:meth:`.Lattice.plot` method:

.. doctest::

    >>> from radtools import lattice_example
    >>> lattice = lattice_example("BCT")
    >>> lattice.plot("brillouin")
    >>> lattice.plot("kpath")

To show the plot use :py:meth:`.Lattice.show` method:

.. doctest::

    >>> lattice.show() # doctest: +SKIP

More examples of the plots with the code snippets can be found in the 
:ref:`guide_crystal_bravais-lattices` guide.

Full list of the available plotting methods can be found in the
:ref:`api_lattice` reference.

K points
========

Path in reciprocal space and k points for plotting and calculation are partially 
implemented in a separate class :py:class:`.Kpoints`.

High symmetry points, path are the part of the :py:class:`.Lattice` class, but once 
they are set an instance of :py:class:`.Kpoints` class is expected to be created in order 
to access labels of the plot, flatten coordinates, whole list of k points, etc.

.. doctest::

    >>> from radtools import Lattice
    >>> lattice = Lattice(1, 1, 1, 90, 90, 90)
    >>> lattice.add_kpoint("Gamma", [0, 0, 0])
    >>> lattice.add_kpoint("X", [0.5, 0, 0])
    >>> lattice.add_kpoint("M", [0.5, 0.5, 0])
    >>> lattice.add_kpoint("CP", [0.5, 0.5, 0.5], plot_name="Custom plot name")
    >>> lattice.path = "Gamma-X|M-CP-X"
    >>> # n = 100 by default, it could be changed on the kp instance.
    >>> kp = lattice.get_kpoints(n=100)

.. note::

    For each Bravais lattice class there is a predefined path and set of 
    kpoints in reciprocal space. See :ref:`guide_crystal_bravais-lattices` for more details.

For the full guide on how to use :py:class:`.Kpoints` class see :ref:`guide_crystal_kpoints`.












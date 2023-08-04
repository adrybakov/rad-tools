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

Bravais lattice type is lazily identified when it is needed:

.. doctest::
    
        >>> from radtools import Lattice
        >>> lattice = Lattice(1, 1, 1, 90, 90, 90)
        >>> lattice.type()
        'CUB'

Identification procedure is realised in the :py:func:`.lepage` function. 
For the algorithm description and reference see :ref:`rad-tools_lepage`.

.. note::
    
    Lattice identification is not a trivial task. 
    The algorithm is based on the assumption that the lattice`s unit cell is primitive.


Reference attributes
====================

You can use the following attributes for the information about the lattice 
based on the Bravais type:

.. doctest::

    >>> from radtools import Lattice
    >>> lattice = Lattice([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
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

    >>> from radtools import Lattice
    >>> lattice = Lattice([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
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
    1.0
    >>> lattice.parameters
    (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)

* Reciprocal space

.. doctest::

    >>> from radtools import Lattice
    >>> from math import pi
    >>> lattice = Lattice([[2*pi, 0, 0], [0, 2*pi, 0], [0, 0, 2*pi]])
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
    >>> print(f"{lattice.k_a:.4f}")
    1.0000
    >>> print(f"{lattice.k_b:.4f}")
    1.0000
    >>> print(f"{lattice.k_c:.4f}")
    1.0000
    >>> lattice.k_alpha
    90.0
    >>> lattice.k_beta
    90.0
    >>> lattice.k_gamma
    90.0

Variation of the lattice
========================
Some of the Bravais lattice types have several variations.

To check the variation of the lattice use :py:attr:`.Lattice.variation` attribute:

.. doctest::

    >>> from radtools import lattice_example, Lattice
    >>> lattice = lattice_example("BCT")
    >>> lattice.variation
    'BCT1'
    >>> lattice = Lattice(1, 1, 1, 90, 90, 90)
    >>> lattice.variation
    'CUB'

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

Path in reciprocal space and k points for plotting and calculation are 
implemented in a separate class :py:class:`.Kpoints`. It is expected to be accessed 
through the :py:attr:`.Lattice.kpoints` attribute. Note that it is store in the lattice itself,
so the you can work with kpoints through the lattice instance or through the kpoints instance:

.. doctest::

    >>> from radtools import Lattice
    >>> lattice = Lattice(1, 1, 1, 90, 90, 90)
    >>> lattice.kpoints.add_hs_point("CP", [0.5, 0.5, 0.5], label="Custom label")
    >>> lattice.kpoints.path = "G-X|M-CP-X"
    >>> lattice.kpoints.path_string
    'G-X|M-CP-X'
    >>> kp = lattice.kpoints
    >>> kp.path_string
    'G-X|M-CP-X'
    >>> kp.path = "G-X|M-X"
    >>> kp.path_string
    'G-X|M-X'
    >>> lattice.kpoints.path_string
    'G-X|M-X'

.. note::

    For each Bravais lattice type there is a predefined path and set of 
    kpoints in reciprocal space. See :ref:`guide_crystal_bravais-lattices` for more details.

For the full guide on how to use :py:class:`.Kpoints` class see :ref:`guide_crystal_kpoints`.












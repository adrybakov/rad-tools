.. _guide_crystal_lattice:

.. currentmodule:: radtools

*******
Lattice
*******

For the full reference see :ref:`api_lattice`

Every Bravais lattice is an instance of the :py:class:`.Lattice` class.
For the guide about Bravais lattices see :ref:`guide_crystal_bravais-lattices`.
This page describes the :py:class:`.Lattice` class and its methods.

Import 
======

    >>> # Exact import
    >>> from radtools.crystal.lattice import Lattice
    >>> # Explicit import
    >>> from radtools.crystal import Lattice
    >>> # Recommended import
    >>> from radtools import Lattice

For the examples in this page we need additional import and some predefined variables:

.. doctest::

    >>> from radtools import lattice_example

Creation
========

Lattice can be created in three different ways:

* From ``cell`` matrix:

.. doctest::

    >>> cell = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    >>> lattice = Lattice(cell)
    >>> lattice.cell
    array([[1, 0, 0],
           [0, 1, 0],
           [0, 0, 1]])

When a lattice created from the cell orientation of the cell is respected,
however the lattice vectors may be renamed. 
See :ref:`library_lattice-standardization` for details.

Creation may change the angles and the lengths of the cell vectors.
It preserve the volume, right- or left- handedness, lattice type and variation
of the cell.

The lattice vector`s lengths are preserved as a set.

The angles between the lattice vectors are preserved as a set with possible
changes of the form: :math:`angle \rightarrow 180 - angle`.

The returned cell may not be the same as the input one, but it is translational
equivalent.

* From three lattice vectors :math:`\vec{a}_1`, :math:`\vec{a}_2`, :math:`\vec{a}_3`:

.. doctest::
    
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
    
        >>> lattice = Lattice(1, 1, 1, 90, 90, 90)
        >>> lattice.type()
        'CUB'

Identification procedure is implemented in the :py:func:`.lepage` function. 
For the algorithm description and reference see :ref:`library_lepage`.

.. note::
    
    Lattice identification is not a trivial task and may be time consuming.
    The algorithm is based on the assumption that the lattice`s unit cell is primitive.

    By default the lattice type is identified during the creation of the lattice 
    (It is required for the lattice standardization). Therefore, the creation of the
    lattice may be time consuming. To avoid this, you can disable the standardization
    of the cell via the ``standardize=False`` argument:

    .. doctest::
    
        >>> lattice = Lattice(1, 1, 1, 90, 90, 90, standardize=False)

    Note that the predefined paths and k points for the lattice are not guaranteed to 
    be correct and reproducible if the lattice is not standardized.

Variation of the lattice
========================
Some of the Bravais lattice types have several variations.

To check the variation of the lattice use :py:attr:`.Lattice.variation` attribute:

.. doctest::

    >>> lattice = lattice_example("BCT")
    >>> lattice.variation
    'BCT1'
    >>> lattice = Lattice(1, 1, 1, 90, 90, 90)
    >>> lattice.variation
    'CUB'

Reference attributes
====================

You can use the following attributes for the information about the lattice 
based on the Bravais type:

.. doctest::

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
    >>> round(lattice.k_a, 4)
    1.0000
    >>> round(lattice.k_b, 4)
    1.0000
    >>> round(lattice.k_c, 4)
    1.0000
    >>> lattice.k_alpha
    90.0
    >>> lattice.k_beta
    90.0
    >>> lattice.k_gamma
    90.0

.. hint::

    Not all properties of the lattice are listed here (for examples the once 
    for the conventional cell are not even mentioned). 
    See :ref:`api_lattice` for the full list of properties.

Plotting of the lattice
=======================

Lattice primitive, conventional unit cells as  well as the Wigner-Seitz cell 
and Brillouin zone can be plotted using :py:meth:`.Lattice.plot` method:

.. doctest::

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
through the :py:attr:`.Lattice.kpoints` attribute. Note that you can work with 
kpoints from the instance of the :py:class:`.Lattice`, since the instance of the 
:py:class:`.Kpoints` class is created when the property is accessed for the first 
time and stored internally for the future:

.. doctest::

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
    The unit cell has to be standardized to use the predefined paths and kpoints.

For the full guide on how to use :py:class:`.Kpoints` class see :ref:`guide_crystal_kpoints`.












.. _guide_crystal_bravais-lattices:

****************
Bravais lattices
****************

.. currentmodule:: radtools

For the full reference see :ref:`api_crystal`

For the description of each Bravais lattice type see :ref:`library_bravais-lattices`.

Bravais lattice notation and standardization follows Setyawan and Curtarolo [1]_.

Each Bravais lattice is an instance of :py:class:`.Lattice` class.


For each Bravais lattice system there is a function defined, which constructs
the instance of :py:class:`.Lattice` class from the parameters. For the names of the 
constructors and corresponding parameters see the :ref:`tables below <library_bravais-lattices>` 
(for full reference see :ref:`api_bravais-lattices`). Before the main table we present
an example of the usage of the constructor for the cubic lattice.

Import
======

.. doctest::

    >>> # Exact import
    >>> from radtools.crystal.bravais_lattice.constructor import CUB
    >>> # Explicit import
    >>> from radtools.crystal import CUB
    >>> # Recommended import
    >>> from radtools import CUB

Creation
========

.. doctest::

    >>> lattice = CUB(1)
    >>> lattice.parameters
    (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)

Constructor can be used to get the cell instead of the lattice:

    >>> cell = CUB(1, return_cell=True)
    >>> cell
    array([[1, 0, 0],
           [0, 1, 0],
           [0, 0, 1]])

Predefined examples
===================

For each type and variation a predefined example of the lattice is available. 
It could be accessed in a following way:

.. doctest::

    >>> import radtools as rad
    >>> cubic_example = rad.lattice_example("cub")
    >>> cubic_example.variation
    'CUB'

.. hint::

    Capitalization of the name of the lattice example is not important:
    ``CUB``, ``cub`` and ``Cub`` are equivalent.


References
==========
.. [1] Setyawan, W. and Curtarolo, S., 2010. 
    High-throughput electronic band structure calculations: Challenges and tools. 
    Computational materials science, 49(2), pp.299-312. 

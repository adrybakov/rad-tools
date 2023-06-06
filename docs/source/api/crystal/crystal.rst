.. _crystal-api:

*******
Crystal
*******

.. currentmodule:: radtools

.. autosummary::
    :toctree: generated/

    Crystal

Properties
==========

.. hint:: 
    All properties and methods of :ref:`lattice-api` are accessible directly from the crystal instance.
    The following interfaces are identical:

    .. doctest::

        >>> import radtools as rad
        >>> c = rad.Crystal()
        >>> c.a
        1.0
        >>> c.lattice.a
        1.0

.. autosummary::

    Crystal.lattice

Atom methods
============

.. autosummary::

    Crystal.add_atom
    Crystal.remove_atom

Positioning
===========

.. autosummary::

    Crystal.get_atom_coordinates
    Crystal.get_distance
    Crystal.get_vector

Defining the type
=================

.. autosummary::

    Crystal.identify
    Crystal.find_primitive_cell
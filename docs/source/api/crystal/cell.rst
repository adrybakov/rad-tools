.. _api_cell:

*****
Cell*
*****

.. currentmodule:: radtools

.. note::
    Cell is not a class, but a submodule,
    which follows the style of the class with the @classmethods.
    It does not create the cell object, but group the functions,
    which corresponds to the cell logic. It is intended to be used as:

    .. doctest::

        >>> from radtools import Cell
        >>> Cell.params([[1,0,0],[0,1,0],[0,0,1]])
        (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)

    which is equivalent to:

    .. doctest::

        >>> from radtools.crystal.cell import params
        >>> params([[1,0,0],[0,1,0],[0,0,1]])
        (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)

Functions
=========

.. autosummary::
    :toctree: generated/

    Cell.reciprocal
    Cell.from_params
    Cell.params
    Cell.primitive

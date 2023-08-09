.. _guide_magnons_diagonalization:

***************
Diagonalization
***************

For the full reference see :py:func:`.solve_via_colpa`

.. currentmodule:: radtools

There are two main methods for the diagonalization of bosonic Hamiltonian.
In RAD-tools we implemented one from [1]_.

Import
======

.. doctest::

    >>> # Exact import
    >>> from radtools.magnons.diagonalization import solve_via_colpa
    >>> # Explicit import
    >>> from radtools.magnons import solve_via_colpa
    >>> # Recommended import
    >>> from radtools import solve_via_colpa

Usage
=====

One need to pass grand-dynamical matrix to receive the eigenvalues and
transformation matrix.

.. doctest::

    >>> solve_via_colpa([[1, 0], [0, 1]])
    (array([1., 1.]), array([[ 1., -0.],
           [-0.,  1.]]))

References
==========
.. [1] Colpa, J.H.P., 1978.
    Diagonalization of the quadratic boson hamiltonian.
    Physica A: Statistical Mechanics and its Applications,
    93(3-4), pp.327-353.
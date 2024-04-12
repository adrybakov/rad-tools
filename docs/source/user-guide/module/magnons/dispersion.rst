.. _guide_magnons_dispersion:

****************
MagnonDispersion
****************

For the full reference see :ref:`api_magnon-dispersion`.

.. currentmodule:: radtools

:py:class:`.MagnonDispersion` is a class for calculating magnon dispersion.
See the :ref:`article in library <library_magnon-dispersion-method>` for more information about the method behind the calculations.

Import
======

.. doctest::

    >>> # Exact import
    >>> from radtools.magnons.dispersion import MagnonDispersion
    >>> # Explicit import
    >>> from radtools.magnons import MagnonDispersion
    >>> # Recommended import
    >>> from radtools import MagnonDispersion

For the examples in this page we need additional import and some predefined variables:

.. doctest::

    >>> from radtools import SpinHamiltonian
    >>> # Ferromagnetic cubic chain
    >>> hamiltonian = SpinHamiltonian(notation="SpinW")
    >>> hamiltonian.add_atom("Fe", position=(0, 0, 0), spin=[0, 0, 1])
    >>> hamiltonian.add_bond("Fe", "Fe", (1, 0, 0), iso=-1)
    >>> hamiltonian.add_bond("Fe", "Fe", (0, 1, 0), iso=-1)
    >>> hamiltonian.add_bond("Fe", "Fe", (0, 0, 1), iso=-1)
    >>> kp = hamiltonian.kpoints
    >>> kp.n = 0

Creation
========

Magnon dispersion is build from some :ref:`guide_spinham`.

.. note::

    :py:class:`Atom`\ 's of the spin Hamiltonian have to have :py:attr:`.Atom.spin_vector`
    attribute defined.

.. doctest::

    >>> dispersion = MagnonDispersion(hamiltonian)


Calculation
===========

To compute the dispersion simply call the dispersion object with the desired k-points.
K-points are an instance of :py:class:`Kpoints` class or list of k-points, where each
k-point is a list of three numbers.

.. doctest::

    >>> kp.path_string
    'G-X-M-G-R-X|M-R'
    >>> dispersion(kp)
    array([[ 0.,  8.,  8., 16., 16.,  0.,  0., 24., 24.,  8., 16., 24.]])

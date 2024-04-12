.. _guide_io_tb2j:

.. currentmodule:: radtools

*****************
|TB2J|_ interface
*****************

Main output of TB2J is the "exchange.out" file with the exchange parameters.

Function :py:func:`load_tb2j_model` reads this file
and construct :py:class:`SpinHamiltonian` from it.

Import
======

.. doctest::

    >>> # Exact import
    >>> from radtools.io.tb2j import load_tb2j_model
    >>> # Explicit import
    >>> from radtools.io import load_tb2j_model
    >>> # Recommended import
    >>> from radtools import load_tb2j_model

Usage
=====

.. doctest::

    >>> hamiltonian = load_tb2j_model('exchange.out')  # doctest: +SKIP

.. _api_hamiltonian:

*******************
ExchangeHamiltonian
*******************

.. currentmodule:: radtools
 
Class
=====

.. autosummary::
    :toctree: generated/

    ExchangeHamiltonian

Energy
======

.. autosummary::

    ExchangeHamiltonian.ferromagnetic_energy

Structure
=========

.. hint:: 
    All properties and methods of :ref:`api_crystal` are accessible directly from the hamiltonian instance.
    As a consequence all properties and methods of :ref:`api_lattice` are accessible directly from the hamiltonian instance.
    Unless they are defined directly for the :py:class:`.ExchangeHamiltonian` and listed here.

.. autosummary::

    ExchangeHamiltonian.crystal
    ExchangeHamiltonian.cell_list
    ExchangeHamiltonian.number_spins_in_unit_cell   
    ExchangeHamiltonian.space_dimensions

Manipulation with the model
===========================

Adding elements
---------------

.. autosummary::

    ExchangeHamiltonian.add_atom
    ExchangeHamiltonian.add_bond

Removing elements
-----------------

.. autosummary::

    ExchangeHamiltonian.remove_atom
    ExchangeHamiltonian.remove_bond

Filtering the model
-------------------

.. autosummary::

    ExchangeHamiltonian.filter
    ExchangeHamiltonian.filtered
    ExchangeHamiltonian.force_symmetry
    ExchangeHamiltonian.forced_symmetry

Notation
========

.. autosummary::

    ExchangeHamiltonian.notation
    ExchangeHamiltonian.set_interpretation

Individual properties
---------------------

.. autosummary::

    ExchangeHamiltonian.double_counting
    ExchangeHamiltonian.spin_normalized
    ExchangeHamiltonian.factor_one_half
    ExchangeHamiltonian.factor_two
    ExchangeHamiltonian.minus_sign

Deprecated parameters
=====================

Usage of this parameters is deprecated. Usage of corresponding 
:ref:`api_crystal` parameters is encouraged. Table of corresponding parameters:

==================== =======================
Deprecated parameter Corresponding parameter
==================== =======================
get_distance         crystal.get_distance
get_bond_vector      crystal.get_vector
get_atom_coordinates crystal.get_atom_coordinates
==================== =======================


.. autosummary::

    ExchangeHamiltonian.cell
    ExchangeHamiltonian.a
    ExchangeHamiltonian.b
    ExchangeHamiltonian.c
    ExchangeHamiltonian.len_a
    ExchangeHamiltonian.len_b
    ExchangeHamiltonian.len_c
    ExchangeHamiltonian.b1
    ExchangeHamiltonian.b2
    ExchangeHamiltonian.b3
    ExchangeHamiltonian.unit_cell_volume
    ExchangeHamiltonian.get_distance
    ExchangeHamiltonian.get_bond_vector
    ExchangeHamiltonian.get_atom_coordinates
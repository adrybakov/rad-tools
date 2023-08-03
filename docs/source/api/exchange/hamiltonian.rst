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

.. hint:: 
    All properties and methods of :ref:`api_crystal` are accessible directly 
    from the hamiltonian instance. As a consequence all properties and methods 
    of :ref:`api_lattice` are accessible directly from the hamiltonian instance.

Energy
======

.. autosummary::
    :toctree: generated/

    ExchangeHamiltonian.ferromagnetic_energy

Structure
=========

.. autosummary::
    :toctree: generated/

    ExchangeHamiltonian.cell_list
    ExchangeHamiltonian.number_spins_in_unit_cell   
    ExchangeHamiltonian.space_dimensions

Manipulation with the model
===========================

Adding elements
---------------

.. autosummary::
    :toctree: generated/

    ExchangeHamiltonian.add_atom
    ExchangeHamiltonian.add_bond

Removing elements
-----------------

.. autosummary::
    :toctree: generated/

    ExchangeHamiltonian.remove_atom
    ExchangeHamiltonian.remove_bond

Filtering the model
-------------------

.. autosummary::
    :toctree: generated/

    ExchangeHamiltonian.filter
    ExchangeHamiltonian.filtered
    ExchangeHamiltonian.force_symmetry
    ExchangeHamiltonian.forced_symmetry

Notation
========

.. autosummary::
    :toctree: generated/

    ExchangeHamiltonian.notation
    ExchangeHamiltonian.set_interpretation

Individual properties
---------------------

.. autosummary::
    :toctree: generated/

    ExchangeHamiltonian.double_counting
    ExchangeHamiltonian.spin_normalized
    ExchangeHamiltonian.factor_one_half
    ExchangeHamiltonian.factor_two
    ExchangeHamiltonian.minus_sign

Saving and loading
==================

.. autosummary::
    :toctree: generated/

    ExchangeHamiltonian.dump_txt
    ExchangeHamiltonian.dump_pickle


Crystal getter
==============

.. autosummary::
    :toctree: generated/

    ExchangeHamiltonian.crystal

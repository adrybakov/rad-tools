.. _api_spin-hamiltonian:

*******************
SpinHamiltonian
*******************

.. currentmodule:: radtools
 
Class
=====

.. autosummary::
    :toctree: generated/

    SpinHamiltonian

.. hint:: 
    All properties and methods of :ref:`api_crystal` are accessible directly 
    from the hamiltonian instance. As a consequence all properties and methods 
    of :ref:`api_lattice` are accessible directly from the hamiltonian instance.

Energy
======

.. autosummary::
    :toctree: generated/

    SpinHamiltonian.ferromagnetic_energy

Structure
=========

.. autosummary::
    :toctree: generated/

    SpinHamiltonian.cell_list
    SpinHamiltonian.number_spins_in_unit_cell   
    SpinHamiltonian.space_dimensions

Manipulation with the model
===========================

Adding elements
---------------

.. autosummary::
    :toctree: generated/

    SpinHamiltonian.add_atom
    SpinHamiltonian.add_bond

Removing elements
-----------------

.. autosummary::
    :toctree: generated/

    SpinHamiltonian.remove_atom
    SpinHamiltonian.remove_bond

Filtering the model
-------------------

.. autosummary::
    :toctree: generated/

    SpinHamiltonian.filter
    SpinHamiltonian.filtered
    SpinHamiltonian.form_model
    SpinHamiltonian.forced_symmetry

Notation
========

.. autosummary::
    :toctree: generated/

    SpinHamiltonian.notation
    SpinHamiltonian.set_interpretation

Individual properties
---------------------

.. autosummary::
    :toctree: generated/

    SpinHamiltonian.double_counting
    SpinHamiltonian.spin_normalized
    SpinHamiltonian.factor_one_half
    SpinHamiltonian.factor_two
    SpinHamiltonian.minus_sign

Saving and loading
==================

.. autosummary::
    :toctree: generated/

    SpinHamiltonian.dump_txt
    SpinHamiltonian.dump_pickle


Crystal getter
==============

.. autosummary::
    :toctree: generated/

    SpinHamiltonian.crystal

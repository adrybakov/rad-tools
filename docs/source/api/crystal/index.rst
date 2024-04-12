.. _api_crystal-module:

.. currentmodule:: radtools

*******
crystal
*******

.. versionadded:: 0.7

Classes
=======

.. toctree::
    :maxdepth: 1

    crystal
    lattice
    atom
    kpoints
    cell

.. _api_crystal-plotting:

Plotting
========

.. toctree::
    :maxdepth: 1

    plotly-backend
    matplotlib-backend

Identification routines
=======================

.. autosummary::
    :toctree: generated/

    niggli
    lepage

Properties
==========

.. autosummary::
    :toctree: generated/

    dipole_dipole_energy
    dipole_dipole_interaction

.. _api_bravais-lattices:

Bravais lattices
================

Examples
--------

.. autosummary::
    :toctree: generated/

    lattice_example

Constructors
------------

.. autosummary::
    :toctree: generated/

    CUB
    FCC
    BCC
    TET
    BCT
    ORC
    ORCF
    ORCI
    ORCC
    HEX
    RHL
    MCL
    MCLC
    TRI

High-symmetry k points
----------------------

.. autosummary::
    :toctree: generated/

    CUB_hs_points
    FCC_hs_points
    BCC_hs_points
    TET_hs_points
    BCT_hs_points
    ORC_hs_points
    ORCF_hs_points
    ORCI_hs_points
    ORCC_hs_points
    HEX_hs_points
    RHL_hs_points
    MCL_hs_points
    MCLC_hs_points
    TRI_hs_points

Cell`s standardization
----------------------

.. autosummary::
    :toctree: generated/

    standardize_cell
    CUB_standardize_cell
    FCC_standardize_cell
    BCC_standardize_cell
    TET_standardize_cell
    BCT_standardize_cell
    ORC_standardize_cell
    ORCF_standardize_cell
    ORCI_standardize_cell
    ORCC_standardize_cell
    HEX_standardize_cell
    RHL_standardize_cell
    MCL_standardize_cell
    MCLC_standardize_cell
    TRI_standardize_cell

Lattice variations
------------------

.. autosummary::
    :toctree: generated/

    BCT_variation
    ORCF_variation
    RHL_variation
    MCLC_variation
    TRI_variation

Constants
=========

.. currentmodule:: radtools.crystal

.. autosummary::
    :toctree: generated/

    ABS_TOL
    REL_TOL
    MIN_LENGTH
    MAX_LENGTH
    ABS_TOL_ANGLE
    REL_TOL_ANGLE
    MIN_ANGLE
    ATOM_TYPES
    PEARSON_SYMBOLS
    BRAVAIS_LATTICE_NAMES
    BRAVAIS_LATTICE_VARIATIONS
    TRANSFORM_TO_CONVENTIONAL
    DEFAULT_K_PATHS
    HS_PLOT_NAMES

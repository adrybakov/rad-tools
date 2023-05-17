.. _rad-tools_lattice:

*******
Lattice
*******

.. automodule:: rad_tools.crystal.lattice

General 3D lattice
==================
.. currentmodule:: rad_tools

General lattice is describe by the class :py:class:`.Lattice`.

.. autosummary::
    :toctree: generated/

    Lattice


Bravais lattices
================

For each type of Bravais lattice a class defined, for some classes there are 
several variants of lattice, each of which are treated under same class 
(see :py:attr:`variant` for each class).

For each type and variant a predefined example of the lattice is available. It could be accessed in a following way:

.. doctest::

    >>> import rad_tools as rad_tools
    >>> cubic_example = rad.cub



Cubic lattice system
--------------------

Pre-defined examples: ``cub``, ``fcc``, ``bcc``.

.. autosummary::
    :toctree: generated/

    CUB
    FCC
    BCC

Tetragonal lattice system
-------------------------

Pre-defined examples: ``tet``, ``bct1``, ``bct2``.

.. autosummary::
    :toctree: generated/

    TET
    BCT

Orthorhombic lattice system
---------------------------

Pre-defined examples: ``orc``, ``orcf1``, ``orcf2``, ``orcf3``, ``orci``, ``orcc``.

.. autosummary::
    :toctree: generated/
    
    ORC
    ORCF
    ORCI
    ORCC

Hexagonal lattice system
------------------------

Pre-defined examples: ``hex``.

.. autosummary::
    :toctree: generated/
    
    HEX

Rhombohedral lattice system
---------------------------

Pre-defined examples: ``rhl1``, ``rhl2``

.. autosummary::
    :toctree: generated/
    
    RHL

Monoclinic lattice system
-------------------------

Pre-defined examples: ``mcl``, ``mclc1``, ``mclc2``, ``mclc3``, ``mclc4``, ``mclc5``.

.. autosummary::
    :toctree: generated/

    MCL
    MCLC

Triclinic lattice system
------------------------

Pre-defined examples: ``tri1a``, ``tri1b``, ``tri2a``, ``tri2b``.

.. autosummary::
    :toctree: generated/

    TRI

Construction of the Brillouin zone 
==================================

Wigner-Seitz cell is constructed in a same way.

.. automodule:: rad_tools.crystal.decomposition



Functions
---------

.. currentmodule:: rad_tools.crystal

.. autosummary::
    :toctree: generated/

    deduct_zone
    get_lattice_points
    define_planes
    define_corners
    define_edges






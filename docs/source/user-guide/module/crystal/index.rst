.. _rad-tools_crystal:

*******
Crystal
*******

.. versionadded:: 0.7

.. currentmodule:: radtools

Crystals are defined as a combination of the the :py:class:`.Lattice` and the set of  
:py:class:`.Atom`\ s.

.. toctree::
    :maxdepth: 2
    
    lattice
    atom


When arbitrary crystal structure is considered the :py:attr:`.Crystal.lattice` is an instance of 
the general :py:class:`.Lattice` class, where unit cell is interpreted as a primitive one. 
The type of the Bravais lattice is not defined as well.
The process of the definition of the crystal type is divided in two steps: 

* Define primitive cell of the :py:class:`.Crystal`
* Define the type of the Bravais lattice.

After those two steps primitive cel is well defined and the standard k path is accessible.

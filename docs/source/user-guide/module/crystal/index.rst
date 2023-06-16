.. _rad-tools_crystal:

*******
Crystal
*******

.. versionadded:: 0.7

For the full reference see :ref:`api_crystal-module`

.. currentmodule:: radtools

Crystal is a base class, which store the structure. 
It is defined as a combination of the the :py:class:`.Lattice` and the set of  
:py:class:`.Atom`\ s.

.. toctree::
    :maxdepth: 2
    
    lattice
    atom


When arbitrary crystal structure is considered the :py:attr:`.Crystal.lattice` is an instance of 
the general :py:class:`.Lattice` class, where unit cell is interpreted as a primitive one
and the type of the Bravais lattice is not defined.
The process of the definition of the crystal type is divided in two steps: 

* Define primitive cell of the :py:class:`.Crystal`
    :py:meth:`.Crystal.find_primitive_cell`
    Currently not implemented, any cell is interpreted as primitive.
* Define the type of the Bravais lattice.
    :py:meth:`.Crystal.identify`

This routine is implemented in :py:meth:`.Crystal.identify` method. It calls for the 
:py:meth:`.Crystal.find_primitive_cell` internally.



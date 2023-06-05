.. _rad-tools_lattice:

*******
Lattice
*******

.. automodule:: radtools.crystal.lattice

General 3D lattice
==================
.. currentmodule:: radtools

General lattice is describe by the class :py:class:`.Lattice`. 


Bravais lattices
================

Bravais lattice notation follows Setyawan and Curtarolo [1]_.

For each type of Bravais lattice a class defined, for some classes there are 
several variations of lattice, each of which are treated under same class 
(see :py:attr:`.Lattice.variation`).

For each type and variation a predefined example of the lattice is available. 
It could be accessed in a following way:

.. doctest::

    >>> import radtools as rad
    >>> cubic_example = rad.lattice_example("cub")
    >>> cubic_example.variation
    'CUB'

Each Bravais Lattice is created by the parameters 
:math:`a`, :math:`b`, :math:`c`, :math:`\alpha`, :math:`\beta`, :math:`\gamma`,
which corresponds to the conventional lattice. However, attributes of the class 
``a``, ``b``, ``c``, ``alpha``, ``beta``, ``gamma`` 
return parameters of the primitive lattice. Conventional lattice may be accessed through
the attributes ``conv_cell``, ``conv_a``, ``conv_b``, ``conv_c``, 
``conv_alpha``, ``conv_beta``, ``conv_gamma`` .


.. currentmodule:: radtools

Creation
--------
Bravais lattice can be created by the instantiation of corresponding class.
For the predefined bravais lattice examples one may use 
:py:func:`.lattice_example` function. 

From the set of lattice parameters (:math:`a`, :math:`c`, :math:`c`,
:math:`\alpha`, :math:`\beta`, :math:`\gamma`) or unit ``cell`` 
(interpreted as primitive cell) one may create the bravais lattice with 
:py:func:`.bravais_lattice_from_param` and :py:func:`bravais_lattice_from_cell`.

Cubic lattice system
--------------------
Predefined examples: ``cub``, ``fcc``, ``bcc``.

.. toctree::
    :maxdepth: 1
    
    bravais-lattices/cub/index
    bravais-lattices/fcc/index
    bravais-lattices/bcc/index

Tetragonal lattice system
-------------------------

Predefined examples: ``tet``, ``bct1``, ``bct2``.

.. toctree::
    :maxdepth: 1
    
    bravais-lattices/tet/index
    bravais-lattices/bct/index

Orthorhombic lattice system
---------------------------

Predefined examples: ``orc``, ``orcf1``, ``orcf2``, ``orcf3``, ``orci``, ``orcc``.

.. toctree::
    :maxdepth: 1
    
    bravais-lattices/orc/index
    bravais-lattices/orcf/index
    bravais-lattices/orci/index
    bravais-lattices/orcc/index

Hexagonal lattice system
------------------------

Predefined examples: ``hex``.

.. toctree::
    :maxdepth: 1
    
    hex

Rhombohedral lattice system
---------------------------

Predefined examples: ``rhl1``, ``rhl2``

.. toctree::
    :maxdepth: 1
    
    rhl

Monoclinic lattice system
-------------------------

Predefined examples: ``mcl``, ``mclc1``, ``mclc2``, ``mclc3``, ``mclc4``, ``mclc5``.

.. toctree::
    :maxdepth: 1

    mcl
    mclc

Triclinic lattice system
------------------------

Predefined examples: ``tri1a``, ``tri1b``, ``tri2a``, ``tri2b``.

.. toctree::
    :maxdepth: 1

    tri


References
==========
.. [1] Setyawan, W. and Curtarolo, S., 2010. 
    High-throughput electronic band structure calculations: Challenges and tools. 
    Computational materials science, 49(2), pp.299-312. 



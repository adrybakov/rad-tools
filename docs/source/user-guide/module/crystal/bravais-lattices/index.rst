.. _guide_crystal_bravais-lattices:

****************
Bravais lattices
****************

.. currentmodule:: radtools

For the full reference see :ref:`api_crystal`

Bravais lattice notation and standardization follows Setyawan and Curtarolo [1]_.

Each Bravais lattice is an instance of :py:class:`.Lattice` class.


For each Bravais lattice system there is a function defined, which constructs
the instance of :py:class:`.Lattice` class from the parameters. For the names of the 
constructors and corresponding parameters see the :ref:`tables below <table_bravais-lattices>` 
(for full reference see :ref:`api_bravais-lattices`). Before the main table we present
an example of the usage of the constructor for the cubic lattice.

Import
======

.. doctest::

    >>> # Exact import
    >>> from radtools.crystal.bravais_lattice.constructor import CUB
    >>> # Explicit import
    >>> from radtools.crystal import CUB
    >>> # Recommended import
    >>> from radtools import CUB

Creation
========

.. doctest::

    >>> lattice = CUB(1)
    >>> lattice.parameters
    (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)

Constructor can be used to get the cell instead of the lattice:

    >>> cell = CUB(1, return_cell=True)
    >>> cell
    array([[1, 0, 0],
           [0, 1, 0],
           [0, 0, 1]])

Predefined examples
===================

For each type and variation a predefined example of the lattice is available. 
It could be accessed in a following way:

.. doctest::

    >>> import radtools as rad
    >>> cubic_example = rad.lattice_example("cub")
    >>> cubic_example.variation
    'CUB'

.. hint::

    Capitalization of the name of the lattice example is not important:
    ``CUB``, ``cub`` and ``Cub`` are equivalent.

.. _table_bravais-lattices:

Individual docs
===============

For the documentation of each Bravais lattice system see the following pages:

Cubic lattice system
--------------------

.. toctree::
    :hidden:
    
    cub/index
    fcc/index
    bcc/index

=================  ==========  ===============  ================
Name               Examples    Parameters       Constructor
=================  ==========  ===============  ================
:ref:`guide_cub`   ``cub``     :math:`a`        :py:func:`.CUB`
:ref:`guide_fcc`   ``fcc``     :math:`a`        :py:func:`.FCC`
:ref:`guide_bcc`   ``bcc``     :math:`a`        :py:func:`.BCC`
=================  ==========  ===============  ================

Tetragonal lattice system
-------------------------

.. toctree::
    :hidden:
    
    tet/index
    bct/index

=================  ==========  ===============  ================
Name               Examples    Parameters       Constructor
=================  ==========  ===============  ================
:ref:`guide_tet`   ``tet``     :math:`a`,       :py:func:`.TET`         
                               :math:`c`   
:ref:`guide_bct`   ``bct``,    :math:`a`,       :py:func:`.BCT`
                   ``bct1``,   :math:`c`
                   ``bct2`` 
=================  ==========  ===============  ================

Orthorhombic lattice system
---------------------------

.. toctree::
    :hidden:
    
    orc/index
    orcf/index
    orci/index
    orcc/index

=================  ==========  ===============  ================
Name               Examples    Parameters       Constructor
=================  ==========  ===============  ================
:ref:`guide_orc`   ``orc``     :math:`a`,       :py:func:`.ORC`         
                               :math:`b`,
                               :math:`c`   
:ref:`guide_orcf`  ``orcf``,   :math:`a`,       :py:func:`.ORCF`         
                   ``orcf1``,  :math:`b`,
                   ``orcf2``,  :math:`c`
                   ``orcf3``   
:ref:`guide_orci`  ``orci``    :math:`a`,       :py:func:`.ORCI`         
                               :math:`b`,
                               :math:`c`  
:ref:`guide_orcc`  ``orcc``    :math:`a`,       :py:func:`.ORCC`         
                               :math:`b`,
                               :math:`c`  
=================  ==========  ===============  ================

Hexagonal lattice system
------------------------

.. toctree::
    :hidden:
    
    hex/index

=================  ==========  ===============  ================
Name               Examples    Parameters       Constructor
=================  ==========  ===============  ================
:ref:`guide_hex`   ``hex``     :math:`a`,       :py:func:`.HEX`         
                               :math:`c`   
=================  ==========  ===============  ================

Rhombohedral lattice system
---------------------------

.. toctree::
    :hidden:
    
    rhl/index

=================  ==========  ===============  ================
Name               Examples    Parameters       Constructor
=================  ==========  ===============  ================
:ref:`guide_rhl`   ``rhl``,    :math:`a`,       :py:func:`.RHL`         
                   ``rhl1``,   :math:`c`
                   ``rhl2``    
=================  ==========  ===============  ================

Monoclinic lattice system
-------------------------

.. toctree::
    :hidden:

    mcl/index
    mclc/index

=================  ==========  ===============  ================
Name               Examples    Parameters       Constructor
=================  ==========  ===============  ================
:ref:`guide_mcl`   ``mcl``     :math:`a`,       :py:func:`.MCL`         
                               :math:`b`,
                               :math:`c`,  
                               :math:`\alpha` 
:ref:`guide_mclc`  ``mclc``,   :math:`a`,       :py:func:`.MCLC`         
                   ``mclc1``,  :math:`b`,
                   ``mclc2``,  :math:`c`,
                   ``mclc3``,  :math:`\alpha`
                   ``mclc4``,
                   ``mclc5``   
=================  ==========  ===============  ================

Triclinic lattice system
------------------------

Predefined examples: ``tri1a``, ``tri1b``, ``tri2a``, ``tri2b``.

.. toctree::
    :hidden:

    tri/index

=================  ==========  ===============  ================
Name               Examples    Parameters       Constructor
=================  ==========  ===============  ================
:ref:`guide_tri`   ``tri1a``,  :math:`a`,       :py:func:`.TRI`         
                   ``tri1b``,  :math:`b`,
                   ``tri2a``,  :math:`c`,  
                   ``tri2b``   :math:`\alpha`, 
                               :math:`\beta`,
                               :math:`\gamma`
=================  ==========  ===============  ================

References
==========
.. [1] Setyawan, W. and Curtarolo, S., 2010. 
    High-throughput electronic band structure calculations: Challenges and tools. 
    Computational materials science, 49(2), pp.299-312. 

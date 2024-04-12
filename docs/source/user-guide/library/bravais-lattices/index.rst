.. _library_bravais-lattices:

****************
Bravais lattices
****************

.. currentmodule:: radtools

For the full reference see :ref:`api_crystal`

For the module guide on Bravais lattices see :ref:`guide_crystal_bravais-lattices`

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

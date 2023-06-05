.. _lattice-tet:

********************
Tetragonal (TET, tP)
********************

Tetragonal lattice is described by the class :py:class:`.TET`.

It is defined by two parameters: :math:`a` and :math:`c` 
with primitive and conventional lattice:

.. math::

    \boldsymbol{a}_1 = (a, 0, 0)

    \boldsymbol{a}_2 = (0, a, 0)

    \boldsymbol{a}_3 = (0, 0, c)

Variations
==========

There are no variations for tetragonal lattice. 
One example is predefined: ``tet`` with :math:`a = \pi` and :math:`c = 1.5\pi`

Example structure
=================

**Default kpath**: :math:`\Gamma-X-M-\Gamma-Z-R-A-Z\vert X-R\vert M-A`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: tet_brillouin.png 
            :target: ../../../../../_images/tet_brillouin.png 
      - .. literalinclude:: tet_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: tet_real.png 
            :target: ../../../../../_images/tet_real.png 
      - .. literalinclude:: tet_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: tet_wigner-seitz.png 
            :target: ../../../../../_images/tet_wigner-seitz.png 
      - .. literalinclude:: tet_wigner-seitz.py
            :language: py


Edge cases
==========

If :math:`a = c`, then the lattice is :ref:`lattice-cub`.
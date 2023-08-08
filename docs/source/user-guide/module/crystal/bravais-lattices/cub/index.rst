.. _guide_cub:

***********
Cubic (CUB)
***********

**Pearson symbol**: cP

**Constructor**:  :py:func:`.CUB`.

It is defined by one parameter: :math:`a` with primitive and conventional lattice:

.. math::

    \boldsymbol{a}_1 = (a, 0, 0)

    \boldsymbol{a}_2 = (0, a, 0)

    \boldsymbol{a}_3 = (0, 0, a)

Variations
==========

There are no variations for cubic lattice. 
One example is predefined: ``cub`` with :math:`a = \pi`.

Example structure
=================

**Default kpath**: :math:`\mathrm{\Gamma-X-M-\Gamma-R-X\vert M-R}`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: cub_brillouin.png 
            :target: ../../../../../_images/cub_brillouin.png 
      - .. literalinclude:: cub_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: cub_real.png 
            :target: ../../../../../_images/cub_real.png 
      - .. literalinclude:: cub_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: cub_wigner-seitz.png 
            :target: ../../../../../_images/cub_wigner-seitz.png 
      - .. literalinclude:: cub_wigner-seitz.py
            :language: py

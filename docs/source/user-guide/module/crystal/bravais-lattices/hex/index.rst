.. _guide_hex:

***************
Hexagonal (HEX)
***************

**Pearson symbol**: hP

**Constructor**:  :py:func:`.HEX`.

It is defined by two parameter: :math:`a` and :math:`c` 
with primitive and conventional lattice:

.. math::

    \boldsymbol{a}_1 = (a/2, -a\sqrt{3}, 0)

    \boldsymbol{a}_2 = (a/2, a\sqrt{3}, 0)

    \boldsymbol{a}_3 = (0, 0, c)

Variations
==========

There are no variations for hexagonal lattice. 
One example is predefined: ``hex`` with :math:`a = \pi` and :math:`c = 2\pi`.

Example structure
=================

**Default kpath**: :math:`\mathrm{\Gamma-M-K-\Gamma-A-L-H-A\vert L-M\vert K-H}`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: hex_brillouin.png 
            :target: ../../../../../_images/hex_brillouin.png 
      - .. literalinclude:: hex_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: hex_real.png 
            :target: ../../../../../_images/hex_real.png 
      - .. literalinclude:: hex_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: hex_wigner-seitz.png 
            :target: ../../../../../_images/hex_wigner-seitz.png 
      - .. literalinclude:: hex_wigner-seitz.py
            :language: py

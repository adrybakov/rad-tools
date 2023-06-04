.. _lattice-orc:

**********************
Orthorhombic (ORC, oP)
**********************

Orthorombic lattice is described by the class :py:class:`.ORC`.


It is defined by three parameter: :math:`a`, :math:`b` and :math:`c` 
with primitive and conventional lattice:

.. math::

    \boldsymbol{a}_1 = (a, 0, 0)

    \boldsymbol{a}_2 = (0, b, 0)

    \boldsymbol{a}_3 = (0, 0, c)

Order of parameters: :math:`a < b < c`

Variations
==========

There are no variations for orthorhombic lattice. 
One example is predefined: ``orc`` with 
:math:`a = \pi`, :math:`b  = 1.5\pi` and :math:`c = 2\pi`.

Example structure
=================

**Default kpath**: :math:`\Gamma-X-S-Y-\Gamma-Z-U-R-T-Z\vert Y-T\vert U-X\vert S-R`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orc_brillouin.png 
            :target: ../../../../../_images/orc_brillouin.png 
      - .. literalinclude:: orc_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orc_real.png 
            :target: ../../../../../_images/orc_real.png 
      - .. literalinclude:: orc_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orc_wigner-seitz.png 
            :target: ../../../../../_images/orc_wigner-seitz.png 
      - .. literalinclude:: orc_wigner-seitz.py
            :language: py


Ordering of lattice parameters
==============================
TODO
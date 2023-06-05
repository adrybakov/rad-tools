.. _lattice-orcc:

************************************
Base-centred orthorhombic (ORCC, oS)
************************************

Base-centered orthorombic lattice is described by the class :py:class:`.ORCC`.


It is defined by three parameter: :math:`a`, :math:`b` and :math:`c` 
with conventional lattice:

.. math::

    \boldsymbol{a}_1 = (a, 0, 0)

    \boldsymbol{a}_2 = (0, b, 0)

    \boldsymbol{a}_3 = (0, 0, c)

And primitive lattice:

.. math::

    \begin{matrix}
        &\boldsymbol{a}_1 = (a/2, -b/2, 0) \\
        &\boldsymbol{a}_2 = (a/2, b/2, 0) \\
        &\boldsymbol{a}_3 = (0, 0, c) \\
    \end{matrix}

Order of parameters: :math:`a < b`

Variations
==========

There are no variations for base-centered orthorombic. 
One example is predefined: ``orcc`` with 
:math:`a = \pi`, :math:`b  = 1.3\pi` and :math:`c = 1.7\pi`.

Example structure
=================

**Default kpath**: :math:`\Gamma-X-S-R-A-Z-\Gamma-Y-X_1-A_1-T-Y\vert Z-T`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orcc_brillouin.png 
            :target: ../../../../../_images/orcc_brillouin.png 
      - .. literalinclude:: orcc_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orcc_real.png 
            :target: ../../../../../_images/orcc_real.png 
      - .. literalinclude:: orcc_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orcc_wigner-seitz.png 
            :target: ../../../../../_images/orcc_wigner-seitz.png 
      - .. literalinclude:: orcc_wigner-seitz.py
            :language: py


Ordering of lattice parameters
==============================
TODO

Edge cases
==========
If :math:`a = b`, then the lattice is :ref:`lattice-tet`.

If :math:`a = b = \sqrt{2} c`, then the lattice is :ref:`lattice-cub`.


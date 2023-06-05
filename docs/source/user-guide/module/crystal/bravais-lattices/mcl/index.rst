.. _lattice-mcl:

****************
Monoclinic (MCL)
****************

**Pearson symbol**: mP

Monoclinic lattice is described by the class :py:class:`.MCL`.

It is defined by four parameter: :math:`a`, :math:`b`, :math:`c` and :math:`\alpha` 
with primitive and conventional lattice:

.. math::


    \boldsymbol{a}_1 = (a, 0, 0)

    \boldsymbol{a}_2 = (0, b, 0)

    \boldsymbol{a}_3 = (0, c\cos(\alpha), c\sin(\alpha))


Order of parameters: :math:`a \le c`, :math:`b \le c`, :math:`\alpha < 90^{\circ}`.


Variations
==========

There are two variations for monoclinic lattice. 
One example is predefined: ``mcl`` with 
MCL(pi, 1.3 * pi, 1.6 * pi, alpha=75)
:math:`a = \pi`, :math:`b = 1.3 \pi` :math:`c = 1.6 \pi` and :math:`\alpha = 75^{\circ}`.


Example structure
=================


**Default kpath**: :math:`\Gamma-Y-H-C-E-M_1-A-X-H_1\vert M-D-Z\vert Y-D`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mcl_brillouin.png 
            :target: ../../../../../_images/mcl_brillouin.png 
      - .. literalinclude:: mcl_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mcl_real.png 
            :target: ../../../../../_images/mcl_real.png 
      - .. literalinclude:: mcl_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mcl_wigner-seitz.png 
            :target: ../../../../../_images/mcl_wigner-seitz.png 
      - .. literalinclude:: mcl_wigner-seitz.py
            :language: py

Ordering of lattice parameters
==============================
TODO

Edge cases
==========

If (:math:`\alpha = 60^{\circ}` or :math:`\alpha = 120^{\circ}`) and :math:`b = c`, 
then the lattice is :ref:`lattice-hex`.

If (:math:`\alpha = 30^{\circ}` or :math:`\alpha = 150^{\circ}`
or :math:`\alpha = 45^{\circ}` or :math:`\alpha = 145^{\circ}`) and :math:`b = c`, 
then the lattice is :ref:`lattice-orcc`.

If (:math:`\alpha = 60^{\circ}` or :math:`\alpha = 120^{\circ}`) and :math:`a \ne b = c/2`, 
then the lattice is :ref:`lattice-orc`.

If :math:`a \ne b \ne c` and :math:`\alpha = 90^{\circ}`, then the lattice is :ref:`lattice-orc`.

If (:math:`\alpha = 60^{\circ}` or :math:`\alpha = 120^{\circ}`) and :math:`a = b = c/2`, 
then the lattice is :ref:`lattice-tet`.

If (:math:`a = b \ne c` or :math:`a = c \ne b` or :math:`b = c \ne a`) and :math:`\alpha = 90^{\circ}`, then the lattice is :ref:`lattice-tet`.

If :math:`a = b = c` and :math:`\alpha = 90^{\circ}`, then the lattice is :ref:`lattice-cub`.

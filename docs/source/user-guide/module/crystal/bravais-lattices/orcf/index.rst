.. _guide_orcf:

********************************
Face-centred orthorhombic (ORCF)
********************************

**Pearson symbol**: oF

Face-centered orthorombic lattice is described by the class :py:class:`.ORCF`.

It is defined by three parameters: :math:`a`, :math:`b` and :math:`c` 
with conventional lattice:

.. math::

        \boldsymbol{a}_1 = (a, 0, 0)

        \boldsymbol{a}_2 = (0, b, 0)

        \boldsymbol{a}_3 = (0, 0, c)

And primitive lattice:

.. math::

    \boldsymbol{a}_1 = (0, b/2, c/2)

    \boldsymbol{a}_2 = (a/2, 0, c/2)

    \boldsymbol{a}_3 = (a/2, b/2, 0)

Variations
==========

There are tree  variations of face-centered orthorombic lattice. 

For the examples of variations
:math:`a` is set to :math:`1`; :math:`b` and :math:`c` fulfil the conditions:

* :math:`b = \dfrac{c}{\sqrt{c^2 - 1}}`

* :math:`c > \sqrt{2}`


First condition defines in ORCF\ :sub:`3` lattice and ensures
ordering of lattice parameters :math:`b > a`. 
Ordering :math:`c > b` is forced by second condition. 

For ORCF\ :sub:`1` and ORCF\ :sub:`2` lattices :math:`a < 1` and :math:`a > 1` is chosen.
While :math:`b` and :math:`c` are the same as for ORCF\ :sub:`3` lattice.

At the end all three parameters are multiplied by :math:`\pi`.

ORCF\ :sub:`1`
--------------

:math:`\dfrac{1}{a^2} > \dfrac{1}{b^2} + \dfrac{1}{c^2}`.

Predefined example: ``orcf1`` with 
:math:`a = 0.7\pi`, :math:`b = 5\pi/4` and :math:`c = 5\pi/3`.

ORCF\ :sub:`2`
--------------

:math:`\dfrac{1}{a^2} < \dfrac{1}{b^2} + \dfrac{1}{c^2}`.

Predefined example: ``orcf2`` with 
:math:`a = 1.2\pi`, :math:`b = 5\pi/4` and :math:`c = 5\pi/3`.

ORCF\ :sub:`3`
--------------

:math:`\dfrac{1}{a^2} = \dfrac{1}{b^2} + \dfrac{1}{c^2}`.

Predefined example: ``orcf3`` with 
:math:`a = \pi`, :math:`b = 5\pi/4` and :math:`c = 5\pi/3`.

Example structures
==================

ORCF\ :sub:`1`
--------------

**Default kpath**: :math:`\Gamma-Y-T-Z-\Gamma-X-A_1-Y\vert T-X_1\vert X-A-Z\vert L-\Gamma`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orcf1_brillouin.png 
            :target: ../../../../../_images/orcf1_brillouin.png 
      - .. literalinclude:: orcf1_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orcf1_real.png 
            :target: ../../../../../_images/orcf1_real.png 
      - .. literalinclude:: orcf1_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orcf1_wigner-seitz.png 
            :target: ../../../../../_images/orcf1_wigner-seitz.png 
      - .. literalinclude:: orcf1_wigner-seitz.py
            :language: py

ORCF\ :sub:`2`
--------------

**Default kpath**: :math:`\Gamma-Y-C-D-X-\Gamma-Z-D_1-H-C\vert C_1-Z\vert X-H_1\vert H-Y\vert L-\Gamma`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orcf2_brillouin.png 
            :target: ../../../../../_images/orcf2_brillouin.png 
      - .. literalinclude:: orcf2_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orcf2_real.png 
            :target: ../../../../../_images/orcf2_real.png 
      - .. literalinclude:: orcf2_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orcf2_wigner-seitz.png 
            :target: ../../../../../_images/orcf2_wigner-seitz.png 
      - .. literalinclude:: orcf2_wigner-seitz.py
            :language: py


ORCF\ :sub:`3`
--------------

**Default kpath**: :math:`\Gamma-Y-T-Z-\Gamma-X-A_1-Y\vert X-A-Z\vert L-\Gamma`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orcf3_brillouin.png 
            :target: ../../../../../_images/orcf3_brillouin.png 
      - .. literalinclude:: orcf3_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orcf3_real.png 
            :target: ../../../../../_images/orcf3_real.png 
      - .. literalinclude:: orcf3_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orcf3_wigner-seitz.png 
            :target: ../../../../../_images/orcf3_wigner-seitz.png 
      - .. literalinclude:: orcf3_wigner-seitz.py
            :language: py

Ordering of lattice parameters
==============================
TODO

Edge cases
==========
If :math:`a = b \ne c` or :math:`a = c \ne b` or :math:`b = c \ne a`, 
then the lattice is :ref:`guide_bct`.

If :math:`a = b = c`, then the lattice is :ref:`guide_fcc`.

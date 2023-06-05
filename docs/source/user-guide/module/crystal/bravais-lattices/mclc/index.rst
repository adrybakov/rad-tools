.. _lattice-mclc:

******************************
Base-centred monoclinic (MCLC)
******************************

**Pearson symbol**: mS

Base-centered monoclinic lattice is described by the class :py:class:`.MCLC`.

It is defined by four parameter: :math:`a`, :math:`b`, :math:`c` and :math:`\alpha` 
with conventional lattice:

.. math::

    \boldsymbol{a}_1 = (a, 0, 0)

    \boldsymbol{a}_2 = (0, b, 0)

    \boldsymbol{a}_3 = (0, c\cos(\alpha), c\sin(\alpha))

And primitive lattice:

.. math::

    \boldsymbol{a}_1 = (a/2, b/2, 0)

    \boldsymbol{a}_2 = (-a/2, b/2, 0)

    \boldsymbol{a}_3 = (0, c\cos(\alpha), c\sin(\alpha))

Order of parameters: :math:`a \le c`, :math:`b \le c`, :math:`\alpha < 90^{\circ}`.

Variations
==========

There are five variations for base-centered monoclinic lattice.

Reciprocal :math:`\gamma` (:math:`k_{\gamma}`) is defined by the equation (for primitive lattice):

.. math::

    \cos(k_{\gamma}) = \frac{a^2 - b^2\sin^2(\alpha)}{a^2 + b^2\sin^2(\alpha)}

For MCLC\ :sub:`2` :math:`k_{\gamma} = 90`, therefore :math:`a = b \sin(\alpha)`. 
For MCLC\ :sub:`1` we choose :math:`a < b \sin(\alpha)` and 
for MCLC\ :sub:`3`, MCLC\ :sub:`4` and MCLC\ :sub:`5` we choose :math:`a > b \sin(\alpha)`.

For the variations 3-5 we define :math:`a = xb\sin(\alpha)`, where :math:`x > 1`.

Then the condition for MCLC\ :sub:`4` gives:

.. math::

    c = \frac{x^2}{x^2 - 1}b\cos(\alpha)

Where :math:`\cos(\alpha) > 0` (:math:`\alpha < 90^{\circ}`), since :math:`x > 1`.

And the ordering condition :math:`b \le c` gives:

.. math::

    \cos(\alpha) \ge \frac{x^2 - 1}{x^2}

For MCLC\ :sub:`3` (MCLC\ :sub:`5`) we choose parameters in a same way as for MCLC\ :sub:`4`, 
but with :math:`c > \frac{x^2}{x^2 - 1}b\cos(\alpha)` (:math:`c < \frac{x^2}{x^2 - 1}b\cos(\alpha)`)


MCLC\ :sub:`1`
--------------

:math:`k_{\gamma} > 90^{\circ}`,

Predefined example: ``mclc1`` with :math:`a = \pi`, :math:`b = 1.4\cdot\pi`, :math:`c = 1.7\cdot\pi` and :math:`\alpha = 80^{\circ}` 

MCLC\ :sub:`2`
--------------

:math:`k_{\gamma} = 90^{\circ}`,

Predefined example: ``mclc2`` with :math:`a = 1.4\cdot\pi\cdot\sin(75^{\circ})`, :math:`b = 1.4\cdot\pi`, :math:`c = 1.7\cdot\pi` and :math:`\alpha=75^{\circ}` 

MCLC\ :sub:`3`
--------------

:math:`k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} < 1`

Predefined example with :math:`b = \pi`, :math:`x = 1.1`, :math:`\alpha = 78^{\circ}`, which produce:

``mclc4`` with :math:`a = 1.1\cdot\sin(78)\cdot\pi`, :math:`b = \pi`, 
:math:`c = 1.8\cdot 121\cdot\cos(65)\cdot\pi/21` and :math:`\alpha = 78^{\circ}` 

MCLC\ :sub:`4`
--------------

:math:`k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} = 1`

Predefined example with :math:`b = \pi`, :math:`x = 1.2`, :math:`\alpha = 65^{\circ}`, which produce:

``mclc4`` with :math:`a = 1.2\sin(65)\pi`, :math:`b = \pi`, 
:math:`c = 36\cos(65)\pi/11` and :math:`\alpha = 65^{\circ}` 

MCLC\ :sub:`5`
--------------

:math:`k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} > 1`

Predefined example with :math:`b = \pi`, :math:`x = 1.4`, :math:`\alpha = 53^{\circ}`, which produce:

``mclc5`` with :math:`a = 1.4\cdot\sin(53)\cdot\pi`, :math:`b = \pi`, 
:math:`c = 0.9\cdot 11\cdot\cos(53)\cdot\pi/6` and :math:`\alpha = 53^{\circ}`


Example structures
==================

MCLC\ :sub:`1`
--------------

**Default kpath**: :math:`\Gamma-Y-F-L-I\vert I_1-Z-F_1\vert Y-X_1\vert X-\Gamma-N\vert M-\Gamma`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mclc1_brillouin.png 
            :target: ../../../../../_images/mclc1_brillouin.png 
      - .. literalinclude:: mclc1_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mclc1_real.png 
            :target: ../../../../../_images/mclc1_real.png 
      - .. literalinclude:: mclc1_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mclc1_wigner-seitz.png 
            :target: ../../../../../_images/mclc1_wigner-seitz.png 
      - .. literalinclude:: mclc1_wigner-seitz.py
            :language: py

MCLC\ :sub:`2`
--------------

**Default kpath**: :math:`\Gamma-Y-F-L-I\vert I_1-Z-F_1\vert N-\Gamma-M`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mclc2_brillouin.png 
            :target: ../../../../../_images/mclc2_brillouin.png 
      - .. literalinclude:: mclc2_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mclc2_real.png 
            :target: ../../../../../_images/mclc2_real.png 
      - .. literalinclude:: mclc2_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mclc2_wigner-seitz.png 
            :target: ../../../../../_images/mclc2_wigner-seitz.png 
      - .. literalinclude:: mclc2_wigner-seitz.py
            :language: py

MCLC\ :sub:`3`
--------------

**Default kpath**: :math:`\Gamma-Y-F-H-Z-I-F_1\vert H_1-Y_1-X-\Gamma-N\vert M-\Gamma`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mclc3_brillouin.png 
            :target: ../../../../../_images/mclc3_brillouin.png 
      - .. literalinclude:: mclc3_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mclc3_real.png 
            :target: ../../../../../_images/mclc3_real.png 
      - .. literalinclude:: mclc3_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mclc3_wigner-seitz.png 
            :target: ../../../../../_images/mclc3_wigner-seitz.png 
      - .. literalinclude:: mclc3_wigner-seitz.py
            :language: py

MCLC\ :sub:`4`
--------------

**Default kpath**: :math:`\Gamma-Y-F-H-Z-I\vert H_1-Y_1-X-\Gamma-N\vert M-\Gamma`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mclc4_brillouin.png 
            :target: ../../../../../_images/mclc4_brillouin.png 
      - .. literalinclude:: mclc4_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mclc4_real.png 
            :target: ../../../../../_images/mclc4_real.png 
      - .. literalinclude:: mclc4_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mclc4_wigner-seitz.png 
            :target: ../../../../../_images/mclc4_wigner-seitz.png 
      - .. literalinclude:: mclc4_wigner-seitz.py
            :language: py

MCLC\ :sub:`5`
--------------

**Default kpath**: :math:`\Gamma-Y-F-L-I\vert I_1-Z-H-F_1\vert H_1-Y_1-X-\Gamma-N\vert M-\Gamma`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mclc5_brillouin.png 
            :target: ../../../../../_images/mclc5_brillouin.png 
      - .. literalinclude:: mclc5_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mclc5_real.png 
            :target: ../../../../../_images/mclc5_real.png 
      - .. literalinclude:: mclc5_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: mclc5_wigner-seitz.png 
            :target: ../../../../../_images/mclc5_wigner-seitz.png 
      - .. literalinclude:: mclc5_wigner-seitz.py
            :language: py


Ordering of lattice parameters
==============================
TODO

Edge cases
==========
TODO
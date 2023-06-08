.. _lattice-tri:

***************
Triclinic (TRI)
***************

**Pearson symbol**: aP

Triclinic lattice is described by the class :py:class:`.TRI`.

It is defined by six parameters: :math:`a`, :math:`b`, :math:`c` and
:math:`\alpha`, :math:`\beta`, :math:`\gamma`.
with primitive and conventional lattice:

.. math::

    \boldsymbol{a}_1 = (a, 0, 0)

    \boldsymbol{a}_2 = (b\cos(\gamma), b\sin(\gamma), 0)

    \boldsymbol{a}_3 = (c\cos(\beta), \frac{c(\cos(\alpha) - \cos(\beta)\cos(\gamma))}{\sin{\gamma}}, \frac{c}{\sin(\gamma)}\sqrt{\sin^2(\gamma) - \cos^2(\alpha) - \cos^2(\beta) + 2\cos(\alpha)\cos(\beta)\cos(\gamma)})


Variations
==========

There are four variation of triclinic lattice.

TRI\ :sub:`1a`
--------------

:math:`k_{\alpha} > 90^{\circ}, k_{\beta} > 90^{\circ}, k_{\gamma} > 90^{\circ}, k_{\gamma} = \min(k_{\alpha}, k_{\beta}, k_{\gamma})`

TRI\ :sub:`2a`
--------------

:math:`k_{\alpha} > 90^{\circ}, k_{\beta} > 90^{\circ}, k_{\gamma} = 90^{\circ}`

TRI\ :sub:`1b`
--------------

:math:`k_{\alpha} < 90^{\circ}, k_{\beta} < 90^{\circ}, k_{\gamma} < 90^{\circ}, k_{\gamma} = \max(k_{\alpha}, k_{\beta}, k_{\gamma})`

TRI\ :sub:`2b`
--------------

:math:`k_{\alpha} < 90^{\circ}, k_{\beta} < 90^{\circ}, k_{\gamma} = 90^{\circ}`

In definition of the examples we cheated and defined them through reciprocal lattice parameters.

Example structures
==================

TRI\ :sub:`1a`
--------------

**Default kpath**: :math:`\Gamma-Y-F-L-I\vert I_1-Z-F_1\vert Y-X_1\vert X-\Gamma-N\vert M-\Gamma`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: tri1a_brillouin.png 
            :target: ../../../../../_images/tri1a_brillouin.png 
      - .. literalinclude:: tri1a_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: tri1a_real.png 
            :target: ../../../../../_images/tri1a_real.png 
      - .. literalinclude:: tri1a_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: tri1a_wigner-seitz.png 
            :target: ../../../../../_images/tri1a_wigner-seitz.png 
      - .. literalinclude:: tri1a_wigner-seitz.py
            :language: py

TRI\ :sub:`2a`
--------------

**Default kpath**: :math:`\Gamma-Y-F-L-I\vert I_1-Z-F_1\vert Y-X_1\vert X-\Gamma-N\vert M-\Gamma`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: tri2a_brillouin.png 
            :target: ../../../../../_images/tri2a_brillouin.png 
      - .. literalinclude:: tri2a_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: tri2a_real.png 
            :target: ../../../../../_images/tri2a_real.png 
      - .. literalinclude:: tri2a_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: tri2a_wigner-seitz.png 
            :target: ../../../../../_images/tri2a_wigner-seitz.png 
      - .. literalinclude:: tri2a_wigner-seitz.py
            :language: py

TRI\ :sub:`1b`
--------------

**Default kpath**: :math:`\Gamma-Y-F-L-I\vert I_1-Z-F_1\vert Y-X_1\vert X-\Gamma-N\vert M-\Gamma`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: tri1b_brillouin.png 
            :target: ../../../../../_images/tri1b_brillouin.png 
      - .. literalinclude:: tri1b_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: tri1b_real.png 
            :target: ../../../../../_images/tri1b_real.png 
      - .. literalinclude:: tri1b_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: tri1b_wigner-seitz.png 
            :target: ../../../../../_images/tri1b_wigner-seitz.png 
      - .. literalinclude:: tri1b_wigner-seitz.py
            :language: py

TRI\ :sub:`2b`
--------------

**Default kpath**: :math:`\Gamma-Y-F-L-I\vert I_1-Z-F_1\vert Y-X_1\vert X-\Gamma-N\vert M-\Gamma`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: tri2b_brillouin.png 
            :target: ../../../../../_images/tri2b_brillouin.png 
      - .. literalinclude:: tri2b_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: tri2b_real.png 
            :target: ../../../../../_images/tri2b_real.png 
      - .. literalinclude:: tri2b_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: tri2b_wigner-seitz.png 
            :target: ../../../../../_images/tri2b_wigner-seitz.png 
      - .. literalinclude:: tri2b_wigner-seitz.py
            :language: py

Ordering of parameters
======================
TODO
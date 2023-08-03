.. _guide_rhl:

******************
Rhombohedral (RHL)
******************

**Pearson symbol**: hR

Rhombohedral lattice is described by the class :py:class:`.RHL`.

It is defined by two parameter: :math:`a` and :math:`\alpha` 
with primitive and conventional lattice:

.. math::


    \boldsymbol{a}_1 = (a\cos(\alpha / 2), -a\sin(\alpha/2), 0)

    \boldsymbol{a}_2 = (a\cos(\alpha / 2), a\sin(\alpha/2), 0)

    \boldsymbol{a}_3 = \left(\frac{\cos(\alpha)}{\cos(\alpha/2)}, 0, a\sqrt{1 - \frac{\cos^2(\alpha)}{\cos^2(\alpha/2)}}\right)


Variations
==========

There are two variations for rhombohedral lattice.

RHL\ :sub:`1`
-------------

:math:`\alpha < 90^{\circ}`.

Predefined example: ``rhl1`` with :math:`a = \pi` and :math:`\alpha = 70` 

RHL\ :sub:`2`
-------------

:math:`\alpha > 90^{\circ}`.

Predefined example: ``rhl2`` with :math:`a = \pi` and :math:`\alpha = 110` 

Example structure
=================

RHL\ :sub:`1`
-------------

**Default kpath**: :math:`\Gamma-L-B_1\vert B-Z-\Gamma-X\vert Q-F-P_1-Z\vert L-P`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: rhl1_brillouin.png 
            :target: ../../../../../_images/rhl1_brillouin.png 
      - .. literalinclude:: rhl1_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: rhl1_real.png 
            :target: ../../../../../_images/rhl1_real.png 
      - .. literalinclude:: rhl1_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: rhl1_wigner-seitz.png 
            :target: ../../../../../_images/rhl1_wigner-seitz.png 
      - .. literalinclude:: rhl1_wigner-seitz.py
            :language: py

RHL\ :sub:`2`
-------------

**Default kpath**: :math:`\Gamma-P-Z-Q-\Gamma-F-P_1-Q_1-L-Z`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: rhl2_brillouin.png 
            :target: ../../../../../_images/rhl2_brillouin.png 
      - .. literalinclude:: rhl2_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: rhl2_real.png 
            :target: ../../../../../_images/rhl2_real.png 
      - .. literalinclude:: rhl2_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: rhl2_wigner-seitz.png 
            :target: ../../../../../_images/rhl2_wigner-seitz.png 
      - .. literalinclude:: rhl2_wigner-seitz.py
            :language: py


Edge cases
==========
In rhombohedral lattice :math:`a = b = c` and :math:`\alpha = \beta = \gamma`, 
thus three edge cases exist:

If :math:`\alpha = 60^{\circ}`, then the lattice is :ref:`guide_fcc`

If :math:`\alpha \approx 109.47122063^{\circ}` (:math:`\cos(\alpha) = -1/3`), 
then the lattice is :ref:`guide_bcc`.

If :math:`\alpha = 90^{\circ}`, then the lattice is :ref:`guide_cub`.
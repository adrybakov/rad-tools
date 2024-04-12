.. _guide_rhl:

******************
Rhombohedral (RHL)
******************

**Pearson symbol**: hR

**Constructor**:  :py:func:`.RHL`

It is defined by two parameter: :math:`a` and :math:`\alpha`
with primitive and conventional lattice:

.. math::

    \begin{matrix}
    \boldsymbol{a}_1 &=& (a\cos(\alpha / 2), &-a\sin(\alpha/2), &0)\\
    \boldsymbol{a}_2 &=& (a\cos(\alpha / 2), &a\sin(\alpha/2), &0)\\
    \boldsymbol{a}_3 &=& \left(\dfrac{\cos\alpha}{\cos(\alpha/2)}\right.,
    &0, &\left.a\sqrt{1 - \dfrac{\cos^2\alpha}{\cos^2(\alpha/2)}}\right)
    \end{matrix}

Cell standardization
====================

No standardization is performed.

K-path
======

RHL\ :sub:`1`
-------------

:math:`\mathrm{\Gamma-L-B_1\vert B-Z-\Gamma-X\vert Q-F-P_1-Z\vert L-P}`

.. math::

    \begin{matrix}
    \eta = \dfrac{1 + 4\cos\alpha}{2 + 4\cos\alpha} &
    \nu = \dfrac{3-2\eta}{4}
    \end{matrix}

=========================  ==============================  ==============================  ==============================
Point                      :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=========================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`    :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{B}`         :math:`\eta`                    :math:`1/2`                     :math:`1 - \eta`
:math:`\mathrm{B_1}`       :math:`1/2`                     :math:`1-\eta`                  :math:`\eta - 1`
:math:`\mathrm{F}`         :math:`1/2`                     :math:`1/2`                     :math:`0`
:math:`\mathrm{L}`         :math:`1/2`                     :math:`0`                       :math:`0`
:math:`\mathrm{L_1}`       :math:`0`                       :math:`0`                       :math:`-1/2`
:math:`\mathrm{P}`         :math:`\eta`                    :math:`\nu`                     :math:`\nu`
:math:`\mathrm{P_1}`       :math:`1-\nu`                   :math:`1-\nu`                   :math:`1-\eta`
:math:`\mathrm{P_2}`       :math:`\nu`                     :math:`\nu`                     :math:`\eta - 1`
:math:`\mathrm{Q}`         :math:`1-\nu`                   :math:`\nu`                     :math:`0`
:math:`\mathrm{X}`         :math:`\nu`                     :math:`0`                       :math:`-\nu`
:math:`\mathrm{Z}`         :math:`1/2`                     :math:`1/2`                     :math:`1/2`
=========================  ==============================  ==============================  ==============================

RHL\ :sub:`2`
-------------

:math:`\mathrm{\Gamma-P-Z-Q-\Gamma-F-P_1-Q_1-L-Z}`

.. math::

    \begin{matrix}
    \eta = \dfrac{1}{2\tan^2(\alpha/2)} &
    \nu = \dfrac{3-2\eta}{4}
    \end{matrix}

=========================  ==============================  ==============================  ==============================
Point                      :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=========================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`    :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{F}`         :math:`1/2`                     :math:`-1/2`                    :math:`0`
:math:`\mathrm{L}`         :math:`1/2`                     :math:`0`                       :math:`0`
:math:`\mathrm{P}`         :math:`1-\nu`                   :math:`-\nu`                    :math:`1-\nu`
:math:`\mathrm{P_1}`       :math:`\nu`                     :math:`\nu-1`                   :math:`\nu-1`
:math:`\mathrm{Q}`         :math:`\eta`                    :math:`\eta`                    :math:`\eta`
:math:`\mathrm{Q_1}`       :math:`1-\eta`                  :math:`-\eta`                   :math:`-\eta`
:math:`\mathrm{Z}`         :math:`1/2`                     :math:`-1/2`                    :math:`1/2`
=========================  ==============================  ==============================  ==============================

Variations
==========

There are two variations for rhombohedral lattice.

RHL\ :sub:`1`
-------------

:math:`\alpha < 90^{\circ}`.

Predefined example: ``rhl1`` with :math:`a = \pi` and :math:`\alpha = 70^{\circ}`

RHL\ :sub:`2`
-------------

:math:`\alpha > 90^{\circ}`.

Predefined example: ``rhl2`` with :math:`a = \pi` and :math:`\alpha = 110^{\circ}`

Examples
========

RHL\ :sub:`1`
-------------

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: rhl1_brillouin.py
    :language: py

.. raw:: html
    :file: rhl1_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: rhl1_real.py
    :language: py

.. raw:: html
    :file: rhl1_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: rhl1_wigner-seitz.py
    :language: py

.. raw:: html
    :file: rhl1_wigner-seitz.html

RHL\ :sub:`2`
-------------

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: rhl2_brillouin.py
    :language: py

.. raw:: html
    :file: rhl2_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: rhl2_real.py
    :language: py

.. raw:: html
    :file: rhl2_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: rhl2_wigner-seitz.py
    :language: py

.. raw:: html
    :file: rhl2_wigner-seitz.html

Edge cases
==========
In rhombohedral lattice :math:`a = b = c` and :math:`\alpha = \beta = \gamma`,
thus three edge cases exist:

If :math:`\alpha = 60^{\circ}`, then the lattice is :ref:`guide_fcc`

If :math:`\alpha \approx 109.47122063^{\circ}` (:math:`\cos(\alpha) = -1/3`),
then the lattice is :ref:`guide_bcc`.

If :math:`\alpha = 90^{\circ}`, then the lattice is :ref:`guide_cub`.

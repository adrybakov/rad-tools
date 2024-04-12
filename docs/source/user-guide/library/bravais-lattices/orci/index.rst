.. _guide_orci:

********************************
Body-centred orthorhombic (ORCI)
********************************

**Pearson symbol**: oI

**Constructor**:  :py:func:`.ORCI`

It is defined by three parameter: :math:`a`, :math:`b` and :math:`c`
with conventional lattice:

.. math::

    \begin{matrix}
    \boldsymbol{a}_1 &=& (a, &0, &0)\\
    \boldsymbol{a}_2 &=& (0, &b, &0)\\
    \boldsymbol{a}_3 &=& (0, &0, &c)
    \end{matrix}

And primitive lattice:

.. math::

    \begin{matrix}
    \boldsymbol{a}_1 &=& (-a/2, &b/2, &c/2)\\
    \boldsymbol{a}_2 &=& (a/2, &-b/2, &c/2)\\
    \boldsymbol{a}_3 &=& (a/2, &b/2, &-c/2)
    \end{matrix}

Order of parameters: :math:`a < b < c`

Cell standardization
====================

Lengths of the lattice vectors of the conventional cell have to satisfy
:math:`\vert\boldsymbol{a}_1\vert < \vert\boldsymbol{a}_2\vert < \vert\boldsymbol{a}_3\vert`.
All vectors of the primitive cell have the same length, therefore we use
conventional lattice vectors for the standardization.

If this condition is not satisfied, then the lattice is transformed to the standard form:

First we order first two vectors by length:

* If :math:`\vert\boldsymbol{a}_1\vert > \vert\boldsymbol{a}_2\vert`
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow (\boldsymbol{a}_2, \boldsymbol{a}_1, -\boldsymbol{a}_3)

Then we find a correct place for the third vector:

* If :math:`\vert\boldsymbol{a}_1\vert > \vert\boldsymbol{a}_3\vert`
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow (\boldsymbol{a}_3, \boldsymbol{a}_1, \boldsymbol{a}_2)

* If :math:`\vert\boldsymbol{a}_2\vert > \vert\boldsymbol{a}_3\vert`
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow (\boldsymbol{a}_1, -\boldsymbol{a}_3, \boldsymbol{a}_2)

.. note::

    The third lattice vector is multiplied by :math:`-1` in some cases to
    preserve the handedness of the cell.

K-path
======

:math:`\mathrm{\Gamma-X-L-T-W-R-X_1-Z-\Gamma-Y-S-W\vert L_1-Y\vert Y_1-Z}`

.. math::

    \begin{matrix}
    \zeta = \dfrac{1 + a^2/c^2}{4} &
    \eta = \dfrac{1 + b^2/c^2}{4} &
    \delta = \dfrac{b^2 - a^2}{4c^2} &
    \mu = \dfrac{a^2 + b^2}{4c^2}
    \end{matrix}

=========================  ==============================  ==============================  ==============================
Point                      :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=========================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`    :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{L}`         :math:`-\mu`                    :math:`\mu`                     :math:`1/2 - \delta`
:math:`\mathrm{L_1}`       :math:`\mu`                     :math:`-\mu`                    :math:`1/2 + \delta`
:math:`\mathrm{L_2}`       :math:`1/2-\delta`              :math:`1/2+\delta`              :math:`-\mu`
:math:`\mathrm{R}`         :math:`0`                       :math:`1/2`                     :math:`0`
:math:`\mathrm{S}`         :math:`1/2`                     :math:`0`                       :math:`0`
:math:`\mathrm{T}`         :math:`0`                       :math:`0`                       :math:`1/2`
:math:`\mathrm{W}`         :math:`1/4`                     :math:`1/4`                     :math:`1/4`
:math:`\mathrm{X}`         :math:`-\zeta`                  :math:`\zeta`                   :math:`\zeta`
:math:`\mathrm{X_1}`       :math:`\zeta`                   :math:`1-\zeta`                 :math:`-\zeta`
:math:`\mathrm{Y}`         :math:`\eta`                    :math:`-\eta`                   :math:`\eta`
:math:`\mathrm{Y_1}`       :math:`1-\eta`                  :math:`\eta`                    :math:`-\eta`
:math:`\mathrm{Z}`         :math:`1/2`                     :math:`1/2`                     :math:`-1/2`
=========================  ==============================  ==============================  ==============================


Variations
==========

There are no variations for body-centered orthorombic.
One example is predefined: ``orci`` with
:math:`a = \pi`, :math:`b  = 1.3\pi` and :math:`c = 1.7\pi`.

Examples
========

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: orci_brillouin.py
    :language: py

.. raw:: html
    :file: orci_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: orci_real.py
    :language: py

.. raw:: html
    :file: orci_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: orci_wigner-seitz.py
    :language: py

.. raw:: html
    :file: orci_wigner-seitz.html

Edge cases
==========
If :math:`a = b \ne c` or :math:`a = c \ne b` or :math:`b = c \ne a`,
then the lattice is :ref:`guide_bct`.

If :math:`a = b = c`, then the lattice is :ref:`guide_bcc`.

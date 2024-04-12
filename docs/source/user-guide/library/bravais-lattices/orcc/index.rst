.. _guide_orcc:

********************************
Base-centred orthorhombic (ORCC)
********************************

**Pearson symbol**: oS

**Constructor**:  :py:func:`.ORCC`

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
        \boldsymbol{a}_1 &=& (a/2, &-b/2, &0)\\
        \boldsymbol{a}_2 &=& (a/2, &b/2, &0)\\
        \boldsymbol{a}_3 &=& (0, &0, &c)
    \end{matrix}

Order of parameters: :math:`a < b`

Cell standardization
====================

Lengths of the lattice vectors of the conventional cell have to satisfy
:math:`\vert\boldsymbol{a}_1\vert < \vert\boldsymbol{a}_2\vert`.
Length of third vector of the primitive cell has to be different from
the lengths of the first two vectors of the primitive cell.

If these conditions are not satisfied, then the lattice is transformed to the standard form:

Before the standardization we check if all angles are equal to :math:`90^{\circ}`.
If yes, then the standardization for :ref:`guide_tet` is called.

First we find the third lattice vector of the primitive cell.
For this step we use lattice vectors of the primitive cell:

* If :math:`\beta \ne 90^{\circ}`
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_3, \boldsymbol{a}_1, \boldsymbol{a}_2)

* If :math:`\alpha \ne 90^{\circ}`
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_2, \boldsymbol{a}_3, \boldsymbol{a}_1)

Then we order first two lattice vector.
For this step we use lattice vectors of the conventional cell:

* If :math:`\vert\boldsymbol{a}_1\vert > \vert\boldsymbol{a}_2\vert`
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (-\boldsymbol{a}_2, \boldsymbol{a}_1, \boldsymbol{a}_3)

.. note::

    The second lattice vector is multiplied by :math:`-1` in some cases to
    preserve the handedness of the cell.

K-path
======

:math:`\mathrm{\Gamma-X-S-R-A-Z-\Gamma-Y-X_1-A_1-T-Y\vert Z-T}`

.. math::

    \zeta = \dfrac{1 + a^2/b^2}{4}

=========================  ==============================  ==============================  ==============================
Point                      :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=========================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`    :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{A}`         :math:`\zeta`                   :math:`\zeta`                   :math:`1/2`
:math:`\mathrm{A_1}`       :math:`-\zeta`                  :math:`1-\zeta`                 :math:`1/2`
:math:`\mathrm{R}`         :math:`0`                       :math:`1/2`                     :math:`1/2`
:math:`\mathrm{S}`         :math:`0`                       :math:`1/2`                     :math:`0`
:math:`\mathrm{T}`         :math:`-1/2`                    :math:`1/2`                     :math:`1/2`
:math:`\mathrm{X}`         :math:`\zeta`                   :math:`\zeta`                   :math:`0`
:math:`\mathrm{X_1}`       :math:`-\zeta`                  :math:`1-\zeta`                 :math:`0`
:math:`\mathrm{Y}`         :math:`-1/2`                    :math:`1/2`                     :math:`0`
:math:`\mathrm{Z}`         :math:`0`                       :math:`0`                       :math:`1/2`
=========================  ==============================  ==============================  ==============================


Variations
==========

There are no variations for base-centered orthorombic.
One example is predefined: ``orcc`` with
:math:`a = \pi`, :math:`b  = 1.3\pi` and :math:`c = 1.7\pi`.

Examples
========

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: orcc_brillouin.py
    :language: py

.. raw:: html
    :file: orcc_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: orcc_real.py
    :language: py

.. raw:: html
    :file: orcc_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: orcc_wigner-seitz.py
    :language: py

.. raw:: html
    :file: orcc_wigner-seitz.html


Edge cases
==========
If :math:`a = b`, then the lattice is :ref:`guide_tet`.

If :math:`a = b = \sqrt{2} c`, then the lattice is :ref:`guide_cub`.

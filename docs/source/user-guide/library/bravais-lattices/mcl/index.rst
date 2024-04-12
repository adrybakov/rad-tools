.. _guide_mcl:

****************
Monoclinic (MCL)
****************

**Pearson symbol**: mP

**Constructor**:  :py:func:`.MCL`

It is defined by four parameter: :math:`a`, :math:`b`, :math:`c` and :math:`\alpha`
with primitive and conventional lattice:

.. math::

    \begin{matrix}
    \boldsymbol{a}_1 &=& (a, &0, &0)\\
    \boldsymbol{a}_2 &=& (0, &b, &0)\\
    \boldsymbol{a}_3 &=& (0, &c\cos\alpha, &c\sin\alpha)
    \end{matrix}


Order of parameters: :math:`b \le c`, :math:`\alpha < 90^{\circ}`.

Cell standardization
====================

Angles between first and third lattice vectors (:math:`\beta`) and between
first and second lattice vectors (:math:`\gamma`) have to be equal to :math:`90^{\circ}`.

Angle between second and third lattice vectors (:math:`\alpha`) has to be less then :math:`90^{\circ}`.

Length of the second lattice vector has to be less or equal to the length of the third lattice vector.

If these conditions are not satisfied, then the lattice is transformed to the standard form:

First we ensure the :math:`\beta = 90^{\circ}` and :math:`\gamma = 90^{\circ}`:

* If :math:`\beta \ne 90^{\circ}`:
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_2, \boldsymbol{a}_3, \boldsymbol{a}_1)

* If :math:`\gamma \ne 90^{\circ}`:
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_3, \boldsymbol{a}_1, \boldsymbol{a}_2)

Then we ensure the :math:`\alpha < 90^{\circ}`:

* If :math:`\alpha > 90^{\circ}`:
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_1, \boldsymbol{a}_3, -\boldsymbol{a}_2)

Finally, we ensure the :math:`b \le c`:

* If :math:`b > c`:
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (-\boldsymbol{a}_1, \boldsymbol{a}_3, \boldsymbol{a}_2)

.. note::

    First and second lattice vectors are multiplied by :math:`-1` in some cases to
    preserve the handedness of the cell.

K-path
======

:math:`\mathrm{\Gamma-Y-H-C-E-M_1-A-X-H_1\vert M-D-Z\vert Y-D}`

.. math::

    \begin{matrix}
    \eta = \dfrac{1 - b\cos\alpha / c}{2\sin^2\alpha} &
    \nu = \dfrac{1}{2} - \dfrac{\eta c\cos\alpha}{b}
    \end{matrix}

=========================  ==============================  ==============================  ==============================
Point                      :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=========================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`    :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{A}`         :math:`1/2`                     :math:`1/2`                     :math:`0`
:math:`\mathrm{C}`         :math:`0`                       :math:`1/2`                     :math:`1/2`
:math:`\mathrm{D}`         :math:`1/2`                     :math:`0`                       :math:`1/2`
:math:`\mathrm{D_1}`       :math:`1/2`                     :math:`0`                       :math:`-1/2`
:math:`\mathrm{E}`         :math:`1/2`                     :math:`1/2`                     :math:`1/2`
:math:`\mathrm{H}`         :math:`0`                       :math:`\eta`                    :math:`1-\nu`
:math:`\mathrm{H_1}`       :math:`0`                       :math:`1-\eta`                  :math:`\nu`
:math:`\mathrm{H_2}`       :math:`0`                       :math:`\eta`                    :math:`-\nu`
:math:`\mathrm{M}`         :math:`1/2`                     :math:`\eta`                    :math:`1-\nu`
:math:`\mathrm{M_1}`       :math:`1/2`                     :math:`1-\eta`                  :math:`\nu`
:math:`\mathrm{M_2}`       :math:`1/2`                     :math:`\eta`                    :math:`-\nu`
:math:`\mathrm{X}`         :math:`0`                       :math:`1/2`                     :math:`0`
:math:`\mathrm{Y}`         :math:`0`                       :math:`0`                       :math:`1/2`
:math:`\mathrm{Y_1}`       :math:`0`                       :math:`0`                       :math:`-1/2`
:math:`\mathrm{Z}`         :math:`1/2`                     :math:`0`                       :math:`0`
=========================  ==============================  ==============================  ==============================

Variations
==========

There are two variations for monoclinic lattice.
One example is predefined: ``mcl`` with
MCL(pi, 1.3 * pi, 1.6 * pi, alpha=75)
:math:`a = \pi`, :math:`b = 1.3 \pi` :math:`c = 1.6 \pi` and :math:`\alpha = 75^{\circ}`.

Examples
========
Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: mcl_brillouin.py
    :language: py

.. raw:: html
    :file: mcl_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: mcl_real.py
    :language: py

.. raw:: html
    :file: mcl_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: mcl_wigner-seitz.py
    :language: py

.. raw:: html
    :file: mcl_wigner-seitz.html

Edge cases
==========

If (:math:`\alpha = 60^{\circ}` or :math:`\alpha = 120^{\circ}`) and :math:`b = c`,
then the lattice is :ref:`guide_hex`.

If (:math:`\alpha = 30^{\circ}` or :math:`\alpha = 150^{\circ}`
or :math:`\alpha = 45^{\circ}` or :math:`\alpha = 145^{\circ}`) and :math:`b = c`,
then the lattice is :ref:`guide_orcc`.

If (:math:`\alpha = 60^{\circ}` or :math:`\alpha = 120^{\circ}`) and :math:`a \ne b = c/2`,
then the lattice is :ref:`guide_orc`.

If :math:`a \ne b \ne c` and :math:`\alpha = 90^{\circ}`, then the lattice is :ref:`guide_orc`.

If (:math:`\alpha = 60^{\circ}` or :math:`\alpha = 120^{\circ}`) and :math:`a = b = c/2`,
then the lattice is :ref:`guide_tet`.

If (:math:`a = b \ne c` or :math:`a = c \ne b` or :math:`b = c \ne a`) and :math:`\alpha = 90^{\circ}`, then the lattice is :ref:`guide_tet`.

If :math:`a = b = c` and :math:`\alpha = 90^{\circ}`, then the lattice is :ref:`guide_cub`.

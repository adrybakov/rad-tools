.. _guide_orcf:

********************************
Face-centred orthorhombic (ORCF)
********************************

**Pearson symbol**: oF

**Constructor**:  :py:func:`.ORCF`

It is defined by three parameters: :math:`a`, :math:`b` and :math:`c`
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
    \boldsymbol{a}_1 &=& (0, &b/2, &c/2)\\
    \boldsymbol{a}_2 &=& (a/2, &0, &c/2)\\
    \boldsymbol{a}_3 &=& (a/2, &b/2, &0)
    \end{matrix}

Cell standardization
====================

Lengths of the lattice vectors of the conventional cell have to satisfy
:math:`\vert\boldsymbol{a}_1\vert < \vert\boldsymbol{a}_2\vert < \vert\boldsymbol{a}_3\vert`.
That means that lattice vectors of the primitive cell have to satisfy:
:math:`\vert\boldsymbol{a}_1\vert > \vert\boldsymbol{a}_2\vert > \vert\boldsymbol{a}_3\vert`.


If this condition is not satisfied, then the lattice is transformed to the standard form.
We use the primitive lattice vectors for the standardization:

First we order first two vectors by length:

* If :math:`\vert\boldsymbol{a}_1\vert < \vert\boldsymbol{a}_2\vert`
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_2, \boldsymbol{a}_1, -\boldsymbol{a}_3)

Then we find a correct place for the third vector:

* If :math:`\vert\boldsymbol{a}_1\vert < \vert\boldsymbol{a}_3\vert`
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_3, \boldsymbol{a}_1, \boldsymbol{a}_2)

* If :math:`\vert\boldsymbol{a}_2\vert < \vert\boldsymbol{a}_3\vert`
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_1, -\boldsymbol{a}_3, \boldsymbol{a}_2)

.. note::

    The third lattice vector is multiplied by :math:`-1` in some cases to
    preserve the handedness of the cell.

K-path
======

ORCF\ :sub:`1`
--------------

:math:`\mathrm{\Gamma-Y-T-Z-\Gamma-X-A_1-Y\vert T-X_1\vert X-A-Z\vert L-\Gamma}`

.. math::

    \begin{matrix}
    \zeta = \dfrac{1 + a^2/b^2 - a^2/c^2}{4} &
    \eta = \dfrac{1 + a^2/b^2 + a^2/c^2}{4}
    \end{matrix}

=========================  ==============================  ==============================  ==============================
Point                      :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=========================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`    :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{A}`         :math:`1/2`                     :math:`1/2 + \zeta`             :math:`\zeta`
:math:`\mathrm{A_1}`       :math:`1/2`                     :math:`1/2 - \zeta`             :math:`1-\zeta`
:math:`\mathrm{L}`         :math:`1/2`                     :math:`1/2`                     :math:`1/2`
:math:`\mathrm{T}`         :math:`1`                       :math:`1/2`                     :math:`1/2`
:math:`\mathrm{X}`         :math:`0`                       :math:`\eta`                    :math:`\eta`
:math:`\mathrm{X_1}`       :math:`1`                       :math:`1-\eta`                  :math:`1-\eta`
:math:`\mathrm{Y}`         :math:`1/2`                     :math:`0`                       :math:`1/2`
:math:`\mathrm{Z}`         :math:`1/2`                     :math:`1/2`                     :math:`0`
=========================  ==============================  ==============================  ==============================

ORCF\ :sub:`2`
--------------

:math:`\mathrm{\Gamma-Y-C-D-X-\Gamma-Z-D_1-H-C\vert C_1-Z\vert X-H_1\vert H-Y\vert L-\Gamma}`

.. math::

    \begin{matrix}
    \eta = \dfrac{1 + a^2/b^2 - a^2/c^2}{4} &
    \delta = \dfrac{1 + b^2/a^2 - b^2/c^2}{4} &
    \phi = \dfrac{1 + c^2/b^2 - c^2/a^2}{4}
    \end{matrix}

=========================  ==============================  ==============================  ==============================
Point                      :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=========================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`    :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{C}`         :math:`1/2`                     :math:`1/2 - \eta`              :math:`1 - \eta`
:math:`\mathrm{C_1}`       :math:`1/2`                     :math:`1/2 + \eta`              :math:`\eta`
:math:`\mathrm{D}`         :math:`1/2-\delta`              :math:`1/2`                     :math:`1 - \delta`
:math:`\mathrm{D_1}`       :math:`1/2+\delta`              :math:`1/2`                     :math:`\delta`
:math:`\mathrm{L}`         :math:`1/2`                     :math:`1/2`                     :math:`1/2`
:math:`\mathrm{H}`         :math:`1 - \phi`                :math:`1/2 - \phi`              :math:`1/2`
:math:`\mathrm{H_1}`       :math:`\phi`                    :math:`1/2 + \phi`              :math:`1/2`
:math:`\mathrm{X}`         :math:`0`                       :math:`1/2`                     :math:`1/2`
:math:`\mathrm{Y}`         :math:`1/2`                     :math:`0`                       :math:`1/2`
:math:`\mathrm{Z}`         :math:`1/2`                     :math:`1/2`                     :math:`0`
=========================  ==============================  ==============================  ==============================

ORCF\ :sub:`3`
--------------

:math:`\mathrm{\Gamma-Y-T-Z-\Gamma-X-A_1-Y\vert X-A-Z\vert L-\Gamma}`

.. math::

    \begin{matrix}
    \zeta = \dfrac{1 + a^2/b^2 - a^2/c^2}{4} &
    \eta = \dfrac{1 + a^2/b^2 + a^2/c^2}{4}
    \end{matrix}

=========================  ==============================  ==============================  ==============================
Point                      :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=========================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`    :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{A}`         :math:`1/2`                     :math:`1/2 + \zeta`             :math:`\zeta`
:math:`\mathrm{A_1}`       :math:`1/2`                     :math:`1/2 - \zeta`             :math:`1-\zeta`
:math:`\mathrm{L}`         :math:`1/2`                     :math:`1/2`                     :math:`1/2`
:math:`\mathrm{T}`         :math:`1`                       :math:`1/2`                     :math:`1/2`
:math:`\mathrm{X}`         :math:`0`                       :math:`\eta`                    :math:`\eta`
:math:`\mathrm{X_1}`       :math:`1`                       :math:`1-\eta`                  :math:`1-\eta`
:math:`\mathrm{Y}`         :math:`1/2`                     :math:`0`                       :math:`1/2`
:math:`\mathrm{Z}`         :math:`1/2`                     :math:`1/2`                     :math:`0`
=========================  ==============================  ==============================  ==============================

Variations
==========

There are three  variations of face-centered orthorombic lattice.

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

Examples
========

ORCF\ :sub:`1`
--------------

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: orcf1_brillouin.py
    :language: py

.. raw:: html
    :file: orcf1_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: orcf1_real.py
    :language: py

.. raw:: html
    :file: orcf1_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: orcf1_wigner-seitz.py
    :language: py

.. raw:: html
    :file: orcf1_wigner-seitz.html

ORCF\ :sub:`2`
--------------

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: orcf2_brillouin.py
    :language: py

.. raw:: html
    :file: orcf2_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: orcf2_real.py
    :language: py

.. raw:: html
    :file: orcf2_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: orcf2_wigner-seitz.py
    :language: py

.. raw:: html
    :file: orcf2_wigner-seitz.html

ORCF\ :sub:`3`
--------------

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: orcf3_brillouin.py
    :language: py

.. raw:: html
    :file: orcf3_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: orcf3_real.py
    :language: py

.. raw:: html
    :file: orcf3_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: orcf3_wigner-seitz.py
    :language: py

.. raw:: html
    :file: orcf3_wigner-seitz.html

Edge cases
==========
If :math:`a = b \ne c` or :math:`a = c \ne b` or :math:`b = c \ne a`,
then the lattice is :ref:`guide_bct`.

If :math:`a = b = c`, then the lattice is :ref:`guide_fcc`.

.. _guide_orc:

******************
Orthorhombic (ORC)
******************

**Pearson symbol**: oP

**Constructor**:  :py:func:`.ORC`


It is defined by three parameter: :math:`a`, :math:`b` and :math:`c`
with primitive and conventional lattice:

.. math::

    \begin{matrix}
    \boldsymbol{a}_1 &=& (a, &0, &0)\\
    \boldsymbol{a}_2 &=& (0, &b, &0)\\
    \boldsymbol{a}_3 &=& (0, &0, &c)
    \end{matrix}

Order of parameters: :math:`a < b < c`

Cell standardization
====================

Lengths of the lattice vectors have to satisfy
:math:`\vert\boldsymbol{a}_1\vert < \vert\boldsymbol{a}_2\vert < \vert\boldsymbol{a}_3\vert`.

If this condition is not satisfied, then the lattice is transformed to the standard form:

First we order first two vectors by length:

* If :math:`\vert\boldsymbol{a}_1\vert > \vert\boldsymbol{a}_2\vert`
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_2, \boldsymbol{a}_1, -\boldsymbol{a}_3)

Then we find a correct place for the third vector:

* If :math:`\vert\boldsymbol{a}_1\vert > \vert\boldsymbol{a}_3\vert`
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_3, \boldsymbol{a}_1, \boldsymbol{a}_2)

* If :math:`\vert\boldsymbol{a}_2\vert > \vert\boldsymbol{a}_3\vert`
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_1, -\boldsymbol{a}_3, \boldsymbol{a}_2)

.. note::

    The third lattice vector is multiplied by :math:`-1` in some cases to
    preserve the handedness of the cell.

K-path
======

:math:`\mathrm{\Gamma-X-S-Y-\Gamma-Z-U-R-T-Z\vert Y-T\vert U-X\vert S-R}`

=======================  ==============================  ==============================  ==============================
Point                    :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=======================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`  :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{R}`       :math:`1/2`                     :math:`1/2`                     :math:`1/2`
:math:`\mathrm{S}`       :math:`1/2`                     :math:`1/2`                     :math:`0`
:math:`\mathrm{T}`       :math:`0`                       :math:`1/2`                     :math:`1/2`
:math:`\mathrm{U}`       :math:`1/2`                     :math:`0`                       :math:`1/2`
:math:`\mathrm{X}`       :math:`1/2`                     :math:`0`                       :math:`0`
:math:`\mathrm{Y}`       :math:`0`                       :math:`1/2`                     :math:`0`
:math:`\mathrm{Z}`       :math:`0`                       :math:`0`                       :math:`1/2`
=======================  ==============================  ==============================  ==============================

Variations
==========

There are no variations for orthorhombic lattice.
One example is predefined: ``orc`` with
:math:`a = \pi`, :math:`b  = 1.5\pi` and :math:`c = 2\pi`.

Examples
========

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: orc_brillouin.py
    :language: py

.. raw:: html
    :file: orc_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: orc_real.py
    :language: py

.. raw:: html
    :file: orc_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: orc_wigner-seitz.py
    :language: py

.. raw:: html
    :file: orc_wigner-seitz.html

Edge cases
==========
If :math:`a = b \ne c` or :math:`a = c \ne b` or :math:`b = c \ne a`,
then the lattice is :ref:`guide_tet`.

If :math:`a = b = c`, then the lattice is :ref:`guide_cub`.

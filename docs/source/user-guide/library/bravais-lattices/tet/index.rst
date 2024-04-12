.. _guide_tet:

****************
Tetragonal (TET)
****************

**Pearson symbol**: tP

**Constructor**:  :py:func:`.TET`

It is defined by two parameters: :math:`a` and :math:`c`
with primitive and conventional lattice:

.. math::

    \begin{matrix}
    \boldsymbol{a}_1 &=& (a, &0, &0)\\
    \boldsymbol{a}_2 &=& (0, &a, &0)\\
    \boldsymbol{a}_3 &=& (0, &0, &c)
    \end{matrix}

Cell standardization
====================

Length of third lattice vector has to be different from the first two.
If this condition is not satisfied,
then the lattice is transformed to the standard form:

* If :math:`\vert\boldsymbol{a}_1\vert = \vert\boldsymbol{a}_3\vert \ne \vert\boldsymbol{a}_2\vert`
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_3, \boldsymbol{a}_1, \boldsymbol{a}_2)

* If :math:`\vert\boldsymbol{a}_2\vert = \vert\boldsymbol{a}_3\vert \ne \vert\boldsymbol{a}_1\vert`
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_2, \boldsymbol{a}_3, \boldsymbol{a}_1)


K-path
======

:math:`\mathrm{\Gamma-X-M-\Gamma-Z-R-A-Z\vert X-R\vert M-A}`

=======================  ==============================  ==============================  ==============================
Point                    :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=======================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`  :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{A}`       :math:`1/2`                     :math:`1/2`                     :math:`1/2`
:math:`\mathrm{M}`       :math:`1/2`                     :math:`1/2`                     :math:`0`
:math:`\mathrm{R}`       :math:`0`                       :math:`1/2`                     :math:`1/2`
:math:`\mathrm{X}`       :math:`0`                       :math:`1/2`                     :math:`0`
:math:`\mathrm{Z}`       :math:`0`                       :math:`0`                       :math:`1/2`
=======================  ==============================  ==============================  ==============================

Variations
==========

There are no variations for tetragonal lattice.
One example is predefined: ``tet`` with :math:`a = \pi` and :math:`c = 1.5\pi`

Examples
========

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: tet_brillouin.py
    :language: py

.. raw:: html
    :file: tet_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: tet_real.py
    :language: py

.. raw:: html
    :file: tet_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: tet_wigner-seitz.py
    :language: py

.. raw:: html
    :file: tet_wigner-seitz.html


Edge cases
==========

If :math:`a = c`, then the lattice is :ref:`guide_cub`.

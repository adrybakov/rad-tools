.. _guide_cub:

***********
Cubic (CUB)
***********

**Pearson symbol**: cP

**Constructor**:  :py:func:`.CUB`

It is defined by one parameter: :math:`a` with primitive and conventional lattice:

.. math::

    \begin{matrix}
    \boldsymbol{a}_1 &=& (a, &0, &0)\\
    \boldsymbol{a}_2 &=& (0, &a, &0)\\
    \boldsymbol{a}_3 &=& (0, &0, &a)
    \end{matrix}

Cell standardization
====================

No standardization is required.

Variations
==========

There are no variations for cubic lattice.
One example is predefined: ``cub`` with :math:`a = \pi`.

K-path
======

:math:`\mathrm{\Gamma-X-M-\Gamma-R-X\vert M-R}`.

=======================  ==============================  ==============================  ==============================
Point                    :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=======================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`  :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{M}`       :math:`1/2`                     :math:`1/2`                     :math:`0`
:math:`\mathrm{R}`       :math:`1/2`                     :math:`1/2`                     :math:`1/2`
:math:`\mathrm{X}`       :math:`0`                       :math:`1/2`                     :math:`0`
=======================  ==============================  ==============================  ==============================

Examples
========

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: cub_brillouin.py
    :language: py

.. raw:: html
    :file: cub_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: cub_real.py
    :language: py

.. raw:: html
    :file: cub_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: cub_wigner-seitz.py
    :language: py

.. raw:: html
    :file: cub_wigner-seitz.html

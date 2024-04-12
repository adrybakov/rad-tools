.. _guide_hex:

***************
Hexagonal (HEX)
***************

**Pearson symbol**: hP

**Constructor**:  :py:func:`.HEX`

It is defined by two parameter: :math:`a` and :math:`c`
with primitive and conventional lattice:

.. math::

    \begin{matrix}
    \boldsymbol{a}_1 &=& (a/2, &-a\sqrt{3}, &0)\\
    \boldsymbol{a}_2 &=& (a/2, &a\sqrt{3}, &0)\\
    \boldsymbol{a}_3 &=& (0, &0, &c)
    \end{matrix}

Cell standardization
====================

Angle between first two lattice vectors (:math:`\gamma`) has to be equal to :math:`120^{\circ}`.

If this condition is not satisfied, then the lattice is transformed to the standard form:


* If :math:`\beta = 120^{\circ}`
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_3, \boldsymbol{a}_1, \boldsymbol{a}_2)

* If :math:`\alpha = 120^{\circ}`
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_2, \boldsymbol{a}_3, \boldsymbol{a}_1)

K-path
======

:math:`\mathrm{\Gamma-M-K-\Gamma-A-L-H-A\vert L-M\vert K-H}`

=========================  ==============================  ==============================  ==============================
Point                      :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=========================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`    :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{A}`         :math:`0`                       :math:`0`                       :math:`1/2`
:math:`\mathrm{H}`         :math:`1/3`                     :math:`1/3`                     :math:`1/2`
:math:`\mathrm{K}`         :math:`1/3`                     :math:`1/3`                     :math:`0`
:math:`\mathrm{L}`         :math:`1/2`                     :math:`0`                       :math:`1/2`
:math:`\mathrm{M}`         :math:`1/2`                     :math:`0`                       :math:`0`
=========================  ==============================  ==============================  ==============================

Variations
==========

There are no variations for hexagonal lattice.
One example is predefined: ``hex`` with :math:`a = \pi` and :math:`c = 2\pi`.

Examples
========

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: hex_brillouin.py
    :language: py

.. raw:: html
    :file: hex_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: hex_real.py
    :language: py

.. raw:: html
    :file: hex_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: hex_wigner-seitz.py
    :language: py

.. raw:: html
    :file: hex_wigner-seitz.html

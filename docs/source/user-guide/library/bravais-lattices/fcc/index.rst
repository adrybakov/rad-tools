.. _guide_fcc:

************************
Face-centred cubic (FCC)
************************

**Pearson symbol**: cF

**Constructor**:  :py:func:`.FCC`

It is defined by one parameter: :math:`a` with conventional lattice:

.. math::

    \begin{matrix}
    \boldsymbol{a}_1 &=& (a, &0, &0)\\
    \boldsymbol{a}_2 &=& (0, &a, &0)\\
    \boldsymbol{a}_3 &=& (0, &0, &a)
    \end{matrix}

And primitive lattice:

.. math::

    \begin{matrix}
    \boldsymbol{a}_1 &=& (0, &a/2, &a/2)\\
    \boldsymbol{a}_2 &=& (a/2, &0, &a/2)\\
    \boldsymbol{a}_3 &=& (a/2, &a/2, &0)
    \end{matrix}



Cell standardization
====================

No standardization is performed.

Variations
==========

There are no variations for face-centered cubic lattice.
One example is predefined: ``fcc`` with :math:`a = \pi`.

K-path
======

:math:`\mathrm{\Gamma-X-W-K-\Gamma-L-U-W-L-K\vert U-X}`

=======================  ==============================  ==============================  ==============================
Point                    :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=======================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`  :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{K}`       :math:`3/8`                     :math:`3/8`                     :math:`3/4`
:math:`\mathrm{L}`       :math:`1/2`                     :math:`1/2`                     :math:`1/2`
:math:`\mathrm{U}`       :math:`5/8`                     :math:`1/4`                     :math:`5/8`
:math:`\mathrm{W}`       :math:`1/2`                     :math:`1/4`                     :math:`3/4`
:math:`\mathrm{X}`       :math:`1/2`                     :math:`0`                       :math:`1/2`
=======================  ==============================  ==============================  ==============================

Examples
========

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: fcc_brillouin.py
    :language: py

.. raw:: html
    :file: fcc_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: fcc_real.py
    :language: py

.. raw:: html
    :file: fcc_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: fcc_wigner-seitz.py
    :language: py

.. raw:: html
    :file: fcc_wigner-seitz.html

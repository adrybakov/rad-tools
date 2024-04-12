.. _guide_bcc:

*************************
Body-centered cubic (BCC)
*************************

**Pearson symbol**: cI

**Constructor**:  :py:func:`.BCC`

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
    \boldsymbol{a}_1 &=& (-a/2,& a/2,& a/2)\\
    \boldsymbol{a}_2 &=& (a/2, &-a/2,& a/2)\\
    \boldsymbol{a}_3 &=& (a/2, &a/2, &-a/2)
    \end{matrix}

Cell standardization
====================

No standardization is performed.

K-path
======

:math:`\mathrm{\Gamma-H-N-\Gamma-P-H\vert P-N}`

=======================  ==============================  ==============================  ==============================
Point                    :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=======================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`  :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{H}`       :math:`1/2`                     :math:`-1/2`                    :math:`1/2`
:math:`\mathrm{P}`       :math:`1/4`                     :math:`1/4`                     :math:`1/4`
:math:`\mathrm{N}`       :math:`0`                       :math:`0`                       :math:`1/2`
=======================  ==============================  ==============================  ==============================

Variations
==========

There are no variations for body-centered cubic lattice.
One example is predefined: ``bcc`` with :math:`a = \pi`.

Examples
========

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: bcc_brillouin.py
    :language: py

.. raw:: html
    :file: bcc_brillouin.html


Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: bcc_real.py
    :language: py

.. raw:: html
    :file: bcc_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: bcc_wigner-seitz.py
    :language: py

.. raw:: html
    :file: bcc_wigner-seitz.html

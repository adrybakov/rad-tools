.. _guide_bct:

*****************************
Body-centred tetragonal (BCT)
*****************************

**Pearson symbol**: tI

**Constructor**:  :py:func:`.BCT`

It is defined by two parameters: :math:`a` and :math:`c`
with conventional lattice:

.. math::

    \begin{matrix}
    \boldsymbol{a}_1 &=& (a, &0, &0)\\
    \boldsymbol{a}_2 &=& (0, &a, &0)\\
    \boldsymbol{a}_3 &=& (0, &0, &c)
    \end{matrix}

And primitive lattice:

.. math::

    \begin{matrix}
    \boldsymbol{a}_1 &=& (-a/2, &a/2, &c/2)\\
    \boldsymbol{a}_2 &=& (a/2, &-a/2, &c/2)\\
    \boldsymbol{a}_3 &=& (a/2, &a/2, &-c/2)
    \end{matrix}

Cell standardization
====================

Length of third lattice vector of the conventional cell has to be different from the first two.
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

BCT\ :sub:`1`
-------------

:math:`\mathrm{\Gamma-X-M-\Gamma-Z-P-N-Z_1-M\vert X-P}`

.. math::

    \eta = \dfrac{1 + c^2/a^2}{4}

=======================  ==============================  ==============================  ==============================
Point                    :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=======================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`  :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{M}`       :math:`-1/2`                    :math:`1/2`                     :math:`1/2`
:math:`\mathrm{N}`       :math:`0`                       :math:`1/2`                     :math:`0`
:math:`\mathrm{P}`       :math:`1/4`                     :math:`1/4`                     :math:`1/4`
:math:`\mathrm{X}`       :math:`0`                       :math:`0`                       :math:`1/2`
:math:`\mathrm{Z}`       :math:`\eta`                    :math:`\eta`                    :math:`-\eta`
:math:`\mathrm{Z}_1`     :math:`-\eta`                   :math:`1-\eta`                  :math:`\eta`
=======================  ==============================  ==============================  ==============================

BCT\ :sub:`2`
-------------

:math:`\mathrm{\Gamma-X-Y-\Sigma-\Gamma-Z-\Sigma_1-N-P-Y_1-Z\vert X-P}`

.. math::

    \begin{matrix}
    \eta = \dfrac{1 + a^2/c^2}{4} &
    \zeta = \dfrac{a^2}{2c^2}
    \end{matrix}

=========================  ==============================  ==============================  ==============================
Point                      :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=========================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`    :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{N}`         :math:`0`                       :math:`1/2`                     :math:`0`
:math:`\mathrm{P}`         :math:`1/4`                     :math:`1/4`                     :math:`1/4`
:math:`\mathrm{\Sigma}`    :math:`-\eta`                   :math:`\eta`                    :math:`\eta`
:math:`\mathrm{\Sigma_1}`  :math:`\eta`                    :math:`1-\eta`                  :math:`-\eta`
:math:`\mathrm{X}`         :math:`0`                       :math:`0`                       :math:`1/2`
:math:`\mathrm{Y}`         :math:`-\zeta`                  :math:`\zeta`                   :math:`1/2`
:math:`\mathrm{Y}_1`       :math:`1/2`                     :math:`1/2`                     :math:`-\zeta`
:math:`\mathrm{Z}`         :math:`1/2`                     :math:`1/2`                     :math:`-1/2`
=========================  ==============================  ==============================  ==============================


Variations
==========

There are two variations of body-centered tetragonal lattice.

BCT\ :sub:`1`
-------------

:math:`c < a`.

Predefined example: ``bct1`` with :math:`a = 1.5\pi` and :math:`c = \pi`.

BCT\ :sub:`2`
-------------

:math:`c > a`.

Predefined example: ``bct2`` with :math:`a = \pi` and :math:`c = 1.5\pi`.


Examples
========

BCT\ :sub:`1`
-------------

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: bct1_brillouin.py
    :language: py

.. raw:: html
    :file: bct1_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: bct1_real.py
    :language: py

.. raw:: html
    :file: bct1_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: bct1_wigner-seitz.py
    :language: py

.. raw:: html
    :file: bct1_wigner-seitz.html

BCT\ :sub:`2`
-------------


Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: bct2_brillouin.py
    :language: py

.. raw:: html
    :file: bct2_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: bct2_real.py
    :language: py

.. raw:: html
    :file: bct2_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: bct2_wigner-seitz.py
    :language: py

.. raw:: html
    :file: bct2_wigner-seitz.html

Edge cases
==========

If :math:`a = c` then the lattice is :ref:`guide_bcc`.

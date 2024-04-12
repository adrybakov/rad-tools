.. _guide_tri:

***************
Triclinic (TRI)
***************

**Pearson symbol**: aP

**Constructor**:  :py:func:`.TRI`

It is defined by six parameters: :math:`a`, :math:`b`, :math:`c` and
:math:`\alpha`, :math:`\beta`, :math:`\gamma`.
with primitive and conventional lattice:

.. math::

    \begin{matrix}
    \boldsymbol{a}_1 = (a, 0, 0)\\
    \boldsymbol{a}_2 = (b\cos\gamma, b\sin\gamma, 0)\\
    \boldsymbol{a}_3 = (c\cos\beta, \dfrac{c(\cos\alpha - \cos\beta\cos\gamma)}{\sin\gamma}, \dfrac{c}{\sin\gamma}\sqrt{\sin^2\gamma - \cos^2\alpha - \cos^2\beta + 2\cos\alpha\cos\beta\cos\gamma})
    \end{matrix}

Cell standardization
====================

Standardization is performed based on the reciprocal cell.

One of the conditions have to be met:

* All reciprocal cell angles (:math:`k_{\alpha}`, :math:`k_{\beta}`, :math:`k_{\gamma}`) are :math:`> 90^{\circ}` and :math:`k_{\gamma} = \min(k_{\alpha}, k_{\beta}, k_{\gamma})`.

* All reciprocal cell angles (:math:`k_{\alpha}`, :math:`k_{\beta}`, :math:`k_{\gamma}`) are :math:`< 90^{\circ}` and :math:`k_{\gamma} = \max(k_{\alpha}, k_{\beta}, k_{\gamma})`.

* :math:`k_{\gamma} = 90^{\circ}` and other two angles are :math:`> 90^{\circ}`.

* :math:`k_{\gamma} = 90^{\circ}` and other two angles are :math:`< 90^{\circ}`.

If these conditions are not satisfied, then the lattice is transformed to the
standard form:

Last two cases
--------------
First we check for the last two cases:

* If :math:`k_{\alpha} = 90^{\circ}`
    .. math::

        (\boldsymbol{b}_1, \boldsymbol{b}_2, \boldsymbol{b}_3) \rightarrow
        (\boldsymbol{b}_2, \boldsymbol{b}_3, \boldsymbol{b}_1)

* If :math:`k_{\beta} = 90^{\circ}`
    .. math::

        (\boldsymbol{b}_1, \boldsymbol{b}_2, \boldsymbol{b}_3) \rightarrow
        (\boldsymbol{b}_3, \boldsymbol{b}_1, \boldsymbol{b}_2)

If one of the last two conditions were met, then now we have :math:`k_{\gamma} = 90^{\circ}`.
We need to choose appropriate values for the remaining two angles:

* If :math:`k_{\alpha} > 90^{\circ}` and :math:`k_{\beta} < 90^{\circ}` or :math:`k_{\alpha} < 90^{\circ}` and :math:`k_{\beta} > 90^{\circ}`:
    .. math::

        (\boldsymbol{b}_1, \boldsymbol{b}_2, \boldsymbol{b}_3) \rightarrow
        (\boldsymbol{b}_2, -\boldsymbol{b}_1, \boldsymbol{b}_3)

First two cases
---------------

If none of the last two conditions were met, then we check for the first two.
First we ensure that all angles are :math:`> 90^{\circ}` or :math:`< 90^{\circ}`:


* If :math:`k_{\alpha} > 90^{\circ}` and :math:`k_{\beta} > 90^{\circ}` and :math:`k_{\gamma} < 90^{\circ}`:
    .. math::

        (\boldsymbol{b}_1, \boldsymbol{b}_2, \boldsymbol{b}_3) \rightarrow
        (-\boldsymbol{b}_1, -\boldsymbol{b}_2, \boldsymbol{b}_3)

* If :math:`k_{\alpha} > 90^{\circ}` and :math:`k_{\beta} < 90^{\circ}` and :math:`k_{\gamma} > 90^{\circ}`:
    .. math::

        (\boldsymbol{b}_1, \boldsymbol{b}_2, \boldsymbol{b}_3) \rightarrow
        (-\boldsymbol{b}_1, \boldsymbol{b}_2, -\boldsymbol{b}_3)

* If :math:`k_{\alpha} > 90^{\circ}` and :math:`k_{\beta} < 90^{\circ}` and :math:`k_{\gamma} < 90^{\circ}`:
    .. math::

        (\boldsymbol{b}_1, \boldsymbol{b}_2, \boldsymbol{b}_3) \rightarrow
        (\boldsymbol{b}_1, -\boldsymbol{b}_2, -\boldsymbol{b}_3)

* If :math:`k_{\alpha} < 90^{\circ}` and :math:`k_{\beta} > 90^{\circ}` and :math:`k_{\gamma} > 90^{\circ}`:
    .. math::

        (\boldsymbol{b}_1, \boldsymbol{b}_2, \boldsymbol{b}_3) \rightarrow
        (\boldsymbol{b}_1, -\boldsymbol{b}_2, -\boldsymbol{b}_3)

* If :math:`k_{\alpha} < 90^{\circ}` and :math:`k_{\beta} > 90^{\circ}` and :math:`k_{\gamma} < 90^{\circ}`:
    .. math::

        (\boldsymbol{b}_1, \boldsymbol{b}_2, \boldsymbol{b}_3) \rightarrow
        (-\boldsymbol{b}_1, \boldsymbol{b}_2, -\boldsymbol{b}_3)

* If :math:`k_{\alpha} < 90^{\circ}` and :math:`k_{\beta} < 90^{\circ}` and :math:`k_{\gamma} > 90^{\circ}`:
    .. math::

        (\boldsymbol{b}_1, \boldsymbol{b}_2, \boldsymbol{b}_3) \rightarrow
        (-\boldsymbol{b}_1, -\boldsymbol{b}_2, \boldsymbol{b}_3)

As the last step we reorder the reciprocal vectors:

Reordering if all angles are :math:`> 90^{\circ}`:

* If :math:`k_{\alpha} = min(k_{\alpha}, k_{\beta}, k_{\gamma})`, then:
    .. math::

        (\boldsymbol{b}_1, \boldsymbol{b}_2, \boldsymbol{b}_3) \rightarrow
        (\boldsymbol{b}_2, \boldsymbol{b}_3, \boldsymbol{b}_1)

* If :math:`k_{\beta} = min(k_{\alpha}, k_{\beta}, k_{\gamma})`, then:
    .. math::

        (\boldsymbol{b}_1, \boldsymbol{b}_2, \boldsymbol{b}_3) \rightarrow
        (\boldsymbol{b}_3, \boldsymbol{b}_1, \boldsymbol{b}_2)


Reordering if all angles are :math:`< 90^{\circ}`:

* If :math:`k_{\alpha} = max(k_{\alpha}, k_{\beta}, k_{\gamma})`, then:
    .. math::

        (\boldsymbol{b}_1, \boldsymbol{b}_2, \boldsymbol{b}_3) \rightarrow
        (\boldsymbol{b}_2, \boldsymbol{b}_3, \boldsymbol{b}_1)

* If :math:`k_{\beta} = max(k_{\alpha}, k_{\beta}, k_{\gamma})`, then:
    .. math::

        (\boldsymbol{b}_1, \boldsymbol{b}_2, \boldsymbol{b}_3) \rightarrow
        (\boldsymbol{b}_3, \boldsymbol{b}_1, \boldsymbol{b}_2)


K-path
======

TRI\ :sub:`1a`
--------------

:math:`\mathrm{X-\Gamma-Y\vert L-\Gamma-Z\vert N-\Gamma-M\vert R-\Gamma}`

=======================  ==============================  ==============================  ==============================
Point                    :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=======================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`  :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{L}`       :math:`1/2`                     :math:`1/2`                     :math:`0`
:math:`\mathrm{M}`       :math:`0`                       :math:`1/2`                     :math:`1/2`
:math:`\mathrm{N}`       :math:`1/2`                     :math:`0`                       :math:`1/2`
:math:`\mathrm{R}`       :math:`1/2`                     :math:`1/2`                     :math:`1/2`
:math:`\mathrm{X}`       :math:`1/2`                     :math:`0`                       :math:`0`
:math:`\mathrm{Y}`       :math:`0`                       :math:`1/2`                     :math:`0`
:math:`\mathrm{Z}`       :math:`0`                       :math:`0`                       :math:`1/2`
=======================  ==============================  ==============================  ==============================

TRI\ :sub:`2a`
--------------

:math:`\mathrm{X-\Gamma-Y\vert L-\Gamma-Z\vert N-\Gamma-M\vert R-\Gamma}`

=======================  ==============================  ==============================  ==============================
Point                    :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=======================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`  :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{L}`       :math:`1/2`                     :math:`1/2`                     :math:`0`
:math:`\mathrm{M}`       :math:`0`                       :math:`1/2`                     :math:`1/2`
:math:`\mathrm{N}`       :math:`1/2`                     :math:`0`                       :math:`1/2`
:math:`\mathrm{R}`       :math:`1/2`                     :math:`1/2`                     :math:`1/2`
:math:`\mathrm{X}`       :math:`1/2`                     :math:`0`                       :math:`0`
:math:`\mathrm{Y}`       :math:`0`                       :math:`1/2`                     :math:`0`
:math:`\mathrm{Z}`       :math:`0`                       :math:`0`                       :math:`1/2`
=======================  ==============================  ==============================  ==============================

TRI\ :sub:`1b`
--------------

:math:`\mathrm{X-\Gamma-Y\vert L-\Gamma-Z\vert N-\Gamma-M\vert R-\Gamma}`

=======================  ==============================  ==============================  ==============================
Point                    :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=======================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`  :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{L}`       :math:`1/2`                     :math:`-1/2`                    :math:`0`
:math:`\mathrm{M}`       :math:`0`                       :math:`0`                       :math:`1/2`
:math:`\mathrm{N}`       :math:`-1/2`                    :math:`-1/2`                    :math:`1/2`
:math:`\mathrm{R}`       :math:`0`                       :math:`-1/2`                    :math:`1/2`
:math:`\mathrm{X}`       :math:`0`                       :math:`-1/2`                    :math:`0`
:math:`\mathrm{Y}`       :math:`1/2`                     :math:`0`                       :math:`0`
:math:`\mathrm{Z}`       :math:`-1/2`                    :math:`0`                       :math:`1/2`
=======================  ==============================  ==============================  ==============================

TRI\ :sub:`2b`
--------------

:math:`\mathrm{X-\Gamma-Y\vert L-\Gamma-Z\vert N-\Gamma-M\vert R-\Gamma}`

=======================  ==============================  ==============================  ==============================
Point                    :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=======================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`  :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{L}`       :math:`1/2`                     :math:`-1/2`                    :math:`0`
:math:`\mathrm{M}`       :math:`0`                       :math:`0`                       :math:`1/2`
:math:`\mathrm{N}`       :math:`-1/2`                    :math:`-1/2`                    :math:`1/2`
:math:`\mathrm{R}`       :math:`0`                       :math:`-1/2`                    :math:`1/2`
:math:`\mathrm{X}`       :math:`0`                       :math:`-1/2`                    :math:`0`
:math:`\mathrm{Y}`       :math:`1/2`                     :math:`0`                       :math:`0`
:math:`\mathrm{Z}`       :math:`-1/2`                    :math:`0`                       :math:`1/2`
=======================  ==============================  ==============================  ==============================

Variations
==========

There are four variations of triclinic lattice.

TRI\ :sub:`1a`
--------------

:math:`k_{\alpha} > 90^{\circ}, k_{\beta} > 90^{\circ}, k_{\gamma} > 90^{\circ}, k_{\gamma} = \min(k_{\alpha}, k_{\beta}, k_{\gamma})`

TRI\ :sub:`2a`
--------------

:math:`k_{\alpha} > 90^{\circ}, k_{\beta} > 90^{\circ}, k_{\gamma} = 90^{\circ}`

TRI\ :sub:`1b`
--------------

:math:`k_{\alpha} < 90^{\circ}, k_{\beta} < 90^{\circ}, k_{\gamma} < 90^{\circ}, k_{\gamma} = \max(k_{\alpha}, k_{\beta}, k_{\gamma})`

TRI\ :sub:`2b`
--------------

:math:`k_{\alpha} < 90^{\circ}, k_{\beta} < 90^{\circ}, k_{\gamma} = 90^{\circ}`

In definition of the examples we cheated and defined them through reciprocal lattice parameters.

Examples
========

TRI\ :sub:`1a`
--------------

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: tri1a_brillouin.py
    :language: py

.. raw:: html
    :file: tri1a_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: tri1a_real.py
    :language: py

.. raw:: html
    :file: tri1a_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: tri1a_wigner-seitz.py
    :language: py

.. raw:: html
    :file: tri1a_wigner-seitz.html

TRI\ :sub:`2a`
--------------

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: tri2a_brillouin.py
    :language: py

.. raw:: html
    :file: tri2a_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: tri2a_real.py
    :language: py

.. raw:: html
    :file: tri2a_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: tri2a_wigner-seitz.py
    :language: py

.. raw:: html
    :file: tri2a_wigner-seitz.html

TRI\ :sub:`1b`
--------------

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: tri1b_brillouin.py
    :language: py

.. raw:: html
    :file: tri1b_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: tri1b_real.py
    :language: py

.. raw:: html
    :file: tri1b_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: tri1b_wigner-seitz.py
    :language: py

.. raw:: html
    :file: tri1b_wigner-seitz.html

TRI\ :sub:`2b`
--------------

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: tri2b_brillouin.py
    :language: py

.. raw:: html
    :file: tri2b_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: tri2b_real.py
    :language: py

.. raw:: html
    :file: tri2b_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: tri2b_wigner-seitz.py
    :language: py

.. raw:: html
    :file: tri2b_wigner-seitz.html

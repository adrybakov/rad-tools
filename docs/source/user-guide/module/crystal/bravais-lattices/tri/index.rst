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

TODO. At the moment no standardization is performed.

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

There are four variation of triclinic lattice.

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

* Brillouin zone and default kpath

.. literalinclude:: tri1a_brillouin.py
    :language: py

.. figure:: tri1a_brillouin.png 
    :target: ../../../../../_images/tri1a_brillouin.png 

* Primitive and conventional cell

.. literalinclude:: tri1a_real.py
    :language: py

.. figure:: tri1a_real.png 
    :target: ../../../../../_images/tri1a_real.png 

* Wigner-Seitz cell

.. literalinclude:: tri1a_wigner-seitz.py
    :language: py

.. figure:: tri1a_wigner-seitz.png 
    :target: ../../../../../_images/tri1a_wigner-seitz.png 

TRI\ :sub:`2a`
--------------

* Brillouin zone and default kpath

.. literalinclude:: tri2a_brillouin.py
    :language: py

.. figure:: tri2a_brillouin.png 
    :target: ../../../../../_images/tri2a_brillouin.png 

* Primitive and conventional cell

.. literalinclude:: tri2a_real.py
    :language: py

.. figure:: tri2a_real.png 
    :target: ../../../../../_images/tri2a_real.png 

* Wigner-Seitz cell

.. literalinclude:: tri2a_wigner-seitz.py
    :language: py

.. figure:: tri2a_wigner-seitz.png 
    :target: ../../../../../_images/tri2a_wigner-seitz.png 

TRI\ :sub:`1b`
--------------

* Brillouin zone and default kpath

.. literalinclude:: tri1b_brillouin.py
    :language: py

.. figure:: tri1b_brillouin.png 
    :target: ../../../../../_images/tri1b_brillouin.png 

* Primitive and conventional cell

.. literalinclude:: tri1b_real.py
    :language: py

.. figure:: tri1b_real.png 
    :target: ../../../../../_images/tri1b_real.png 

* Wigner-Seitz cell

.. literalinclude:: tri1b_wigner-seitz.py
    :language: py

.. figure:: tri1b_wigner-seitz.png 
    :target: ../../../../../_images/tri1b_wigner-seitz.png 

TRI\ :sub:`2b`
--------------

* Brillouin zone and default kpath

.. literalinclude:: tri2b_brillouin.py
    :language: py

.. figure:: tri2b_brillouin.png 
    :target: ../../../../../_images/tri2b_brillouin.png 

* Primitive and conventional cell

.. literalinclude:: tri2b_real.py
    :language: py

.. figure:: tri2b_real.png 
    :target: ../../../../../_images/tri2b_real.png 

* Wigner-Seitz cell

.. literalinclude:: tri2b_wigner-seitz.py
    :language: py

.. figure:: tri2b_wigner-seitz.png 
    :target: ../../../../../_images/tri2b_wigner-seitz.png 

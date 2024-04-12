.. _guide_mclc:

******************************
Base-centred monoclinic (MCLC)
******************************

**Pearson symbol**: mS

**Constructor**:  :py:func:`.MCLC`

It is defined by four parameter: :math:`a`, :math:`b`, :math:`c` and :math:`\alpha`
with conventional lattice:

.. math::

    \begin{matrix}
    \boldsymbol{a}_1 &=& (a, &0, &0)\\
    \boldsymbol{a}_2 &=& (0, &b, &0)\\
    \boldsymbol{a}_3 &=& (0, &c\cos\alpha, &c\sin\alpha)
    \end{matrix}

And primitive lattice:

.. math::

    \begin{matrix}
    \boldsymbol{a}_1 &=& (a/2, &b/2, &0)\\
    \boldsymbol{a}_2 &=& (-a/2, &b/2, &0)\\
    \boldsymbol{a}_3 &=& (0, &c\cos\alpha, &c\sin\alpha)
    \end{matrix}

Order of parameters: :math:`b \le c`, :math:`\alpha < 90^{\circ}`.

Cell standardization
====================

Length of the third vector of the primitive lattice has to be different from the
length of the other two lattice vectors.

Angle between second and third lattice vectors of conventional lattice (:math:`\alpha`)
has to be less then :math:`90^{\circ}`.

Length of the second lattice vector of the conventional lattice has to be less or
equal to the length of the third lattice vector.

If these conditions are not satisfied, then the lattice is transformed to the standard form:

First we ensure the length of the third vector is different from the length of the other two vectors.
For this step we use vectors of the primitive lattice:

* If :math:`\vert\boldsymbol{a}_1\vert = \vert\boldsymbol{a}_3\vert`:
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_3, \boldsymbol{a}_1, \boldsymbol{a}_2)

* If :math:`\vert\boldsymbol{a}_2\vert = \vert\boldsymbol{a}_3\vert`:
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_2, \boldsymbol{a}_3, \boldsymbol{a}_1)

Then we ensure the :math:`\alpha < 90^{\circ}`. For this step we use vectors of the conventional lattice:

* If :math:`\alpha > 90^{\circ}`:
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (\boldsymbol{a}_1, \boldsymbol{a}_3, -\boldsymbol{a}_2)

Finally, we ensure the :math:`b \le c`. For this step we use vectors of the conventional lattice:

* If :math:`b > c`:
    .. math::

        (\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3) \rightarrow
        (-\boldsymbol{a}_1, \boldsymbol{a}_3, \boldsymbol{a}_2)

.. note::

    First and second lattice vectors are multiplied by :math:`-1` in some cases to
    preserve the handedness of the cell.

K-path
======

MCLC\ :sub:`1`
--------------

:math:`\mathrm{\Gamma-Y-F-L-I\vert I_1-Z-F_1\vert Y-X_1\vert X-\Gamma-N\vert M-\Gamma}`

.. math::

    \begin{matrix}
    \zeta = \dfrac{2 - b\cos\alpha/c}{4\sin^2\alpha} &
    \eta = \dfrac{1}{2} + \dfrac{2\zeta c\cos\alpha}{b} \\
    \psi = \dfrac{3}{4} - \dfrac{a^2}{4b^2\sin^2\alpha} &
    \phi = \psi + \dfrac{(3/4-\psi)b\cos\alpha}{c}
    \end{matrix}

=========================  ==============================  ==============================  ==============================
Point                      :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=========================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`    :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{N}`         :math:`1/2`                     :math:`0`                       :math:`0`
:math:`\mathrm{N_1}`       :math:`0`                       :math:`-1/2`                    :math:`0`
:math:`\mathrm{F}`         :math:`1-\zeta`                 :math:`1-\zeta`                 :math:`1-\eta`
:math:`\mathrm{F_1}`       :math:`\zeta`                   :math:`\zeta`                   :math:`\eta`
:math:`\mathrm{F_2}`       :math:`-\zeta`                  :math:`-\zeta`                  :math:`1-\eta`
:math:`\mathrm{F_3}`       :math:`1-\zeta`                 :math:`-\zeta`                  :math:`1-\eta`
:math:`\mathrm{I}`         :math:`\phi`                    :math:`1-\phi`                  :math:`1/2`
:math:`\mathrm{I_1}`       :math:`1-\phi`                  :math:`\phi-1`                  :math:`1/2`
:math:`\mathrm{L}`         :math:`1/2`                     :math:`1/2`                     :math:`1/2`
:math:`\mathrm{M}`         :math:`1/2`                     :math:`0`                       :math:`1/2`
:math:`\mathrm{X}`         :math:`1-\psi`                  :math:`\psi-1`                  :math:`0`
:math:`\mathrm{X_1}`       :math:`\psi`                    :math:`1-\psi`                  :math:`0`
:math:`\mathrm{X_2}`       :math:`\psi-1`                  :math:`-\psi`                   :math:`0`
:math:`\mathrm{Y}`         :math:`1/2`                     :math:`1/2`                     :math:`0`
:math:`\mathrm{Y_1}`       :math:`-1/2`                    :math:`-1/2`                    :math:`0`
:math:`\mathrm{Z}`         :math:`0`                       :math:`0`                       :math:`1/2`
=========================  ==============================  ==============================  ==============================


MCLC\ :sub:`2`
--------------

:math:`\mathrm{\Gamma-Y-F-L-I\vert I_1-Z-F_1\vert N-\Gamma-M}`

.. math::

    \begin{matrix}
    \zeta = \dfrac{2 - b\cos\alpha/c}{4\sin^2\alpha} &
    \eta = \dfrac{1}{2} + \dfrac{2\zeta c\cos\alpha}{b} \\
    \psi = \dfrac{3}{4} - \dfrac{a^2}{4b^2\sin^2\alpha} &
    \phi = \psi + \dfrac{(3/4-\psi)b\cos\alpha}{c}
    \end{matrix}

=========================  ==============================  ==============================  ==============================
Point                      :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=========================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`    :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{N}`         :math:`1/2`                     :math:`0`                       :math:`0`
:math:`\mathrm{N_1}`       :math:`0`                       :math:`-1/2`                    :math:`0`
:math:`\mathrm{F}`         :math:`1-\zeta`                 :math:`1-\zeta`                 :math:`1-\eta`
:math:`\mathrm{F_1}`       :math:`\zeta`                   :math:`\zeta`                   :math:`\eta`
:math:`\mathrm{F_2}`       :math:`-\zeta`                  :math:`-\zeta`                  :math:`1-\eta`
:math:`\mathrm{F_3}`       :math:`1-\zeta`                 :math:`-\zeta`                  :math:`1-\eta`
:math:`\mathrm{I}`         :math:`\phi`                    :math:`1-\phi`                  :math:`1/2`
:math:`\mathrm{I_1}`       :math:`1-\phi`                  :math:`\phi-1`                  :math:`1/2`
:math:`\mathrm{L}`         :math:`1/2`                     :math:`1/2`                     :math:`1/2`
:math:`\mathrm{M}`         :math:`1/2`                     :math:`0`                       :math:`1/2`
:math:`\mathrm{X}`         :math:`1-\psi`                  :math:`\psi-1`                  :math:`0`
:math:`\mathrm{X_1}`       :math:`\psi`                    :math:`1-\psi`                  :math:`0`
:math:`\mathrm{X_2}`       :math:`\psi-1`                  :math:`-\psi`                   :math:`0`
:math:`\mathrm{Y}`         :math:`1/2`                     :math:`1/2`                     :math:`0`
:math:`\mathrm{Y_1}`       :math:`-1/2`                    :math:`-1/2`                    :math:`0`
:math:`\mathrm{Z}`         :math:`0`                       :math:`0`                       :math:`1/2`
=========================  ==============================  ==============================  ==============================

MCLC\ :sub:`3`
--------------

:math:`\mathrm{\Gamma-Y-F-H-Z-I-F_1\vert H_1-Y_1-X-\Gamma-N\vert M-\Gamma}`

.. math::

    \begin{matrix}
    \mu = \dfrac{1+b^2/a^2}{4} &
    \delta = \dfrac{bc\cos\alpha}{2a^2} &
    \zeta = \mu - \dfrac{1}{4} + \dfrac{1 - b\cos\alpha/c}{4\sin^2\alpha} \\
    \eta = \dfrac{1}{2} + \dfrac{2\zeta c \cos\alpha}{b} &
    \phi = 1 + \zeta - 2\mu &
    \psi = \eta - 2\delta
    \end{matrix}

=========================  ==============================  ==============================  ==============================
Point                      :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=========================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`    :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{F}`         :math:`1-\phi`                  :math:`1-\phi`                  :math:`1-\psi`
:math:`\mathrm{F_1}`       :math:`\phi`                    :math:`\phi-1`                  :math:`\psi`
:math:`\mathrm{F_2}`       :math:`1-\phi`                  :math:`-\phi`                   :math:`1-\psi`
:math:`\mathrm{H}`         :math:`\zeta`                   :math:`\zeta`                   :math:`\eta`
:math:`\mathrm{H_1}`       :math:`1-\zeta`                 :math:`-\zeta`                  :math:`1-\eta`
:math:`\mathrm{H_2}`       :math:`-\zeta`                  :math:`-\zeta`                  :math:`1-\eta`
:math:`\mathrm{I}`         :math:`1/2`                     :math:`-1/2`                    :math:`1/2`
:math:`\mathrm{M}`         :math:`1/2`                     :math:`0`                       :math:`1/2`
:math:`\mathrm{N}`         :math:`1/2`                     :math:`0`                       :math:`0`
:math:`\mathrm{N_1}`       :math:`0`                       :math:`-1/2`                    :math:`0`
:math:`\mathrm{X}`         :math:`1/2`                     :math:`-1/2`                    :math:`0`
:math:`\mathrm{Y}`         :math:`\mu`                     :math:`\mu`                     :math:`\delta`
:math:`\mathrm{Y_1}`       :math:`1-\mu`                   :math:`-\mu`                    :math:`-\delta`
:math:`\mathrm{Y_2}`       :math:`-\mu`                    :math:`-\mu`                    :math:`-\delta`
:math:`\mathrm{Y_3}`       :math:`\mu`                     :math:`\mu-1`                   :math:`\delta`
:math:`\mathrm{Z}`         :math:`0`                       :math:`0`                       :math:`1/2`
=========================  ==============================  ==============================  ==============================


MCLC\ :sub:`4`
--------------

:math:`\mathrm{\Gamma-Y-F-H-Z-I\vert H_1-Y_1-X-\Gamma-N\vert M-\Gamma}`

.. math::

    \begin{matrix}
    \mu = \dfrac{1+b^2/a^2}{4} &
    \delta = \dfrac{bc\cos\alpha}{2a^2} &
    \zeta = \mu - \dfrac{1}{4} + \dfrac{1 - b\cos\alpha/c}{4\sin^2\alpha} \\
    \eta = \dfrac{1}{2} + \dfrac{2\zeta c \cos\alpha}{b} &
    \phi = 1 + \zeta - 2\mu &
    \psi = \eta - 2\delta
    \end{matrix}

=========================  ==============================  ==============================  ==============================
Point                      :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=========================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`    :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{F}`         :math:`1-\phi`                  :math:`1-\phi`                  :math:`1-\psi`
:math:`\mathrm{F_1}`       :math:`\phi`                    :math:`\phi-1`                  :math:`\psi`
:math:`\mathrm{F_2}`       :math:`1-\phi`                  :math:`-\phi`                   :math:`1-\psi`
:math:`\mathrm{H}`         :math:`\zeta`                   :math:`\zeta`                   :math:`\eta`
:math:`\mathrm{H_1}`       :math:`1-\zeta`                 :math:`-\zeta`                  :math:`1-\eta`
:math:`\mathrm{H_2}`       :math:`-\zeta`                  :math:`-\zeta`                  :math:`1-\eta`
:math:`\mathrm{I}`         :math:`1/2`                     :math:`-1/2`                    :math:`1/2`
:math:`\mathrm{M}`         :math:`1/2`                     :math:`0`                       :math:`1/2`
:math:`\mathrm{N}`         :math:`1/2`                     :math:`0`                       :math:`0`
:math:`\mathrm{N_1}`       :math:`0`                       :math:`-1/2`                    :math:`0`
:math:`\mathrm{X}`         :math:`1/2`                     :math:`-1/2`                    :math:`0`
:math:`\mathrm{Y}`         :math:`\mu`                     :math:`\mu`                     :math:`\delta`
:math:`\mathrm{Y_1}`       :math:`1-\mu`                   :math:`-\mu`                    :math:`-\delta`
:math:`\mathrm{Y_2}`       :math:`-\mu`                    :math:`-\mu`                    :math:`-\delta`
:math:`\mathrm{Y_3}`       :math:`\mu`                     :math:`\mu-1`                   :math:`\delta`
:math:`\mathrm{Z}`         :math:`0`                       :math:`0`                       :math:`1/2`
=========================  ==============================  ==============================  ==============================



MCLC\ :sub:`5`
--------------

:math:`\mathrm{\Gamma-Y-F-L-I\vert I_1-Z-H-F_1\vert H_1-Y_1-X-\Gamma-N\vert M-\Gamma}`

.. math::

    \begin{matrix}
    \zeta = \dfrac{b^2}{4a^2} + \dfrac{1 - b\cos\alpha/c}{4\sin^2\alpha} &
    \eta = \dfrac{1}{2} + \dfrac{2\zeta c\cos\alpha}{b} \\
    \mu = \dfrac{\eta}{2} + \dfrac{b^2}{4a^2} - \dfrac{bc\cos\alpha}{2a^2} &
    \nu = 2\mu - \zeta \\
    \omega = \dfrac{(4\nu - 1 - b^2\sin^2\alpha/a^2)c}{2b\cos\alpha} &
    \delta = \dfrac{\zeta c \cos\alpha}{b} + \dfrac{\omega}{2} - \dfrac{1}{4} &
    \rho = 1 - \dfrac{\zeta a^2}{b^2}
    \end{matrix}

=========================  ==============================  ==============================  ==============================
Point                      :math:`\times\boldsymbol{b}_1`  :math:`\times\boldsymbol{b}_2`  :math:`\times\boldsymbol{b}_3`
=========================  ==============================  ==============================  ==============================
:math:`\mathrm{\Gamma}`    :math:`0`                       :math:`0`                       :math:`0`
:math:`\mathrm{F}`         :math:`\nu`                     :math:`\nu`                     :math:`\omega`
:math:`\mathrm{F_1}`       :math:`1-\nu`                   :math:`-\nu`                    :math:`1-\omega`
:math:`\mathrm{F_2}`       :math:`\nu`                     :math:`\nu-1`                   :math:`\omega`
:math:`\mathrm{H}`         :math:`\zeta`                   :math:`\zeta`                   :math:`\eta`
:math:`\mathrm{H_1}`       :math:`1-\zeta`                 :math:`-\zeta`                  :math:`1-\eta`
:math:`\mathrm{H_2}`       :math:`-\zeta`                  :math:`-\zeta`                  :math:`1-\eta`
:math:`\mathrm{I}`         :math:`\rho`                    :math:`1-\rho`                  :math:`1/2`
:math:`\mathrm{I_1}`       :math:`1-\rho`                  :math:`\rho-1`                  :math:`1/2`
:math:`\mathrm{L}`         :math:`1/2`                     :math:`1/2`                     :math:`1/2`
:math:`\mathrm{M}`         :math:`1/2`                     :math:`0`                       :math:`1/2`
:math:`\mathrm{N}`         :math:`1/2`                     :math:`0`                       :math:`0`
:math:`\mathrm{N_1}`       :math:`0`                       :math:`-1/2`                    :math:`0`
:math:`\mathrm{X}`         :math:`1/2`                     :math:`-1/2`                    :math:`0`
:math:`\mathrm{Y}`         :math:`\mu`                     :math:`\mu`                     :math:`\delta`
:math:`\mathrm{Y_1}`       :math:`1-\mu`                   :math:`-\mu`                    :math:`-\delta`
:math:`\mathrm{Y_2}`       :math:`-\mu`                    :math:`-\mu`                    :math:`-\delta`
:math:`\mathrm{Y_3}`       :math:`\mu`                     :math:`\mu-1`                   :math:`\delta`
:math:`\mathrm{Z}`         :math:`0`                       :math:`0`                       :math:`1/2`
=========================  ==============================  ==============================  ==============================


Variations
==========

There are five variations for base-centered monoclinic lattice.

Reciprocal :math:`\gamma` (:math:`k_{\gamma}`) is defined by the equation (for primitive lattice):

.. math::

    \cos(k_{\gamma}) = \frac{a^2 - b^2\sin^2(\alpha)}{a^2 + b^2\sin^2(\alpha)}

For MCLC\ :sub:`2` :math:`k_{\gamma} = 90`, therefore :math:`a = b \sin(\alpha)`.
For MCLC\ :sub:`1` we choose :math:`a < b \sin(\alpha)` and
for MCLC\ :sub:`3`, MCLC\ :sub:`4` and MCLC\ :sub:`5` we choose :math:`a > b \sin(\alpha)`.

For the variations 3-5 we define :math:`a = xb\sin(\alpha)`, where :math:`x > 1`.

Then the condition for MCLC\ :sub:`4` gives:

.. math::

    c = \frac{x^2}{x^2 - 1}b\cos(\alpha)

Where :math:`\cos(\alpha) > 0` (:math:`\alpha < 90^{\circ}`), since :math:`x > 1`.

And the ordering condition :math:`b \le c` gives:

.. math::

    \cos(\alpha) \ge \frac{x^2 - 1}{x^2}

For MCLC\ :sub:`3` (MCLC\ :sub:`5`) we choose parameters in a same way as for MCLC\ :sub:`4`,
but with :math:`c > \frac{x^2}{x^2 - 1}b\cos(\alpha)` (:math:`c < \frac{x^2}{x^2 - 1}b\cos(\alpha)`)

MCLC\ :sub:`1`
--------------

:math:`k_{\gamma} > 90^{\circ}`,

Predefined example: ``mclc1`` with :math:`a = \pi`, :math:`b = 1.4\cdot\pi`, :math:`c = 1.7\cdot\pi` and :math:`\alpha = 80^{\circ}`

MCLC\ :sub:`2`
--------------

:math:`k_{\gamma} = 90^{\circ}`,

Predefined example: ``mclc2`` with :math:`a = 1.4\cdot\pi\cdot\sin(75^{\circ})`, :math:`b = 1.4\cdot\pi`, :math:`c = 1.7\cdot\pi` and :math:`\alpha=75^{\circ}`

MCLC\ :sub:`3`
--------------

:math:`k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} < 1`

Predefined example with :math:`b = \pi`, :math:`x = 1.1`, :math:`\alpha = 78^{\circ}`, which produce:

``mclc4`` with :math:`a = 1.1\cdot\sin(78)\cdot\pi`, :math:`b = \pi`,
:math:`c = 1.8\cdot 121\cdot\cos(65)\cdot\pi/21` and :math:`\alpha = 78^{\circ}`

MCLC\ :sub:`4`
--------------

:math:`k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} = 1`

Predefined example with :math:`b = \pi`, :math:`x = 1.2`, :math:`\alpha = 65^{\circ}`, which produce:

``mclc4`` with :math:`a = 1.2\sin(65)\pi`, :math:`b = \pi`,
:math:`c = 36\cos(65)\pi/11` and :math:`\alpha = 65^{\circ}`

MCLC\ :sub:`5`
--------------

:math:`k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} > 1`

Predefined example with :math:`b = \pi`, :math:`x = 1.4`, :math:`\alpha = 53^{\circ}`, which produce:

``mclc5`` with :math:`a = 1.4\cdot\sin(53)\cdot\pi`, :math:`b = \pi`,
:math:`c = 0.9\cdot 11\cdot\cos(53)\cdot\pi/6` and :math:`\alpha = 53^{\circ}`


Examples
========

MCLC\ :sub:`1`
--------------

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: mclc1_brillouin.py
    :language: py

.. raw:: html
    :file: mclc1_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: mclc1_real.py
    :language: py

.. raw:: html
    :file: mclc1_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: mclc1_wigner-seitz.py
    :language: py

.. raw:: html
    :file: mclc1_wigner-seitz.html

MCLC\ :sub:`2`
--------------

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: mclc2_brillouin.py
    :language: py

.. raw:: html
    :file: mclc2_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: mclc2_real.py
    :language: py

.. raw:: html
    :file: mclc2_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: mclc2_wigner-seitz.py
    :language: py

.. raw:: html
    :file: mclc2_wigner-seitz.html

MCLC\ :sub:`3`
--------------

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: mclc3_brillouin.py
    :language: py

.. raw:: html
    :file: mclc3_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: mclc3_real.py
    :language: py

.. raw:: html
    :file: mclc3_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: mclc3_wigner-seitz.py
    :language: py

.. raw:: html
    :file: mclc3_wigner-seitz.html

MCLC\ :sub:`4`
--------------

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: mclc4_brillouin.py
    :language: py

.. raw:: html
    :file: mclc4_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: mclc4_real.py
    :language: py

.. raw:: html
    :file: mclc4_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: mclc4_wigner-seitz.py
    :language: py

.. raw:: html
    :file: mclc4_wigner-seitz.html

MCLC\ :sub:`5`
--------------

Brillouin zone and default kpath
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: mclc5_brillouin.py
    :language: py

.. raw:: html
    :file: mclc5_brillouin.html

Primitive and conventional cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: mclc5_real.py
    :language: py

.. raw:: html
    :file: mclc5_real.html

Wigner-Seitz cell
^^^^^^^^^^^^^^^^^
.. literalinclude:: mclc5_wigner-seitz.py
    :language: py

.. raw:: html
    :file: mclc5_wigner-seitz.html

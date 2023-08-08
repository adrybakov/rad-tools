.. _guide_orci:

********************************
Body-centred orthorhombic (ORCI)
********************************

**Pearson symbol**: oI

**Constructor**:  :py:func:`.ORCI`.


It is defined by three parameter: :math:`a`, :math:`b` and :math:`c` 
with conventional lattice:

.. math::

    \boldsymbol{a}_1 = (a, 0, 0)

    \boldsymbol{a}_2 = (0, b, 0)

    \boldsymbol{a}_3 = (0, 0, c)

And primitive lattice:

.. math::

    \boldsymbol{a}_1 = (-a/2, b/2, c/2)

    \boldsymbol{a}_2 = (a/2, -b/2, c/2)

    \boldsymbol{a}_3 = (a/2, b/2, -c/2)

Order of parameters: :math:`a < b < c`

Variations
==========

There are no variations for body-centered orthorombic. 
One example is predefined: ``orci`` with 
:math:`a = \pi`, :math:`b  = 1.3\pi` and :math:`c = 1.7\pi`.

Example structure
=================

**Default kpath**: :math:`\mathrm{\Gamma-X-L-T-W-R-X_1-Z-\Gamma-Y-S-W\vert L_1-Y\vert Y_1-Z}`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orci_brillouin.png 
            :target: ../../../../../_images/orci_brillouin.png 
      - .. literalinclude:: orci_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orci_real.png 
            :target: ../../../../../_images/orci_real.png 
      - .. literalinclude:: orci_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: orci_wigner-seitz.png 
            :target: ../../../../../_images/orci_wigner-seitz.png 
      - .. literalinclude:: orci_wigner-seitz.py
            :language: py


Ordering of lattice parameters
==============================
TODO

Edge cases
==========
If :math:`a = b \ne c` or :math:`a = c \ne b` or :math:`b = c \ne a`, 
then the lattice is :ref:`guide_bct`.

If :math:`a = b = c`, then the lattice is :ref:`guide_bcc`.

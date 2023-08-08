.. _guide_fcc:

************************
Face-centred cubic (FCC)
************************

**Pearson symbol**: cF

**Constructor**:  :py:func:`.FCC`.

It is defined by one parameter: :math:`a` with conventional lattice:

.. math::

    \boldsymbol{a}_1 = (a, 0, 0)

    \boldsymbol{a}_2 = (0, a, 0)

    \boldsymbol{a}_3 = (0, 0, a)

And primitive lattice:

.. math::

    \boldsymbol{a}_1 = (0, a/2, a/2)

    \boldsymbol{a}_2 = (a/2, 0, a/2)

    \boldsymbol{a}_3 = (a/2, a/2, 0)

Variations
==========

There are no variations for face-centered cubic lattice. 
One example is predefined: ``fcc`` with :math:`a = \pi`.

Example structure
=================

**Default kpath**: :math:`\mathrm{\Gamma-X-W-K-\Gamma-L-U-W-L-K\vert U-X}`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: fcc_brillouin.png 
            :target: ../../../../../_images/fcc_brillouin.png 
      - .. literalinclude:: fcc_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: fcc_real.png 
            :target: ../../../../../_images/fcc_real.png 
      - .. literalinclude:: fcc_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: fcc_wigner-seitz.png 
            :target: ../../../../../_images/fcc_wigner-seitz.png 
      - .. literalinclude:: fcc_wigner-seitz.py
            :language: py

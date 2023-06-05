.. _lattice-bct:

*****************************
Body-centred tetragonal (BCT)
*****************************

**Pearson symbol**: tI

Body-centered tetragonal lattice is described by the class :py:class:`.BCT`.

It is defined by two parameters: :math:`a` and :math:`c` 
with conventional lattice:

.. math::

    \boldsymbol{a}_1 = (a, 0, 0)

    \boldsymbol{a}_2 = (0, a, 0)

    \boldsymbol{a}_3 = (0, 0, c)

And primitive lattice:

.. math::

    \boldsymbol{a}_1 = (-a/2, a/2, c/2)

    \boldsymbol{a}_2 = (a/2, -a/2, c/2)

    \boldsymbol{a}_3 = (a/2, a/2, -c/2)

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


Example structures
==================

BCT\ :sub:`1`
-------------

**Default kpath**: :math:`\Gamma-X-M-\Gamma-Z-P-N-Z_1-M\vert X-P`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: bct1_brillouin.png 
            :target: ../../../../../_images/bct1_brillouin.png 
      - .. literalinclude:: bct1_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: bct1_real.png 
            :target: ../../../../../_images/bct1_real.png 
      - .. literalinclude:: bct1_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: bct1_wigner-seitz.png 
            :target: ../../../../../_images/bct1_wigner-seitz.png 
      - .. literalinclude:: bct1_wigner-seitz.py
            :language: py

BCT\ :sub:`2`
-------------

**Default kpath**: :math:`\Gamma-X-Y-\Sigma-\Gamma-Z-\Sigma_1-N-P-Y_1-Z\vert X-P`.

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: bct2_brillouin.png 
            :target: ../../../../../_images/bct2_brillouin.png 
      - .. literalinclude:: bct2_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: bct2_real.png 
            :target: ../../../../../_images/bct2_real.png 
      - .. literalinclude:: bct2_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: bct2_wigner-seitz.png 
            :target: ../../../../../_images/bct2_wigner-seitz.png 
      - .. literalinclude:: bct2_wigner-seitz.py
            :language: py

Edge cases
==========

If :math:`a = c` then the lattice is :ref:`lattice-bcc`.
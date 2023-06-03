.. _lattice-cub:

***
CUB
***

Cubic lattice is described by the class :py:class:`.CUB`.

It is defined by one parameter: :math:`a` with primitive and conventional lattice:

.. math::

    \boldsymbol{a}_1 = (a, 0, 0)

    \boldsymbol{a}_2 = (0, a, 0)

    \boldsymbol{a}_3 = (0, 0, a)

.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: cub_brillouin.png 
            :target: ../../../../../_images/cub_brillouin.png 
      - .. literalinclude:: plot_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: cub_real.png 
            :target: ../../../../../_images/cub_real.png 
      - .. literalinclude:: plot_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: cub_wigner-seitz.png 
            :target: ../../../../../_images/cub_wigner-seitz.png 
      - .. literalinclude:: plot_wigner-seitz.py
            :language: py
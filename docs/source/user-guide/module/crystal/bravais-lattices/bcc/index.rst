.. _lattice-bcc:

***
BCC
***

Body-centered cubic lattice is described by the class :py:class:`.BCC`.

It is defined by one parameter: :math:`a` with conventional lattice:

.. math::

    \boldsymbol{a}_1 = (a, 0, 0)

    \boldsymbol{a}_2 = (0, a, 0)

    \boldsymbol{a}_3 = (0, 0, a)

And primitive lattice:

.. math::

    \boldsymbol{a}_1 = (-a/2, a/2, a/2)

    \boldsymbol{a}_2 = (a/2, -a/2, a/2)

    \boldsymbol{a}_3 = (a/2, a/2, -a/2)


.. list-table:: Brillouin zone and default kpath
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: bcc_brillouin.png 
            :target: ../../../../../_images/bcc_brillouin.png 
      - .. literalinclude:: plot_brillouin.py
            :language: py

.. list-table:: Primitive and conventional cell
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: bcc_real.png 
            :target: ../../../../../_images/bcc_real.png 
      - .. literalinclude:: plot_real.py
            :language: py

.. list-table:: Wigner-Seitz cell
    :widths: 70 30
    :header-rows: 1

    * - Picture
      - Code
    * - .. figure:: bcc_wigner-seitz.png 
            :target: ../../../../../_images/bcc_wigner-seitz.png 
      - .. literalinclude:: plot_wigner-seitz.py
            :language: py

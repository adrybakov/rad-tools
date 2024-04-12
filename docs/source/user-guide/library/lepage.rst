.. _library_lepage:

****************
LePage algorithm
****************

The algorithm for the Bravais lattice type identification is bases on the search for the
twofold rotation axes and described in details in [1]_.
Here we recall the algorithm from Table 1 of the original publication with minor modifications as
implemented in :py:func:`.lepage` function.

3D Lattice
==========

Relative coordinates are given with respect to the primitive direct unit cell
as defined in [2]_ unless specified directly.


Step I
------
Compute niggli reduced cell via the call to the :py:func:`.niggli` function.
Define cell and reciprocal cell as Niggli cell and it's reciprocal pair.

Step II
--------
Generate all (125) the Miller indices for direct (:math:`U`) and
reciprocal (:math:`h`) lattice:

.. code-block:: python

    [-2,-2,-2],
    ...,
    [2, 2, 2]

Step III
--------
Run over all miller indices of direct cell and of all miller indices
of reciprocal cell.

If :math:`\vert U \cdot h\vert = 2` or :math:`\vert U \cdot h\vert = 1` compute :math:`\delta`:

.. math::

    \delta = \frac{\vert \boldsymbol{t}\times\boldsymbol{\tau}\vert}{\vert \boldsymbol{t}\cdot\boldsymbol{\tau}\vert}

where (@ - means matrix multiplication)

.. code-block:: python

    t = U @ cell
    tau = h @ reciprocal_cell

If :math:`\delta` is less then a ``limit``, keep the entry:

.. math::

    (U, \frac{\boldsymbol{t}}{\vert\boldsymbol{t}\vert}, \vert U \cdot h\vert, \delta)

Step IV
-------
Filter results to eliminate "twins" of direct indices:

* :math:`[1, 2, 0]` and :math:`[1, 2, 0]`
* :math:`[1, 2, 0]` and :math:`[-1, -2, 0]`
* :math:`[1, 0, -1]` and :math:`[2, 0, -2]`: :math:`[1, 0, -1]` is kept.

At this step the ``axes`` list is computed

Step V
------

Compute pairwise cosines between the twofold axes directions (``angles``):

.. math::

    \cos(\alpha_{ij}) = \frac{\vert\boldsymbol{t}_i\cdot\boldsymbol{t}_j\vert}{\vert\boldsymbol{t}_i\vert\vert\boldsymbol{t}_j\vert}

.. _step-vii:

Step VI
-------
If you came from go to: delete axes with the biggest delta from ``axes`` list.

Compute maximum :math:`\delta`. If it is less then maximum delta, finish and return result.
If not proceed with the consecutive checks for the system types.

Step VII (cubic)
----------------

Check for the cubic system. Cubic system has
9 even-order symmetry axis (in relative coordinates):

.. math::

    \begin{matrix}
        1:& (1, 0, 0) \\
        2:& (0, 1, 0) \\
        3:& (0, 0, 1) \\
        4:& (1, 1, 0) \\
        5:& (1, -1, 0) \\
        6:& (1, 0, 1) \\
        7:& (-1, 0, 1) \\
        8:& (0, 1, 1) \\
        9:& (0, -1, 1) \\
    \end{matrix}

with the following angle matrix:

.. math::

    \begin{matrix}
          & 1          & 2          & 3          & 4          & 5          & 6          & 7          & 8          & 9          \\
        1 & 0^{\circ}  & 90^{\circ} & 90^{\circ} & 45^{\circ} & 45^{\circ} & 45^{\circ} & 45^{\circ} & 90^{\circ} & 90^{\circ} \\
        2 & 90^{\circ} & 0^{\circ}  & 90^{\circ} & 45^{\circ} & 45^{\circ} & 90^{\circ} & 90^{\circ} & 45^{\circ} & 45^{\circ} \\
        3 & 90^{\circ} & 90^{\circ} & 0^{\circ}  & 90^{\circ} & 90^{\circ} & 45^{\circ} & 45^{\circ} & 45^{\circ} & 45^{\circ} \\
        4 & 45^{\circ} & 45^{\circ} & 90^{\circ} & 0^{\circ}  & 90^{\circ} & 60^{\circ} & 60^{\circ} & 60^{\circ} & 60^{\circ} \\
        5 & 45^{\circ} & 45^{\circ} & 90^{\circ} & 90^{\circ} & 0^{\circ}  & 60^{\circ} & 60^{\circ} & 60^{\circ} & 60^{\circ} \\
        6 & 45^{\circ} & 90^{\circ} & 45^{\circ} & 60^{\circ} & 60^{\circ} & 0^{\circ}  & 90^{\circ} & 60^{\circ} & 60^{\circ} \\
        7 & 45^{\circ} & 90^{\circ} & 45^{\circ} & 60^{\circ} & 60^{\circ} & 90^{\circ} & 0^{\circ}  & 60^{\circ} & 60^{\circ} \\
        8 & 90^{\circ} & 45^{\circ} & 45^{\circ} & 60^{\circ} & 60^{\circ} & 60^{\circ} & 60^{\circ} & 0^{\circ}  & 90^{\circ} \\
        9 & 90^{\circ} & 45^{\circ} & 45^{\circ} & 60^{\circ} & 60^{\circ} & 60^{\circ} & 60^{\circ} & 60^{\circ} & 0^{\circ}
    \end{matrix}

If ``angles`` is the same as the cubic angle matrix,
then find three axes with the following set of angles:
:math:`(0^{\circ} \times 1, 90^{\circ}\times 4, 45^{\circ} \times 4)`, put their Miller indices
in the matrix and compute its determinant :math:`\Delta`.

* If :math:`\vert\Delta\vert = 1`, then set system type to :ref:`"CUB" <guide_cub>`.
* If :math:`\vert\Delta\vert = 2`, then set system type to :ref:`"BCC" <guide_bcc>`.
* If :math:`\vert\Delta\vert = 4`, then set system type to :ref:`"FCC" <guide_fcc>`.

Go to :ref:`step-vii`.

Step VIII (hexagonal)
---------------------

Check for the hexagonal system. Hexagonal system has
7 even-order symmetry axis (in relative coordinates):

.. math::

    \begin{matrix}
        1:& (1, 0, 0) \\
        2:& (2, 1, 0) \\
        3:& (1, 1, 0) \\
        4:& (1, 2, 0) \\
        5:& (0, 1, 0) \\
        6:& (-1, 1, 0) \\
        7:& (0, 0, 1)
    \end{matrix}

with the following angle matrix:

.. math::

    \begin{matrix}
          & 1          & 2          & 3          & 4          & 5          & 6          & 7          \\
        1 & 0^{\circ}  & 30^{\circ} & 60^{\circ} & 90^{\circ} & 60^{\circ} & 30^{\circ} & 90^{\circ} \\
        2 & 30^{\circ} & 0^{\circ}  & 30^{\circ} & 60^{\circ} & 90^{\circ} & 60^{\circ} & 30^{\circ} \\
        3 & 60^{\circ} & 30^{\circ} & 0^{\circ}  & 30^{\circ} & 60^{\circ} & 90^{\circ} & 90^{\circ} \\
        4 & 90^{\circ} & 60^{\circ} & 30^{\circ} & 0^{\circ}  & 30^{\circ} & 60^{\circ} & 90^{\circ} \\
        5 & 60^{\circ} & 90^{\circ} & 60^{\circ} & 30^{\circ} & 0^{\circ}  & 30^{\circ} & 90^{\circ} \\
        6 & 30^{\circ} & 60^{\circ} & 90^{\circ} & 60^{\circ} & 30^{\circ} & 0^{\circ}  & 90^{\circ} \\
        7 & 90^{\circ} & 90^{\circ} & 90^{\circ} & 90^{\circ} & 90^{\circ} & 90^{\circ} & 0^{\circ}
    \end{matrix}

If ``angles`` is the same as the hexagonal angle matrix,
then set system type to :ref:`"HEX" <guide_hex>`.

Go to :ref:`step-vii`.

Step IX (tetragonal)
--------------------

Check for the tetragonal system. Tetragonal system has
5 even-order symmetry axis (in relative coordinates):

.. math::

    \begin{matrix}
        1:& (1, 0, 0) \\
        2:& (0, 1, 0) \\
        3:& (0, 0, 1) \\
        4:& (1, 1, 0) \\
        5:& (1, -1, 0)
    \end{matrix}

with the following angle matrix:

.. math::

    \begin{matrix}
          & 1          & 2          & 3          & 4          & 5          \\
        1 & 0^{\circ}  & 90^{\circ} & 90^{\circ} & 45^{\circ} & 45^{\circ} \\
        2 & 90^{\circ} & 0^{\circ}  & 90^{\circ} & 45^{\circ} & 45^{\circ} \\
        3 & 90^{\circ} & 90^{\circ} & 0^{\circ}  & 90^{\circ} & 90^{\circ} \\
        4 & 45^{\circ} & 45^{\circ} & 90^{\circ} & 0^{\circ}  & 90^{\circ} \\
        5 & 45^{\circ} & 45^{\circ} & 90^{\circ} & 90^{\circ} & 0^{\circ}
    \end{matrix}

If ``angles`` is the same as the tetragonal angle matrix,
then find one axes with the following set of angles:
:math:`(0^{\circ} \times 1, 90^{\circ}\times 4)`. Take two axes with minimal length form the remaining four.
Make a matrix from the Miller indices of the three axes
and compute its determinant :math:`\Delta`.

* If :math:`\vert\Delta\vert  = 1`, then set system type to :ref:`"TET" <guide_tet>`.
* If :math:`\vert\Delta\vert  = 2`, then set system type to :ref:`"BCT" <guide_bct>`.

Go to :ref:`step-vii`.

Step X (rhombohedral)
---------------------

Check for the rhombohedral system. Rhombohedral system has
3 even-order symmetry axis (in relative coordinates):

.. math::

    \begin{matrix}
        1:& (1, -1, 0) \\
        2:& (0, 1, -1) \\
        3:& (1, 0, -1) \\
    \end{matrix}

with the following angle matrix:

.. math::

    \begin{matrix}
          & 1          & 2          & 3           \\
        1 & 0^{\circ}  & 60^{\circ} & 60^{\circ} \\
        2 & 60^{\circ} & 0^{\circ}  & 60^{\circ} \\
        3 & 60^{\circ} & 60^{\circ} & 0^{\circ}
    \end{matrix}

If ``angles`` is the same as the rhombohedral angle matrix,
then set system type to :ref:`"RHL" <guide_rhl>`.

Go to :ref:`step-vii`.

Step XI (orthorhombic)
----------------------

Check for the orthorhombic system. Orthorhombic system has
3 even-order symmetry axis (in relative coordinates):

.. math::

    \begin{matrix}
        1:& (1, 0, 0) \\
        2:& (0, 1, 0) \\
        3:& (0, 0, 1) \\
    \end{matrix}

with the following angle matrix:

.. math::

    \begin{matrix}
          & 1          & 2          & 3           \\
        1 & 0^{\circ}  & 90^{\circ} & 90^{\circ}  \\
        2 & 90^{\circ} & 0^{\circ}  & 90^{\circ}  \\
        3 & 90^{\circ} & 90^{\circ} & 0
    \end{matrix}

If ``angles`` is the same as the orthorhombic angle matrix,
then make a matrix from the Miller indices of the three symmetry axes and
compute its determinant :math:`\Delta`.

* If :math:`\vert\Delta\vert  = 1`, then set system type to :ref:`"ORC" <guide_orc>`.
* If :math:`\vert\Delta\vert  = 4`, then set system type to ":ref:`"ORCF" <guide_orcf>`.
* If :math:`\vert\Delta\vert  = 2`, then check for :ref:`"ORCC" <guide_orcc>` vs :ref:`"ORCI" <guide_orci>`.
    Define matrix :math:`C` as the matrix where columns are the Miller indices of
    the three symmetry axes. Compute the vector:

    .. code-block:: python

        v = C @ [1, 1, 1]

    If the elements of v are |coprime|_, then set system type to :ref:`"ORCI" <guide_orci>`,
    otherwise set the system type to :ref:`"ORCC" <guide_orcc>`.

Go to :ref:`step-vii`.


Step XII (monoclinic)
---------------------

Check for the monoclinic system. Monoclinic system has
1 even-order symmetry axis (in relative coordinates
with respect to the conventional lattice as defined in [2]_):

.. math::

    \begin{matrix}
        1:& (1, 0, 0) \\
    \end{matrix}

with the following angle matrix:

.. math::

    \begin{matrix}
          & 1         \\
        1 & 0^{\circ}
    \end{matrix}

If ``angles`` is the same as the monoclinic angle matrix,
then define two shortest translation vectors in the plane
perpendicular to the twofold rotation axis. Put Miller indices of these
two vectors and of twofold axis in a matrix and compute its determinant :math:`\Delta`

* If :math:`\vert\Delta\vert  = 1`, then set system type to :ref:`"MCL" <guide_mcl>`.
* If :math:`\vert\Delta\vert  = 2`, then set system type to :ref:`"MCLC" <guide_mclc>`.


Go to :ref:`step-vii`.

Step XIII (trigonal)
--------------------

If all previous checks failed set system type to :ref:`"TRI" <guide_tri>` and go to :ref:`step-vii`.


References
==========

.. [1] Le Page, Y., 1982.
    The derivation of the axes of the conventional unit cell from
    the dimensions of the Buerger-reduced cell.
    Journal of Applied Crystallography, 15(3), pp.255-259.

.. [2] Setyawan, W. and Curtarolo, S., 2010.
    High-throughput electronic band structure calculations:
    Challenges and tools.
    Computational materials science, 49(2), pp.299-312.

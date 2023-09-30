.. _library_magnon-dispersion-method:

*****************
Magnon dispersion
*****************

Method is based on the [1]_. with matrix diagonalisation via [2]_.

We believe that a mistake was made in the derivation in original paper [1]_, 
which is discussed in much details here: :download:`magnon-dispersion.pdf <magnon-dispersion.pdf>`.

Magnon dispersion is computed via diagonalization of the matrix 
:math:`\boldsymbol{h}(\boldsymbol{k})`.
This matrix is defined in the original paper [1]_ and corrected in our 
:download:`comment <magnon-dispersion.pdf>`:

.. math::

    \boldsymbol{h}(\boldsymbol{k}) = 
    \begin{pmatrix}
    2\boldsymbol{A}(\boldsymbol{k}) - 2\boldsymbol{C} & 2\boldsymbol{B}(\boldsymbol{k}) \\
    2\boldsymbol{B}^{\dagger}(\boldsymbol{k}) & 2\overline{\boldsymbol{A}(-\boldsymbol{k})} - 2\boldsymbol{C} \\
    \end{pmatrix}

Diagonalization of this matrix imposes two conditions on it: Hermicity and
positive definiteness. If both conditions are satisfied, then the :py:class:`.MagnonDispersion`
will return a set of positive eigenfrequencies. If positive definiteness
is not satisfied, we define special strategies for the three cases:

* It is positive semi-definite, but not positive definite. 
    Following [1]_ we add small positive number (:math:`10^{-8}`) to the 
    diagonal of the matrix :math:`\boldsymbol{h}(\boldsymbol{k})` and then diagonalize it.

* If it is negative definite.
    We multiply the matrix by :math:`-1`, diagonalize it
    and multiply the result by :math:`-1`. In that way the set of negative 
    eigenfrequencies is returned, which are not correct magnon energies, but it can give you an idea of the 
    magnetic structure stability (try to apply spin spiral with the 
    :ref:`Q vector <rad-plot-tb2j-magnons_spiral-vector>`, which correspond to the minimum of energy). 

* If it is negative semi-definite, but not negative definite.
    We add small negative number (:math:`-10^{-8}`) to the diagonal of the matrix 
    :math:`\boldsymbol{h}(\boldsymbol{k})` and follow previous case.

In all other cases the :py:class:`.MagnonDispersion` will return ``0`` or ``None``.


References
==========

.. [1] Toth, S. and Lake, B., 2015. 
    Linear spin wave theory for single-Q incommensurate magnetic structures. 
    Journal of Physics: Condensed Matter, 27(16), p.166002.

.. [2] Colpa, J.H.P., 1978. 
    Diagonalization of the quadratic boson hamiltonian. 
    Physica A: Statistical Mechanics and its Applications, 93(3-4), pp.327-353.

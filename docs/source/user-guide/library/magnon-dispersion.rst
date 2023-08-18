.. _library_magnon-dispersion-method:

*****************
Magnon dispersion
*****************

Method is based on the [1]_. with matrix diagonalisation via [2]_.

We believe that in original paper [1]_ there is a mistake in the derivation, 
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

Diagonalization of this matrix impose two condition on it: Hermicity and
positive definiteness. If both conditions are satisfied, then the :py:class:`.MagnonDispersion`
will return set of positive eigenfrequencies. In three cases, when positive definiteness
is not satisfied, we define special strategies for the :py:class:`.MagnonDispersion`:

* It is positive semidefinite, but not positive definite. 
    In that case following [1]_ we add small positive number (:math:`10^{-8}`) to the 
    diagonal of the matrix :math:`\boldsymbol{h}(\boldsymbol{k})` and then diagonalize it.

* If it is negative definite.
    We multiply the matrix by :math:`-1`, diagonalize it
    and multiply the result by :math:`-1`. In that way the result will be a set of negative 
    eigenfrequencies, which are not correct magnon energies, but can give you an idea of the 
    magnetic structure stability. 

* If it is negative semidefinite, but not negative definite.
    In that case we add small negative number (:math:`-10^{-8}`) to the
    diagonal of the matrix :math:`\boldsymbol{h}(\boldsymbol{k})` and then follow previous case.

In all other cases the :py:class:`.MagnonDispersion` will return ``0`` or ``None``.


References
==========

.. [1] Toth, S. and Lake, B., 2015. 
    Linear spin wave theory for single-Q incommensurate magnetic structures. 
    Journal of Physics: Condensed Matter, 27(16), p.166002.

.. [2] Colpa, J.H.P., 1978. 
    Diagonalization of the quadratic boson hamiltonian. 
    Physica A: Statistical Mechanics and its Applications, 93(3-4), pp.327-353.
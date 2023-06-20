.. _rad-tools_parameter:

.. currentmodule:: radtools

******************
Exchange Parameter
******************

For the full reference see :ref:`api_parameter`

:py:class:`ExchangeParameter` is a wrap around :numpy:`ndarray` with predefined properties 
specific for the exchange parameter.

Isotropic, anisotropic and DMI parts of exchange matrix are defined as follows:

* :py:attr:`.iso`:

.. math::

    J_{iso} = \dfrac{\text{tr}(\boldsymbol{J})}{3}

* :py:attr:`.aniso`:

.. math::

    \mathbf{J}_{aniso} = \boldsymbol{J}_{symm} = \dfrac{\boldsymbol{J} + \boldsymbol{J}^T}{2} - \dfrac{1}{3}\text{tr}(\boldsymbol{J})\cdot\mathbf{I}

* :py:attr:`.dmi`:

.. math::

    \vec{\boldsymbol{D}} = (D_x, D_y, D_z)

.. math::

    \mathbf{J}_{asymm} = \dfrac{\boldsymbol{J} - \boldsymbol{J}^T}{2} = 
    \begin{bmatrix}
    0    & D_z  & -D_y \\
    -D_z & 0    & D_x  \\
    D_y  & -D_x & 0    \\
    \end{bmatrix}


where :math:`\boldsymbol{J}` is the full exchange matrix. Isotropic exchange is a scalar,
Anisotropic part is a traceless part of full matrix`s symmetric part 
and antisymmetric part of full matrix is described by Dzyaloshinskii-Moriya vector.


Creation
========

It is created as any instance of a python class. 
The constructor takes the following arguments:

* :py:attr:`.matrix` - the :math:`3 \times 3` matrix of the exchange parameter.
* :py:attr:`.iso` - the isotropic part of the exchange parameter.
* :py:attr:`.aniso` - the anisotropic part of the exchange parameter.
* :py:attr:`.dmi` - the Dzyaloshinskii-Moriya vector.

The matrix is the main argument and it is used to calculate the other parameters.
Other arguments are used to calculate the matrix only if it is not provided.

.. doctest::

    >>> import radtools as rad
    >>> J = rad.ExchangeParameter(iso=1)
    >>> J
    ExchangeParameter(array([[1., 0., 0.],
           [0., 1., 0.],
           [0., 0., 1.]]))  

For the exchange parameter convenient string representation is defined:

.. doctest::

    >>> print(J)
     1.0000 0.0000 0.0000
     0.0000 1.0000 0.0000
     0.0000 0.0000 1.0000

.. doctest::

    >>> J = rad.ExchangeParameter(dmi=(1,1,1))
    >>> print(J)
      0.0000  1.0000 -1.0000
     -1.0000  0.0000  1.0000
      1.0000 -1.0000  0.0000

Every property of the parameter can be accessed as an attribute. 
Below we list examples for each group of properties.

Full exchange matrix
====================
Only full matrix is stored, every other parameter is calculated from it.

.. doctest::

    >>> J = rad.ExchangeParameter(matrix=[[1,2,3],[4,5,6],[7,8,9]])
    >>> J.matrix
    array([[1., 2., 3.],
           [4., 5., 6.],
           [7., 8., 9.]])

Symmetric and asymmetric part can be accessed as:

.. doctest::

    >>> J = rad.ExchangeParameter(matrix=[[1,2,3],[4,5,6],[7,8,9]])
    >>> J.symm_matrix
    array([[1., 3., 5.],
           [3., 5., 7.],
           [5., 7., 9.]])
    >>> J.asymm_matrix
    array([[ 0., -1., -2.],
           [ 1.,  0., -1.],
           [ 2.,  1.,  0.]])

Istoropic exchange
==================

Isotropic exchange is a scalar and it is defined as a trace of the exchange matrix:

.. doctest::

    >>> J = rad.ExchangeParameter(matrix=[[1,2,3],[4,5,6],[7,8,9]])
    >>> J.iso
    5.0

For isotroic exchange matrix form is defined:

.. doctest::

    >>> J = rad.ExchangeParameter(iso=1)
    >>> J.matrix
    array([[1., 0., 0.],
           [0., 1., 0.],
           [0., 0., 1.]])

Anisotropic exchange
====================

.. note::
    
    Anisotropic exchange is a traceless part of symmetric part of exchange matrix.

.. doctest::

    >>> J = rad.ExchangeParameter(matrix=[[1,2,3],[4,5,6],[7,8,9]])
    >>> J.aniso
    array([[-4.,  3.,  5.],
           [ 3.,  0.,  7.],
           [ 5.,  7.,  4.]])

For anisotropic exchange usually reduced to the diagonal part, 
thus diagonal part and its matrix form is explicitly defined:

.. doctest::

    >>> J = rad.ExchangeParameter(matrix=[[1,2,3],[4,5,6],[7,8,9]])
    >>> J.aniso_diagonal
    array([-4.,  0.,  4.])
    >>> J.aniso_diagonal_matrix
    array([[-4.,  0.,  0.],
           [ 0.,  0.,  0.],
           [ 0.,  0.,  4.]])

Dzyaloshinskii-Moriya interaction
=================================

.. note::
    
    Dzyaloshinskii-Moriya interaction is an antisymmetric part of exchange matrix.

.. doctest::

    >>> J = rad.ExchangeParameter(matrix=[[1,2,3],[4,5,6],[7,8,9]])
    >>> J.dmi
    array([-1.,  2., -1.])

Matrix form of DMI is accessible as:

.. doctest::

    >>> J = rad.ExchangeParameter(matrix=[[1,2,3],[4,5,6],[7,8,9]])
    >>> J.dmi_matrix
    array([[ 0., -1., -2.],
           [ 1.,  0., -1.],
           [ 2.,  1.,  0.]])

Two useful values are defined for DMI: its module and relative strength to isotropic exchange:

.. math::

    \frac{\vert \vec{\boldsymbol{D}}\vert}{\vert J\vert}

.. doctest::

    >>> J = rad.ExchangeParameter(matrix=[[1,2,3],[4,5,6],[7,8,9]])
    >>> print(f"{J.dmi_module:.4f}")
    2.4495
    >>> print(f"{J.rel_dmi:.4f}")
    0.4899

Arithmetic operations
=====================

Arithmetic operations are defined for the exchange parameter. 
They act on the full exchange matrix and return :py:class:`.ExchangeParameter` instance.

.. doctest::

    >>> J1 = rad.ExchangeParameter(matrix=[[1,2,3],[4,5,6],[7,8,9]])
    >>> J2 = rad.ExchangeParameter(matrix=[[1,2,3],[4,5,6],[7,8,9]])
    >>> J1 + J2
    ExchangeParameter(array([[ 2.,  4.,  6.],
           [ 8., 10., 12.],
           [14., 16., 18.]]))
    >>> J1 - J2
    ExchangeParameter(array([[0., 0., 0.],
           [0., 0., 0.],
           [0., 0., 0.]]))
    >>> 2 * J2
    ExchangeParameter(array([[ 2.,  4.,  6.],
           [ 8., 10., 12.],
           [14., 16., 18.]]))
    >>> J2 / 2
    ExchangeParameter(array([[0.5, 1. , 1.5],
           [2. , 2.5, 3. ],
           [3.5, 4. , 4.5]]))
    >>> J1 // 3
    ExchangeParameter(array([[0., 0., 1.],
           [1., 1., 2.],
           [2., 2., 3.]]))
    >>> J1 % 3
    ExchangeParameter(array([[1., 2., 0.],
           [1., 2., 0.],
           [1., 2., 0.]]))
    >>> -J1
    ExchangeParameter(array([[-1., -2., -3.],
           [-4., -5., -6.],
           [-7., -8., -9.]]))
    >>> abs(-J1)
    ExchangeParameter(array([[1., 2., 3.],
           [4., 5., 6.],
           [7., 8., 9.]]))
    >>> J1 == J2
    True
    >>> J1 != J2
    False

Matrix multiplication works with :py:class:`.ExchangeParameter` instance as well, but 
it returns :numpy:`ndarray` instance:

.. doctest::

    >>> import numpy as np
    >>> J1 = rad.ExchangeParameter(matrix=[[1,2,3],[4,5,6],[7,8,9]])
    >>> a = np.eye(3)
    >>> J1 @ a
    array([[1., 2., 3.],
           [4., 5., 6.],
           [7., 8., 9.]])


|NumPy|_ interface
==================

|array_interface|_ is defined for the exchange parameter. 
It aims any numpy function onto the full exchange matrix.

Only transpose is defined directly in order to return :py:class:`.ExchangeParameter` instance:

.. doctest::

    >>> J = rad.ExchangeParameter(matrix=[[1,2,3],[4,5,6],[7,8,9]])
    >>> J.T
    ExchangeParameter(array([[1., 4., 7.],
           [2., 5., 8.],
           [3., 6., 9.]]))

Any other numpy function should work as expected, however it is not tested:

.. doctest::

    >>> J = rad.ExchangeParameter(matrix=[[1,2,3],[4,5,6],[7,8,9]])
    >>> np.sum(J)
    45.0
    >>> np.diag(J)
    array([1., 5., 9.])
    >>> np.iscomplex(J)
    array([[False, False, False],
           [False, False, False],
           [False, False, False]])















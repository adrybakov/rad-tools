.. _guide_spinham_parameter:

.. currentmodule:: radtools

*****************
ExchangeParameter
*****************

For the full reference see :ref:`api_parameter`

:py:class:`ExchangeParameter` is a wrap around :numpy:`ndarray` with a number
of predefined properties, which are specific for the exchange parameter.

The full exchange parameter :math:`\boldsymbol{J}` is a :math:`3\times3` matrix 
of real numbers. It is usually separated into three parts: 

* Isotropic (Heisenberg) exchange (:py:attr:`.iso`):

.. math::

    J_{iso} = \dfrac{\text{tr}(\boldsymbol{J})}{3}

* Symmetric anisotropic exchange (:py:attr:`.aniso`):

.. math::

    \mathbf{J}_{aniso} = \boldsymbol{J}_{symm} = \dfrac{\boldsymbol{J} + \boldsymbol{J}^T}{2} - \dfrac{1}{3}\text{tr}(\boldsymbol{J})\cdot\mathbf{I}

* Antisymmetric anisotropic (Dzyaloshinskii-Moriya) exchange (:py:attr:`.dmi`):

.. math::

    \vec{\boldsymbol{D}} = (D_x, D_y, D_z)

.. math::

    \mathbf{J}_{asymm} = \dfrac{\boldsymbol{J} - \boldsymbol{J}^T}{2} = 
    \begin{bmatrix}
    0    & D_z  & -D_y \\
    -D_z & 0    & D_x  \\
    D_y  & -D_x & 0    \\
    \end{bmatrix}


Isotropic exchange is a scalar, anisotropic part is a traceless symmetric part of the 
full matrix and antisymmetric part of full matrix can be described 
via Dzyaloshinskii-Moriya vector.

Import
======
The following import statements are equivalent:

.. doctest::

    >>> # Exact import
    >>> from radtools.spinham.parameter import ExchangeParameter
    >>> # Explicit import
    >>> from radtools.spinham import ExchangeParameter
    >>> # Recommended import
    >>> from radtools import ExchangeParameter

Creation
========

The constructor of :py:class:`.ExchangeParameter` takes the following arguments:

* :py:attr:`.matrix` - the :math:`3 \times 3` matrix of the exchange parameter.
* :py:attr:`.iso` - the isotropic part of the exchange parameter.
* :py:attr:`.aniso` - the anisotropic part of the exchange parameter.
* :py:attr:`.dmi` - the Dzyaloshinskii-Moriya vector.

If the ``matrix`` is provided, all other parameters are ignored and corresponding
values are calculated from the matrix. If ``matrix`` is not provided,
the other parameters are used to calculate the matrix.

.. doctest::

    >>> J = ExchangeParameter(iso=1)
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
    >>> J = ExchangeParameter(dmi=(1, 1, 1))
    >>> print(J)
      0.0000  1.0000 -1.0000
     -1.0000  0.0000  1.0000
      1.0000 -1.0000  0.0000

It supports string formatting. The format is passed to :py:func:`.print_2d_array` function internally:

.. doctest::

    >>> print(f"{J:.2f}")
      0.00  1.00 -1.00
     -1.00  0.00  1.00
      1.00 -1.00  0.00
    >>> print(f"{J:.2e}")
      0.00e+00  1.00e+00 -1.00e+00
     -1.00e+00  0.00e+00  1.00e+00
      1.00e+00 -1.00e+00  0.00e+00

The full matrix of the exchange parameter is stored as a 
:math:`3 \times 3` :numpy:`ndarray`. Every other parameter is calculated from it.
The examples below are grouped by the type of exchange interaction.

Full exchange matrix
====================

.. doctest::

    >>> J = ExchangeParameter(matrix=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> J.matrix
    array([[1., 2., 3.],
           [4., 5., 6.],
           [7., 8., 9.]])

Symmetric and asymmetric part can be accessed as:

.. doctest::

    >>> J.symm_matrix
    array([[1., 3., 5.],
           [3., 5., 7.],
           [5., 7., 9.]])
    >>> J.asymm_matrix
    array([[ 0., -1., -2.],
           [ 1.,  0., -1.],
           [ 2.,  1.,  0.]])

Isotropic exchange
==================

Isotropic exchange is a scalar and it is defined as a trace of the exchange matrix 
divided by 3:

.. doctest::

    >>> J.iso
    5.0

For isotropic exchange matrix form is defined:

.. doctest::

    >>> J.iso_matrix
    array([[5., 0., 0.],
           [0., 5., 0.],
           [0., 0., 5.]])

Anisotropic exchange
====================

.. note::
    
    Anisotropic exchange is a traceless part of symmetric part of full exchange matrix.

.. doctest::

    >>> J.aniso
    array([[-4.,  3.,  5.],
           [ 3.,  0.,  7.],
           [ 5.,  7.,  4.]])

Anisotropic exchange is often reduced to the diagonal part, 
thus diagonal part and its matrix form are explicitly defined:

.. doctest::

    >>> J.aniso_diagonal
    array([-4.,  0.,  4.])
    >>> J.aniso_diagonal_matrix
    array([[-4.,  0.,  0.],
           [ 0.,  0.,  0.],
           [ 0.,  0.,  4.]])

Dzyaloshinskii-Moriya interaction
=================================

.. note::
    
    Dzyaloshinskii-Moriya interaction is an antisymmetric part of full exchange matrix.

.. doctest::

    >>> J.dmi
    array([-1.,  2., -1.])

Matrix form of DMI is accessible via:

.. doctest::

    >>> J.dmi_matrix
    array([[ 0., -1., -2.],
           [ 1.,  0., -1.],
           [ 2.,  1.,  0.]])

Two useful values are defined for DMI: its module and relative strength to isotropic exchange:

.. math::

    \frac{\vert \vec{\boldsymbol{D}}\vert}{\vert J\vert}

.. doctest::

    >>> round(J.dmi_module, 4)
    2.4495
    >>> round(J.rel_dmi, 4)
    0.4899

Arithmetic operations
=====================

Arithmetic operations are defined for the exchange parameter. 
They act on the full exchange matrix and return :py:class:`.ExchangeParameter` instance.

.. doctest::

    >>> J1 = ExchangeParameter(matrix=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> J2 = ExchangeParameter(matrix=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
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

Matrix multiplication works with instances of :py:class:`.ExchangeParameter` as well, but 
it returns instance of :numpy:`ndarray`:

.. doctest::

    >>> import numpy as np
    >>> a = np.eye(3)
    >>> J @ a
    array([[1., 2., 3.],
           [4., 5., 6.],
           [7., 8., 9.]])


|NumPy|_ interface
==================

|array_interface|_ is defined for the exchange parameter. 
It aims any numpy function on the full exchange matrix.

Transpose is defined explicitly in order to return :py:class:`.ExchangeParameter` instance:

.. doctest::

    >>> J.T
    ExchangeParameter(array([[1., 4., 7.],
           [2., 5., 8.],
           [3., 6., 9.]]))

Any other numpy function should work as expected, however it is not tested:

.. doctest::

    >>> np.sum(J)
    45.0
    >>> np.diag(J)
    array([1., 5., 9.])
    >>> np.iscomplex(J)
    array([[False, False, False],
           [False, False, False],
           [False, False, False]])















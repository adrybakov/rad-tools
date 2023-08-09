.. _guide_geometry:

********
Geometry
********

For the full reference see :ref:`api_geometry`.

.. currentmodule:: radtools


Import
======

.. doctest::

    >>> # Exact import
    >>> from radtools.geometry import volume
    >>> # Recommended import
    >>> from radtools import volume

Volume
======

Volume can be computed from a several types of input data:

* Cell

:math:`3\times3` array, rows are the cell vectors, columns are the :math:`xyz` coordinates.

.. doctest::

    >>> cell = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    >>> volume(cell)
    1.0

* Three vectors

.. doctest::

    >>> a = [1, 0, 0]
    >>> b = [0, 1, 0]
    >>> c = [0, 0, 1]
    >>> volume(a, b, c)
    1.0

* Six cell parameters

.. doctest::

    >>> a = 1
    >>> b = 1
    >>> c = 1
    >>> alpha = 90
    >>> beta = 90
    >>> gamma = 90
    >>> volume(a, b, c, alpha, beta, gamma)
    1.0

Angle
=====

Computes the smallest angle between two vectors.

.. doctest::

    >>> from radtools import angle
    >>> a = [1, 0, 0]
    >>> b = [0, 1, 0]
    >>> # In degrees by default
    >>> angle(a, b)
    90.0
    >>> # In radians
    >>> round(angle(a, b, radians=True), 4)
    1.5708
    >>> # For the zero vector the angle is not defined
    >>> angle(a, [0, 0, 0])
    Traceback (most recent call last):
    ...
    ValueError: Angle is ill defined (zero vector).

Coordinates
===========

Computes relative coordinates of a vector with respect to some basis.
Basis is defined by three vectors, which are passed to the function in a form of a 
:math:`3\times3` array. Rows of the cell are the basis vectors, columns are the 
:math:`xyz` coordinates.

.. doctest::

    >>> from radtools import absolute_to_relative
    >>> basis = [[2, 0, 0], [0, 4, 0], [0, 0, 8]]
    >>> vector = [1, 1, 1]
    >>> absolute_to_relative(basis, vector)
    array([0.5  , 0.25 , 0.125])

Parallelepiped
==============

Check if the parallelepiped can be formed from six parameters.
:math:`a`, :math:`b`, :math:`c` are the lengths of the three vectors forming the
parallelepiped. :math:`\alpha`, :math:`\beta`, :math:`\gamma` are the angles between
the vectors.

.. doctest::

    >>> from radtools import parallelepiped_check
    >>> a = 1
    >>> b = 1
    >>> c = 1
    >>> alpha = 90
    >>> beta = 90
    >>> gamma = 90
    >>> parallelepiped_check(a, b, c, alpha, beta, gamma)
    True
    >>> parallelepiped_check(a, b, c, alpha, beta, 181)
    False

Orthonormal basis
=================

Spans an orthonormal basis from one vector. The input vector is normalized and
forms :math:`\vec{e}_3`. The orthonormal basis is computed by rotation of the 
standard basis

.. math::

    \begin{aligned}
        \vec{e}_1^{st} &= (1, 0, 0) \\
        \vec{e}_2^{st} &= (0, 1, 0) \\
        \vec{e}_3^{st} &= (0, 0, 1) \\
    \end{aligned}

with the rotation around the axis perpendicular to both :math:`\vec{e}_3^{st}` and
:math:`\vec{e}_3` by the angle :math:`\theta` between them.

.. doctest::

    >>> from radtools import span_orthonormal_set
    >>> vector = [1, 1, 1]
    >>> span_orthonormal_set(vector)
    array([[ 0.78867513, -0.21132487, -0.57735027],
           [-0.21132487,  0.78867513, -0.57735027],
           [ 0.57735027,  0.57735027,  0.57735027]])
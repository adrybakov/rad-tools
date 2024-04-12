.. _guide_numerical:

******************
Numerical accuracy
******************

For the full reference see :ref:`api_numerical`.

.. currentmodule:: radtools


Import
======

.. doctest::

    >>> # Exact import
    >>> from radtools.numerical import compare_numerically
    >>> # Recommended import
    >>> from radtools import compare_numerically

Comparison
==========

Compares two numbers within a given tolerance.

Tolerance for the comparison is given as:

.. math::

    \varepsilon = atol + rtol \cdot y

By default :math:`atol` and :math:`rtol` are equal to :py:const:`.crystal.REL_TOL`
and :py:const:`.crystal.ABS_TOL` respectively. You can pass exact values for eps.

The conditions are translated as (source [1]_):

===============  ==========================================================================
Condition        Numerical condition
===============  ==========================================================================
:math:`x < y`    :math:`x < y - \varepsilon`
:math:`x > y`    :math:`y < x - \varepsilon`
:math:`x \le y`  not :math:`(y < x - \varepsilon)`
:math:`x \ge y`  not :math:`(x < y - \varepsilon)`
:math:`x = y`    not :math:`(x < y - \varepsilon` or :math:`y < x - \varepsilon)`
:math:`x \ne y`  :math:`x < y - \varepsilon` or :math:`y < x - \varepsilon`
===============  ==========================================================================

.. doctest::

    >>> compare_numerically(1.0, "==", 1.0 + 1e-10)
    True
    >>> compare_numerically(1.0, "==", 1.0 + 2e-4)
    False
    >>> compare_numerically(1.0, "!=", 1.0 + 2e-4, eps=1e-3)
    False


References
==========
.. [1] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
    Numerically stable algorithms for the computation of reduced unit cells.
    Acta Crystallographica Section A: Foundations of Crystallography,
    60(1), pp.1-6.

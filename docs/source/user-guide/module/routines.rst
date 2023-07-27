.. _rad-tools_routines:

.. currentmodule:: radtools

*****************
Isolated routines
*****************

For the full reference see :ref:`api_utils`.

A number of routines are defined in this module. 
Majority of them are used internally and are not intended for direct use.

In this guide we describe only the routines, which are intended for direct use.

Data presentation
=================

:py:func:`.print_2d_array` - print 2D array in a nice table format.

It provides custom formatting, colour highlighting:

.. hint::

    For colour highlighting pass ``highlight=True`` to the function.
    For colour highlighting to work, the terminal must support |ANSI|_ escape sequences.

    It highlights positive values in red, negative values in blue, zero values in green.
    Complex and real parts of complex numbers are highlighted separately.


Real-valued array
-----------------

.. doctest::

    >>> import radtools as rad
    >>> array = [[1, 2], [3, 4], [5, 6]]
    >>> rad.print_2d_array(array)
    ┌──────┬──────┐
    │ 1.00 │ 2.00 │
    ├──────┼──────┤
    │ 3.00 │ 4.00 │
    ├──────┼──────┤
    │ 5.00 │ 6.00 │
    └──────┴──────┘
    >>> rad.print_2d_array([[0, 1., -0.],[1, 1, 1]])
    ┌───────┬───────┬───────┐
    │  0.00 │  1.00 │ -0.00 │
    ├───────┼───────┼───────┤
    │  1.00 │  1.00 │  1.00 │
    └───────┴───────┴───────┘

Custom formatting
-----------------

.. doctest::

    >>> import radtools as rad
    >>> rad.print_2d_array(array, fmt="10.2f")
    ┌────────────┬────────────┐
    │       1.00 │       2.00 │
    ├────────────┼────────────┤
    │       3.00 │       4.00 │
    ├────────────┼────────────┤
    │       5.00 │       6.00 │
    └────────────┴────────────┘
    >>> rad.print_2d_array(array, fmt=".2f")
    ┌──────┬──────┐
    │ 1.00 │ 2.00 │
    ├──────┼──────┤
    │ 3.00 │ 4.00 │
    ├──────┼──────┤
    │ 5.00 │ 6.00 │
    └──────┴──────┘

No borders
----------

.. doctest::

    >>> import radtools as rad
    >>> rad.print_2d_array(array, borders=False)
     1.00 2.00
     3.00 4.00
     5.00 6.00

Shift
-----

.. doctest::

    >>> import radtools as rad
    >>> rad.print_2d_array(array, shift=3)
       ┌──────┬──────┐
       │ 1.00 │ 2.00 │
       ├──────┼──────┤
       │ 3.00 │ 4.00 │
       ├──────┼──────┤
       │ 5.00 │ 6.00 │
       └──────┴──────┘

Scientific notation
-------------------

.. doctest::

    >>> import radtools as rad
    >>> array = [[1, 2], [3, 4], [52414345345, 6]]
    >>> rad.print_2d_array(array, fmt="10.2E")
    ┌────────────┬────────────┐
    │   1.00E+00 │   2.00E+00 │
    ├────────────┼────────────┤
    │   3.00E+00 │   4.00E+00 │
    ├────────────┼────────────┤
    │   5.24E+10 │   6.00E+00 │
    └────────────┴────────────┘

Complex-valued array
--------------------

.. doctest::

    >>> import radtools as rad
    >>> array = [[1, 2 + 1j], [3, 4], [52, 6]]
    >>> rad.print_2d_array(array)
    ┌───────┬────────────────┐
    │  1.00 │  2.00 + i1.00  │
    ├───────┼────────────────┤
    │  3.00 │  4.00          │
    ├───────┼────────────────┤
    │ 52.00 │  6.00          │
    └───────┴────────────────┘
    >>> rad.print_2d_array(array, fmt="4.2E")
    ┌──────────┬──────────────────────┐
    │ 1.00E+00 │ 2.00E+00 + i1.00E+00 │
    ├──────────┼──────────────────────┤
    │ 3.00E+00 │ 4.00E+00             │
    ├──────────┼──────────────────────┤
    │ 5.20E+01 │ 6.00E+00             │
    └──────────┴──────────────────────┘
    >>> array = [[1, 2 - 1j], [3, 4], [52, 6]]
    >>> rad.print_2d_array(array)
    ┌───────┬────────────────┐
    │  1.00 │  2.00 - i1.00  │
    ├───────┼────────────────┤
    │  3.00 │  4.00          │
    ├───────┼────────────────┤
    │ 52.00 │  6.00          │
    └───────┴────────────────┘

Complex-valued array with real part equal to zero
-------------------------------------------------

.. doctest::

    >>> import radtools as rad
    >>> array = [[1, 1j], [3, 4], [52, 6]]
    >>> rad.print_2d_array(array)
    ┌───────┬────────────────┐
    │  1.00 │       + i1.00  │
    ├───────┼────────────────┤
    │  3.00 │  4.00          │
    ├───────┼────────────────┤
    │ 52.00 │  6.00          │
    └───────┴────────────────┘

Empty arrays
------------

.. doctest::

    >>> import radtools as rad
    >>> rad.print_2d_array([])
    None
    >>> rad.print_2d_array([[]])
    None


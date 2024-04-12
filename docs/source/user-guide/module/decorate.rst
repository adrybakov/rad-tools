.. _guide_decorate:

******************
Decoration of data
******************

For the full reference see :ref:`api_decorate`.

.. currentmodule:: radtools

This module contains functions to decorate data for plotting or printing.
It does not fall into the scope of intended functionality of the package,
but is included for convenience.

Import
======

.. doctest::

    >>> # Exact import
    >>> from radtools.decorate.array import print_2d_array
    >>> # Explicit import
    >>> from radtools.decorate import print_2d_array
    >>> # Recommended import
    >>> from radtools import print_2d_array

2D arrays
=========

This function can take any numerical 2D or 1D array and print it in a nice
format. It is useful for debugging and printing of results.


It provides custom formatting, colour highlighting:

.. hint::

    For colour highlighting pass ``highlight=True`` to the function.
    For colour highlighting to work, the terminal must support |ANSI|_ escape sequences.

    It highlights positive values in red, negative values in blue, zero values in green.
    Complex and real parts of complex numbers are highlighted separately.


Real-valued array
-----------------

.. doctest::

    >>> array = [[1, 2], [3, 4], [5, 6]]
    >>> print_2d_array(array)
    ┌──────┬──────┐
    │ 1.00 │ 2.00 │
    ├──────┼──────┤
    │ 3.00 │ 4.00 │
    ├──────┼──────┤
    │ 5.00 │ 6.00 │
    └──────┴──────┘
    >>> print_2d_array([[0, 1., -0.],[1, 1, 1]])
    ┌──────┬──────┬──────┐
    │ 0.00 │ 1.00 │ 0.00 │
    ├──────┼──────┼──────┤
    │ 1.00 │ 1.00 │ 1.00 │
    └──────┴──────┴──────┘

Custom formatting
-----------------

.. doctest::

    >>> print_2d_array(array, fmt="10.2f")
    ┌────────────┬────────────┐
    │       1.00 │       2.00 │
    ├────────────┼────────────┤
    │       3.00 │       4.00 │
    ├────────────┼────────────┤
    │       5.00 │       6.00 │
    └────────────┴────────────┘
    >>> print_2d_array(array, fmt=".2f")
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

    >>> print_2d_array(array, borders=False)
     1.00 2.00
     3.00 4.00
     5.00 6.00

Shift
-----

.. doctest::

    >>> print_2d_array(array, shift=3)
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

    >>> array = [[1, 2], [3, 4], [52414345345, 6]]
    >>> print_2d_array(array, fmt="10.2E")
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

    >>> array = [[1, 2 + 1j], [3, 4], [52, 6]]
    >>> print_2d_array(array)
    ┌───────┬──────────────┐
    │  1.00 │ 2.00 + i1.00 │
    ├───────┼──────────────┤
    │  3.00 │ 4.00         │
    ├───────┼──────────────┤
    │ 52.00 │ 6.00         │
    └───────┴──────────────┘
    >>> print_2d_array(array, fmt="4.2E")
    ┌──────────┬──────────────────────┐
    │ 1.00E+00 │ 2.00E+00 + i1.00E+00 │
    ├──────────┼──────────────────────┤
    │ 3.00E+00 │ 4.00E+00             │
    ├──────────┼──────────────────────┤
    │ 5.20E+01 │ 6.00E+00             │
    └──────────┴──────────────────────┘
    >>> array = [[1, 2 - 1j], [3, 4], [52, 6]]
    >>> print_2d_array(array)
    ┌───────┬──────────────┐
    │  1.00 │ 2.00 - i1.00 │
    ├───────┼──────────────┤
    │  3.00 │ 4.00         │
    ├───────┼──────────────┤
    │ 52.00 │ 6.00         │
    └───────┴──────────────┘

Complex-valued array with real part equal to zero
-------------------------------------------------

.. doctest::

    >>> array = [[1, 1j], [3, 4], [52, 6]]
    >>> print_2d_array(array)
    ┌───────┬──────────────┐
    │  1.00 │      + i1.00 │
    ├───────┼──────────────┤
    │  3.00 │ 4.00         │
    ├───────┼──────────────┤
    │ 52.00 │ 6.00         │
    └───────┴──────────────┘

Empty cells
-----------

.. doctest::

    >>> array = [[1, 2], [3, 4], [5, 6]]
    >>> array[1][1] = None
    >>> print_2d_array(array)
    ┌──────┬──────┐
    │ 1.00 │ 2.00 │
    ├──────┼──────┤
    │ 3.00 │      │
    ├──────┼──────┤
    │ 5.00 │ 6.00 │
    └──────┴──────┘
    >>> array[1][0] = None
    >>> print_2d_array(array)
    ┌──────┬──────┐
    │ 1.00 │ 2.00 │
    ├──────┼──────┤
    │      │      │
    ├──────┼──────┤
    │ 5.00 │ 6.00 │
    └──────┴──────┘
    >>> array[0][1] = None
    >>> array[2][1] = None
    >>> print_2d_array(array)
    ┌──────┬──┐
    │ 1.00 │  │
    ├──────┼──┤
    │      │  │
    ├──────┼──┤
    │ 5.00 │  │
    └──────┴──┘

Empty arrays
------------

.. doctest::

    >>> print_2d_array([])
    None
    >>> print_2d_array([[]])
    None

Long numbers
------------

.. doctest::

    >>> array = [[1, 2], [3, 4], [52414345345, 6]]
    >>> print_2d_array(array)
    ┌────────────────┬──────┐
    │           1.00 │ 2.00 │
    ├────────────────┼──────┤
    │           3.00 │ 4.00 │
    ├────────────────┼──────┤
    │ 52414345345.00 │ 6.00 │
    └────────────────┴──────┘

Headers and footers
-------------------

.. doctest::

    >>> array = [[1, 2], [3, 4], [5, 6]]
    >>> print_2d_array(array, header_row=["a", "b"])
    ┌──────┬──────┐
    │    a │    b │
    ├──────┼──────┤
    │ 1.00 │ 2.00 │
    ├──────┼──────┤
    │ 3.00 │ 4.00 │
    ├──────┼──────┤
    │ 5.00 │ 6.00 │
    └──────┴──────┘
    >>> print_2d_array(array, footer_row=["c", "d"])
    ┌──────┬──────┐
    │ 1.00 │ 2.00 │
    ├──────┼──────┤
    │ 3.00 │ 4.00 │
    ├──────┼──────┤
    │ 5.00 │ 6.00 │
    ├──────┼──────┤
    │    c │    d │
    └──────┴──────┘
    >>> print_2d_array(array, header_column=["a", "b", "c"])
    ┌───┬──────┬──────┐
    │ a │ 1.00 │ 2.00 │
    ├───┼──────┼──────┤
    │ b │ 3.00 │ 4.00 │
    ├───┼──────┼──────┤
    │ c │ 5.00 │ 6.00 │
    └───┴──────┴──────┘
    >>> print_2d_array(array, footer_column=["a", "b", "c"])
    ┌──────┬──────┬───┐
    │ 1.00 │ 2.00 │ a │
    ├──────┼──────┼───┤
    │ 3.00 │ 4.00 │ b │
    ├──────┼──────┼───┤
    │ 5.00 │ 6.00 │ c │
    └──────┴──────┴───┘
    >>> print_2d_array(array, header_column=["a", "B", "c"], header_row = ["", "A", "B"])
    ┌───┬──────┬──────┐
    │   │    A │    B │
    ├───┼──────┼──────┤
    │ a │ 1.00 │ 2.00 │
    ├───┼──────┼──────┤
    │ B │ 3.00 │ 4.00 │
    ├───┼──────┼──────┤
    │ c │ 5.00 │ 6.00 │
    └───┴──────┴──────┘
    >>> print_2d_array(array, header_column=["a", "B", "c"], header_row = ["corner", "A", "B"])
    ┌────────┬──────┬──────┐
    │ corner │    A │    B │
    ├────────┼──────┼──────┤
    │      a │ 1.00 │ 2.00 │
    ├────────┼──────┼──────┤
    │      B │ 3.00 │ 4.00 │
    ├────────┼──────┼──────┤
    │      c │ 5.00 │ 6.00 │
    └────────┴──────┴──────┘
    >>> print_2d_array(array, header_column=["a", "B", "c"], header_row = ["", "A", "B",""], footer_row = ["","c","d",""], footer_column=["c","d","e"])
    ┌───┬──────┬──────┬───┐
    │   │    A │    B │   │
    ├───┼──────┼──────┼───┤
    │ a │ 1.00 │ 2.00 │ c │
    ├───┼──────┼──────┼───┤
    │ B │ 3.00 │ 4.00 │ d │
    ├───┼──────┼──────┼───┤
    │ c │ 5.00 │ 6.00 │ e │
    ├───┼──────┼──────┼───┤
    │   │    c │    d │   │
    └───┴──────┴──────┴───┘
    >>> array = [[1, 2], [3, 4], [5, 6+1j]]
    >>> print_2d_array(array, header_row=["a", "b"])
    ┌──────┬──────────────┐
    │    a │            b │
    ├──────┼──────────────┤
    │ 1.00 │ 2.00         │
    ├──────┼──────────────┤
    │ 3.00 │ 4.00         │
    ├──────┼──────────────┤
    │ 5.00 │ 6.00 + i1.00 │
    └──────┴──────────────┘
    >>> print_2d_array(array, header_row=["a", "b"], fmt="^.2f")
    ┌──────┬──────────────┐
    │  a   │      b       │
    ├──────┼──────────────┤
    │ 1.00 │ 2.00         │
    ├──────┼──────────────┤
    │ 3.00 │ 4.00         │
    ├──────┼──────────────┤
    │ 5.00 │ 6.00 + i1.00 │
    └──────┴──────────────┘


Axis lines
==========

Two shortcuts with the common linestyle are defined for convenience:

.. doctest::

    >>> from radtools import plot_hlines, plot_vlines # doctest: +SKIP
    >>> import matplotlib.pyplot as plt # doctest: +SKIP
    >>> fig, ax = plt.subplots() # doctest: +SKIP
    >>> plot_hlines(ax, [1, 2, 3]) # doctest: +SKIP
    >>> plot_vlines(ax, [1, 2, 3]) # doctest: +SKIP
    >>> # You can pass any keyword arguments to the underlying matplotlib function
    >>> plot_hlines(ax, [1, 2, 3], color="red", linewidth=2) # doctest: +SKIP

Colormap
========

Another shortcut is defined for the custom colormap:

.. doctest::

    >>> from radtools import custom_cmap # doctest: +SKIP
    >>> colormap = custom_cmap("#e0218a", "#2e1b25") # doctest: +SKIP

Logo and stats
==============

Three functions are defined for printing the logo and stats of the package:

Logo
----

.. doctest::

    >>> from radtools import logo
    >>> print(logo()) # doctest: +SKIP
    ██████╗  █████╗ ██████╗       ████████╗ █████╗  █████╗ ██╗      ██████╗
    ██╔══██╗██╔══██╗██╔══██╗      ╚══██╔══╝██╔══██╗██╔══██╗██║     ██╔════╝
    ██████╔╝███████║██║  ██║█████╗   ██║   ██║  ██║██║  ██║██║     ╚═█████╗
    ██╔══██╗██╔══██║██║  ██║╚════╝   ██║   ██║  ██║██║  ██║██║       ╚══██║
    ██║  ██║██║  ██║██████╔╝         ██║   ╚█████╔╝╚█████╔╝███████╗██████╔╝
    ╚═╝  ╚═╝╚═╝  ╚═╝╚═════╝          ╚═╝    ╚════╝  ╚════╝ ╚══════╝╚═════╝
                                                                ▄   ▄
                          Version: 0.8.4                        █▀█▀█
                   Documentation: rad-tools.org                 █▄█▄█
                  Release date: 20 September 2023                ███   ▄▄
        Git hash: 5b4dc9b04aecb7efebb9cd710c02e0ad7fb68e44       ████ █  █
                        Licence: GNU GPLv3                       ████    █
                                                                 ▀▀▀▀▀▀▀▀

The same behaviour can be achieved within the console (you may need to use ``python3``):

.. code-block::

    python -m radtools

One-line summary
----------------

.. doctest::

    >>> from radtools import stamp_line
    >>> print(stamp_line()) # doctest: +SKIP
    on 9 August 2023 at 16:50:17 by rad-tools 0.8.0
    >>> print(stamp_line(doclink=True)) # doctest: +SKIP
    on 9 August 2023 at 16:51:0 by rad-tools 0.8.0 Documentation: rad-tools.org

License
-------

.. doctest::

    >>> from radtools import license
    >>> print(license()) # doctest: +SKIP

This code prints the full text of the GNU GPLv3 license.

The same behaviour can be achieved within the console (you may need to use ``python3``):

.. code-block::

    python -m radtools --license

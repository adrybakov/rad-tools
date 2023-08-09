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
    ┌──────┬──────┬──────┐
    │ 0.00 │ 1.00 │ 0.00 │
    ├──────┼──────┼──────┤
    │ 1.00 │ 1.00 │ 1.00 │
    └──────┴──────┴──────┘

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
                          Version: 0.8.0                        █▀█▀█      
                    Release date: 2 August 2023                 █▄█▄█      
        Git hash: 9c9b087fa02be0cafdadaef6ec1c7926fe36e3d6       ███   ▄▄  
                   Documentation: rad-tools.org                  ████ █  █ 
                       Licence: MIT License                      ████    █ 
                                                                 ▀▀▀▀▀▀▀▀  


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
    MIT License

    Copyright (c) 2022-2023 Andrey Rybakov

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

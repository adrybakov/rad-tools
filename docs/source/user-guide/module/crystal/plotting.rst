.. _guide_crystal_plotting:

.. currentmodule:: radtools

********
Plotting
********

For the full reference see :ref:`api_crystal-plotting`

More examples of the plots with the code snippets can be found in the 
:ref:`library_bravais-lattices` guide.

Import 
======

    >>> # Exact import
    >>> from radtools.crystal.lattice import Lattice
    >>> # Explicit import
    >>> from radtools.crystal import Lattice
    >>> # Recommended import
    >>> from radtools import Lattice

For the examples in this page we need additional import and some predefined variables:

.. doctest::

    >>> from radtools import lattice_example

Creation
========


Plotting of the lattice
=======================

Lattice primitive, conventional unit cells as  well as the Wigner-Seitz cell 
and Brillouin zone can be plotted using :py:meth:`.Lattice.plot` method:

.. doctest::

    >>> lattice = lattice_example("BCT")
    >>> lattice.plot("brillouin")
    >>> lattice.plot("kpath")

To show the plot use :py:meth:`.Lattice.show` method:

.. doctest::

    >>> lattice.show() # doctest: +SKIP


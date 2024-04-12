.. _guide_crystal_plotting:

.. currentmodule:: radtools

********
Plotting
********

For the full reference see :ref:`api_crystal-plotting`

More examples of the plots with the code snippets can be found in the
:ref:`library_bravais-lattices` guide.

.. note::

    The plotting is written for the :py:class:`.Lattice` and consequently
    available for :py:class`.Crystal` and :py:class:`.SpinHamiltonian` classes.

Import
======

.. doctest::

    >>> # Exact import
    >>> from radtools.crystal.lattice_plotter import MatplotlibBackend, PlotlyBackend
    >>> # Explicit import
    >>> from radtools.crystal import MatplotlibBackend, PlotlyBackend
    >>> # Recommended import
    >>> from radtools import MatplotlibBackend, PlotlyBackend


For the examples in this page we need additional import and some predefined variables:

.. doctest::

    >>> from radtools import HEX
    >>> lattice = HEX(2, 1)



Creation
========

Each backend is created as any other object in Python:

.. doctest::

    >>> mb = MatplotlibBackend()
    >>> pb = PlotlyBackend()


Plotting with Plotly
====================

|plotly|_ is considered to be the main plotting backend.

The creation of the plotly backend can take one optional argument: ``fig``:

.. doctest::

    >>> from plotly import graph_objects as go
    >>> fig = go.Figure()
    >>> pb = PlotlyBackend(fig=fig)

The typical workflow consist in plotting one or several lattices and then
showing or saving the figure:

.. doctest::

    >>> pb.plot(lattice, kind="brillouin", label="Brillouin zone", color="red")
    >>> pb.plot(lattice, kind="kpath", label="K-path", color="black")
    >>> pb.save('hex_lattice.html')
    >>> pb.show() # doctest: +SKIP

Plotting with Matplotlib
========================

|matplotlib|_ is considered to be the secondary plotting backend.

The creation of the matplotlib backend can take one optional arguments: ``fig and ``ax``:

.. doctest::

    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure(figsize=(6, 6))
    >>> ax = fig.add_subplot(projection="3d")
    >>> mb = MatplotlibBackend(fig=fig, ax=ax)

The typical workflow consist in plotting one or several lattices and then
showing or saving the figure:

.. doctest::

    >>> mb.plot(lattice, kind="brillouin", label="Brillouin zone", color="red")
    >>> mb.plot(lattice, kind="kpath", label="K-path", color="black")
    >>> mb.save('hex_lattice.png')
    >>> mb.show() # doctest: +SKIP

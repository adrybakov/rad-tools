.. _module-guide:

************
Module guide
************

Module guide is directed at the users, who wants to use the package as python library.
The main objective of the module guide is to
familiarize you with the package and provide examples of usage.

If you want to use command-line interface check out the :ref:`scripts-guide`.

Data structure and import
=========================

RAD-tools is a collection of submodules, which combine the classes and functions
with the same topic. For the full reference of the submodules, please refer to the
:ref:`api`.

All public methods of each submodule are exposed to the main entry
point (``radtools``):

.. doctest::

    >>> from radtools import print_2d_array

Explicit imports are supported as well:

.. doctest::

    >>> from radtools.decorate import print_2d_array

As well as the exact imports:

.. doctest::

        >>> from radtools.decorate.array import print_2d_array


Submodules
==========

For the detailed guide on each submodule please refer to the pages below

.. toctree::
    :maxdepth: 1

    geometry
    numerical
    decorate
    score

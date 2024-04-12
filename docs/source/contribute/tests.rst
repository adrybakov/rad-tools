.. _contribute_tests:

*******
Testing
*******

In RAD-tools we rely on |pytest|_, |hypothesis|_ and |doctest|_ for testing.

Unit tests
==========

All unit tests are located in the "utests/" directory.
To run the tests, you can use the following command (provided that the |GNU-make|_
command is available)::

    make test

Alternatively, you can run::

    pytest -s

Structure of the "utest/" directory loosely follows the structure of the "src/radtools/"
directory.

For example if you've added a new function named "rotate()" to the
geometry moduyle from the "src/radtools/geometry.py" file,
then you should add a new test function named "test_rotate()" to the file
"utest/test_geometry.py".

Documentation tests
===================

Across the documentation there are many examples of how to use RAD-tools with code snippets.
These code snippets are tested using |doctest|_, which ensure that the documentation
correctly reflects an actual behavior of the code. In order to run doctests you need
to build the :ref:`documentation <contribute_documentation>` and then run the doctests
using the following command (provided that the |GNU-make|_ command is available)::

    make doctest

Alternatively, you can run::

    sphinx-build -b doctest "docs/source" "docs/build"

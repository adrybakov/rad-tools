.. _contribute_documentation:

*************
Documentation
*************

The documentation of RAD-tools is build with |sphinx|_.

.. hint::

  The best way to get a feeling about how the documentation of RAD-tools is structured is
  to read the source code in the "docs/source/" directory and compare it's content and
  structure with this webpage. If you have any doubts we encourage you to
  :ref:`contact us <support>`.

Building the documentation
==========================

To build documentation simply run (provided that |GNU-make|_ command is available)::

  make html

Alternatively, you can run the following command::

  sphinx-build -M html "docs/source" "docs/build"

Documentation structure
=======================

The documentation falls into two big parts:

* API ("docs/source/api/" directory)

  Semi-automatically generated documentation of the source code, it is mostly build
  based on the docstrings of the source code, using |sphinx-autodoc|_ and
  |sphinx-autosummary|_.

  It's content loosely follows the structure of the source code folders. The methods and
  classes are listed manually for the better presentation of the documentation, the rest
  is automatically generated. Please consult existing files to get a feeling about the
  structure.

* User guide ("docs/source/user-guide/" directory)

  Hand-written |reStructuredText|_ files with the usage examples and the explanation of
  the RAD-tools's functionality.

  We separate the user guide into several parts:

  - "docs/source/user-guide/module/" directory
    Description of the python module. The majority of the doctests are written here.

  - "docs/source/user-guide/scripts/" directory
    Description of the scripts that are included in the RAD-tools.

  At the moment we have a few rules formulated for the documentation of the python module:

  - Each function and class has to have "Import" section.
  - Each class has to have "Creation" section.
  - Each page of the module guide have to have the link to the corresponding page of the
    API.

The rest of the documentation is located in the "docs/source/" directory and it includes,
among other things:

* "docs/source/conf.py" file

  The configuration file for |sphinx|_.

* "docs/source/index.rst" file

  The main page of the documentation. It includes the table of contents and the
  introduction to the RAD-tools.

* "docs/source/support.rst" file

  The page with the information about how to get support for the users of RAD-tools.

* "docs/source/release-notes/" directory

  The release notes for each version of RAD-tools.

* "docs/source/contribute/" directory

  Root folder for the documentation of how to contribute to RAD-tools.


Docstrings
==========

All public classes and functions have to have a docstring.
The docstring has to be written in the |numpydoc|_ style guide.

To get a feeling about the style you can read examples in the source code of RAD-tools.

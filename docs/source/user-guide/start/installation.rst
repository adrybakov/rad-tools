.. _rad-tools_installation:

************
Installation
************

Requirement for RAD-tools installation are:

* |Python|_ itself (>=3.6)
* |NumPy|_
* |SciPy|_
* |matplotlib|_
* |tqdm|_
* |termcolor|_

Most likely you already have Python installed on your machine
(if not check these links: |Python-installation|).

RAD-tools can be installed with :ref:`pip <installation-pip>` 
or from :ref:`source <installation-source>`.

Check Python
============

The easiest way to check if you have python installed 
is to run the following command in your terminal:

.. code-block:: console

   python

If you see something like that, then you have python installed:

.. code-block:: console

   Python 3.10.9 (main, Dec 15 2022, 18:25:35) [Clang 14.0.0 (clang-1400.0.29.202)] on darwin
   Type "help", "copyright", "credits" or "license" for more information.
   >>> 

In most cases ``python`` command launches python3, 
however if it launches python2, 
then you may need to use ``python3`` instead 
(and ``pip3`` instead of ``pip`` in the following).

.. hint::
   Use ``exit()`` to close python console.

The packages (|NumPy|_, |SciPy|_, |matplotlib|_, |tqdm|_, |termcolor|_) 
are installed automatically during the installation of RAD-tools, 
so you do not have to worry about them.

.. _installation-pip:

pip
===

To install RAD-tools, run (you may need to use ``pip3``):

.. code-block:: console

   pip install rad-tools

.. _installation-source:

Installing from source
======================

* Clone the project to your local computer:

.. code-block:: python

   git clone git@github.com:adrybakov/rad-tools.git

* Change the directory:

.. code-block:: python

   cd rad-tools

* To install RAD-tools, run (you may need to use ``pip3``):

.. code-block:: console

   pip install rad-tools

Additionally you may run the unit tests, which requires |pytest| to be installed:

.. code-block:: console

   make test

.. note::
   pytest requires Python 3.7+

.. hint::
   To install pytest, run (you may need to use ``pip3``):

   .. code-block:: console

      pip install pytest


Update
======

If you want to update the package to the latest available version (|version|)
type the following in your terminal (you may need to use ``pip3``):

.. code-block:: console

   pip install rad-tools --upgrade

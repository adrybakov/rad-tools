*********
RAD-TOOLS
*********
Various scripts for condense matter calculations and post-processing.

.. image:: https://badge.fury.io/py/rad-tools.svg
    :target: https://badge.fury.io/py/rad-tools
    
.. image:: https://readthedocs.org/projects/rad-tools/badge/?version=stable
    :target: https://rad-tools.org/en/stable/?badge=stable
    :alt: Documentation Status
   
.. image:: https://static.pepy.tech/personalized-badge/rad-tools?period=total&units=international_system&left_color=black&right_color=blue&left_text=Downloads
 :target: https://pepy.tech/project/rad-tools

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
   
.. image:: https://img.shields.io/github/license/adrybakov/rad-tools
   :alt: GitHub

The package covers post-processing scenarios for the results of 
`Quantum Espresso <https://www.quantum-espresso.org>`_, 
`TB2J <https://tb2j.readthedocs.io/en/latest/>`_ 
and `Wannier90 <http://www.wannier.org/>`_, as well as provides some custom scripts.

It is expected to be used in two ways:

* As Python library

* Via scripts (i.e. usage from command line)

For the detailed description check
`documentation. <https://rad-tools.org>`_

Installation
============

Requirement for RAD-tools installation are:

* Python itself (>=3.8)
* NumPy
* SciPy
* matplotlib
* tqdm
* termcolor

RAD-tools can be installed with ``pip`` or from source.

pip
---

To install RAD-tools, run (you may need to use ``pip3``):

.. code-block:: console

   pip install rad-tools

Installing from source
----------------------

* Clone the project to your local computer:

.. code-block:: python

   git clone git@github.com:adrybakov/rad-tools.git

* Change the directory:

.. code-block:: python

   cd rad-tools

* To install RAD-tools, run (you may need to use ``pip3``):

.. code-block:: console

   pip install rad-tools

Additionally you may run the unit tests, 
which requires pytest (requires Python 3.7+) to be installed:

.. code-block:: console

   make test

To install pytest, run (you may need to use ``pip3``):

.. code-block:: console

   pip install pytest

.. _announcement:

********************
Change of the domain
********************

**DATE: 20.10.2025**

As of today RAD_tools will be primarily served from the `rad-tools.readthedocs.io <https://rad-tools.readthedocs.io>`_ domain.

`rad-tools.org <https://rad-tools.ord>`_ domain will redirect to the new one and **will be disabled** in the future.

******************************
Release of magnopy and wulfric
******************************

**DATE: 20.05.2025**

Everything related to the Bravais lattices has been separated in a dedicated python package 
called |wulfric|_. 

Everything related to the spin Hamiltonian and magnons is being superseded by the package
called |magnopy|_. 

Starting from RAD-tools v1.0 both those parts are being removed. Across documentation
links to the appropriate documentation pages of either |wulfric|_ or |magnopy|_ has been
left when possible.

The source code of RAD-tools is modified in a way that every function that is removed
still can be called, but will print a message with redirection to the appropriate package.

Those redirect messages in documentation and in the source code will be removed after a
year of the transition period.

.. _guide_io_tb2j:

.. currentmodule:: radtools

*****************
|TB2J|_ interface
*****************

Main output of TB2J is the "exchange.out" file with the exchange parameters.

Function :py:func:`read_tb2j_model` reads this file 
and construct :py:class:`SpinHamiltonian` from it.
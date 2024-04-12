.. _guide_spinham:

****************
Spin Hamiltonian
****************

For the full reference see :ref:`api_spinham`.

.. currentmodule:: radtools

Spin Hamiltonian module describe the logic around the spin Hamiltonian.
:py:class:`.SpinHamiltonian` is the main class of this module. It stores bonds between
:py:class:`.Atom` :math:`i`, which is located in the unit cell :math:`(0, 0, 0)` and
:py:class:`.Atom` :math:`j`, which is located in the unit cell :math:`(i, j, k)` unit cell.

The main building block of the spin Hamiltonian is exchange parameter
(:math:`\boldsymbol{J}`), which is implemented in a separate class
:py:class:`.ExchangeParameter`. For the guide on the classes see:

.. toctree::
    :caption: Individual guides on
    :maxdepth: 1

    hamiltonian
    parameter

The notation of the spin Hamiltonian is an important issue. Without clear notation
exchange parameters does not make much sense. RAD-tools can support any
notation out of the most common ones. :py:class:`.ExchangeParameter` is not aware of the
notation. All logic for the notation is implemented in the :py:class:`.SpinHamiltonian`.

Detailed guide on the notation is available in the :ref:`library_spinham_notation` page

We encourage you to read this guide once with full attention.

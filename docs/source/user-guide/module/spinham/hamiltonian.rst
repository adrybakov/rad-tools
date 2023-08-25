.. _guide_spin-hamiltonian:

***************
SpinHamiltonian
***************

For the full reference see :ref:`api_spin-hamiltonian`.

.. currentmodule:: radtools

Main idea of the spin Hamiltonian structure could be expressed as 

.. code-block:: text

    (atom1, atom2, (i,j,k)) -> exchange_parameter

where ``atom1`` and ``atom2`` are instances of :py:class:`.Atom`, ``(i,j,k)`` is a tuple of
integers, which defines the unit cell (in relative coordinates). ``exchange_parameter`` is an instance of
:py:class:`.ExchangeParameter`. The :py:class:`.SpinHamiltonian` supports 
dictionary-like access to the exchange parameters 
(see :ref:`examples <guide_spin-hamiltonian_dictionary>`).

Import
======

.. doctest::

    >>> # Exact import
    >>> from radtools.spinham.hamiltonian import SpinHamiltonian
    >>> # Explicit import
    >>> from radtools.spinham import SpinHamiltonian
    >>> # Recommended import
    >>> from radtools import SpinHamiltonian

For the examples in this page we need additional import and some predefined variables:

.. doctest::

    >>> from radtools import ExchangeParameter, Atom, Crystal, lattice_example
    >>> Cr = Atom("Cr", spin=1.5)
    >>> Cr1 = Atom("Cr", position=(0.5, 0, 0))
    >>> Cr2 = Atom("Cr", position=(1, 0, 0))
    >>> J1 = ExchangeParameter(iso=1)
    >>> J2 = ExchangeParameter(iso=2)
    >>> tet_lattice = lattice_example("TET")
    >>> crystal = Crystal(a1 = [1, 0, 0], a2 = [0, 1, 0], a3 = [0, 0, 1], atoms=[Cr])

Creation
========

.. hint::
    
    For the creation of the spin Hamiltonian from the |TB2J|_ file see 
    :py:func:`.load_tb2j_model`.

.. doctest::

    >>> hamiltonian = SpinHamiltonian()

:py:class:`.SpinHamiltonian` is a child of :py:class:`.Crystal`. Constructor can take an
instance of the :py:class:`.Crystal` class or keyword arguments, which will be passed to the
:py:class:`.Crystal` constructor:

.. doctest::

    >>> hamiltonian = SpinHamiltonian(crystal)
    >>> # Is equivalent to
    >>> hamiltonian = SpinHamiltonian(a1 = [1, 0, 0], a2 = [0, 1, 0], a3 = [0, 0, 1], atoms=[Cr])

For the full list of constructor parameters see 
:py:class:`.SpinHamiltonian` documentation.

After the creation of the :py:class:`.SpinHamiltonian` you 
need to add atoms and bonds to it.

Adding atoms
============

Atoms can be passed to the constructor of the :py:class:`.SpinHamiltonian` or added later.
Atom addition is inherited from the :py:class:`.Crystal`:

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> len(hamiltonian.atoms)
    0
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(1, 0, 0))
    >>> hamiltonian.add_atom("Br")
    >>> len(hamiltonian.atoms)
    2

Adding bonds
============

The bond is added to the Hamiltonian via :py:meth:`.SpinHamiltonian.add_bond` method:

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.add_bond(Cr, Cr, (1, 0, 0), iso=1)
    >>> hamiltonian["Cr", "Cr", (1, 0, 0)].iso
    1.0

.. note::

    The bond is added to the Hamiltonian, but no atom is explicitly added to the crystal.
    You can skip the atom addition if you are using :py:class:`.Atom` instances, 
    not the names in the :py:meth:`.SpinHamiltonian.add_bond` method.
    The atom will be silently added to the system.
    When atom`s name is used, atom object is extracted from the spin Hamiltonian, thus
    it has to be explicitly added first:

    .. doctest::

        >>> hamiltonian = SpinHamiltonian()
        >>> hamiltonian.add_bond("Cr", "Cr", (1, 0, 0), iso=1)
        Traceback (most recent call last):
        ...
        ValueError: No match found for name = Cr, index = None
        >>> hamiltonian.add_atom(Cr)
        >>> hamiltonian.add_bond("Cr", "Cr", (1, 0, 0), iso=1)
        >>> hamiltonian["Cr", "Cr", (1, 0, 0)].iso
        1.0

Equivalent way to  add bond:

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian[Cr, Cr, (1, 0, 0)] = J1
    >>> hamiltonian[Cr, "Cr", (2, 0, 0)] = J2
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    1.0
    >>> hamiltonian[Cr, Cr, (2, 0, 0)].iso
    2.0

.. note::

    If the atom is already added to the Hamiltonian, then you can use both atom`s name and 
    the instance of the atom class. 
    If you have more then one atom with the same name in the Hamiltonian, 
    then you have to use :py:meth:`.Atom.fullname` instead of the name:

    .. doctest::

        >>> hamiltonian = SpinHamiltonian()
        >>> hamiltonian.add_atom(Cr1)
        >>> hamiltonian.add_atom(Cr2)
        >>> hamiltonian[Cr1, Cr2, (1, 0, 0)] = J2
        >>> hamiltonian[Cr1, Cr2, (1, 0, 0)].iso
        2.0
        >>> hamiltonian["Cr", "Cr", (1, 0, 0)] = J1
        Traceback (most recent call last):
        ...
        ValueError: Multiple matches found for name = Cr, index = None
        >>> hamiltonian["Cr__1", "Cr__2", (1, 0, 0)] = J1
        >>> hamiltonian["Cr__1", "Cr__2", (1, 0, 0)].iso
        1.0

Removing atoms
==============

Removing an atom from the :py:class:`.SpinHamiltonian` is not the same as removing
it from the :py:class:`.Crystal`. When atom is removed from the :py:class:`.SpinHamiltonian`
all bonds, which are connected to it are removed as well. 

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(1, 0, 0))
    >>> hamiltonian.add_bond("Cr", "Cr", (1, 0, 0), iso=1)
    >>> len(hamiltonian)
    1
    >>> len(hamiltonian.atoms)
    1
    >>> hamiltonian.remove_atom("Cr")
    >>> len(hamiltonian)
    0
    >>> len(hamiltonian.atoms)
    0

Removing bonds
==============

To remove the bond use :py:meth:`.SpinHamiltonian.remove_bond` method:

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.add_bond(Cr, Cr, (1, 0, 0), iso=1)
    >>> (Cr, Cr, (1, 0, 0)) in hamiltonian
    True
    >>> hamiltonian.remove_bond("Cr", Cr, (1, 0, 0))
    >>> (Cr, Cr, (1, 0, 0)) in hamiltonian
    False

.. hint::

    Note how we used atom name to remove the bond.

Equivalent way to delete bond mimics the dictionary behaviour:

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.add_bond(Cr, Cr, (1, 0, 0), iso=1)
    >>> (Cr, Cr, (1, 0, 0)) in hamiltonian
    True
    >>> del hamiltonian[Cr, Cr, (1, 0, 0)]
    >>> (Cr, Cr, (1, 0, 0)) in hamiltonian
    False


Magnetic atoms
==============

Atoms, which has at least one bond attached to them are called magnetic atoms.
They could be accessed with :py:attr:`.SpinHamiltonian.magnetic_atoms` property:

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(0, 0, 0))
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(1, 0, 0))
    >>> hamiltonian.magnetic_atoms
    []
    >>> hamiltonian.add_bond("Cr__1", "Cr__1", (0, 1, 0), iso=1)
    >>> for atom in hamiltonian.magnetic_atoms:
    ...     print(atom.fullname)
    ...
    Cr__1
    >>> hamiltonian.add_bond("Cr__1", "Cr__2", (0, 0, 0), iso=1)
    >>> for atom in hamiltonian.magnetic_atoms:
    ...     print(atom.fullname)
    ...
    Cr__1
    Cr__2


Atom name vs instance
=====================

In many situations atom name (:py:attr:`.Atom.name`) or 
fullname (:py:attr:`.Atom.fullname`) can be used with the :py:class:`.SpinHamiltonian`.

You can use name if you have only one atom with this name in the Hamiltonian. Otherwise
you have to use fullname. For example:

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(0, 0, 0))
    >>> # Can use name, since there is only one atom with this name
    >>> hamiltonian.add_bond("Cr", "Cr", (0, 1, 0), iso=1)
    >>> for atom1, atom2, (i, j, k), J in hamiltonian:
    ...     print(atom1.fullname, atom2.fullname, (i, j, k), J.iso)
    ...
    Cr__1 Cr__1 (0, 1, 0) 1.0
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(1, 0, 0))
    >>> # Have to use fullname, since there are two atoms with the same name
    >>> hamiltonian.add_bond("Cr__1", "Cr__1", (0, 1, 0), iso=1)
    >>> hamiltonian.add_bond("Cr__1", "Cr__2", (0, 0, 0), iso=2)
    >>> for atom1, atom2, (i,j,k), J in hamiltonian:
    ...     print(atom1.fullname, atom2.fullname, (i, j, k), J.iso)
    ...
    Cr__1 Cr__1 (0, 1, 0) 1.0
    Cr__1 Cr__2 (0, 0, 0) 2.0

If you want to get an :py:class:`.Atom` instance from the :py:class:`.SpinHamiltonian`
you can use :py:meth:`.SpinHamiltonian.get_atom` method:

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(0, 0, 0))
    >>> # Can use name, since there is only one atom with this name
    >>> hamiltonian.get_atom("Cr").fullname
    'Cr__1'
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(1, 0, 0))
    >>> # Have to use fullname, since there are two atoms with the same name
    >>> hamiltonian.get_atom("Cr__1").fullname
    'Cr__1'
    >>> hamiltonian.get_atom("Cr__2").fullname
    'Cr__2'


Filtering
=========

Spin Hamiltonian could be filtered by distance, template or set of (i, j, k) tuples 
(R vectors).

Use :py:meth:`.SpinHamiltonian.filter` to filter the instance of the
:py:class:`.SpinHamiltonian` and :py:meth:`.SpinHamiltonian.filtered` to get
the filtered copy of the :py:class:`.SpinHamiltonian`:

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.add_atom("Cr", position=(0.25, 0.25, 0))
    >>> hamiltonian.add_atom("Cr", position=(0.75, 0.75, 0))
    >>> bonds = [
    ...     (12, "Cr__1", "Cr__2", (0, 0, 0)),
    ...     (12, "Cr__2", "Cr__1", (0, 0, 0)),
    ...     (12, "Cr__1", "Cr__1", (1, 0, 0)),
    ...     (12, "Cr__1", "Cr__1", (-1, 0, 0)),
    ...     (12, "Cr__2", "Cr__2", (1, 0, 0)),
    ...     (12, "Cr__2", "Cr__2", (-1, 0, 0)),
    ...     (12, "Cr__1", "Cr__1", (0, 2, 0)),
    ...     (12, "Cr__1", "Cr__1", (0, -2, 0)),
    ...     (12, "Cr__2", "Cr__2", (0, 2, 0)),
    ...     (12, "Cr__2", "Cr__2", (0, -2, 0)),
    ...     (12, "Cr__2", "Cr__1", (2, 2, 0)),
    ...     (12, "Cr__1", "Cr__2", (-2, -2, 0)),
    ... ]
    >>> for iso, atom1, atom2, R in bonds:
    ...     hamiltonian.add_bond(atom1, atom2, R, iso=iso)
    >>> for atom1, atom2, R, J in hamiltonian:
    ...     print(atom1.fullname, atom2.fullname, R, J.iso)
    ...
    Cr__1 Cr__2 (0, 0, 0) 12.0
    Cr__2 Cr__1 (0, 0, 0) 12.0
    Cr__1 Cr__1 (1, 0, 0) 12.0
    Cr__1 Cr__1 (-1, 0, 0) 12.0
    Cr__2 Cr__2 (1, 0, 0) 12.0
    Cr__2 Cr__2 (-1, 0, 0) 12.0
    Cr__1 Cr__1 (0, 2, 0) 12.0
    Cr__1 Cr__1 (0, -2, 0) 12.0
    Cr__2 Cr__2 (0, 2, 0) 12.0
    Cr__2 Cr__2 (0, -2, 0) 12.0
    Cr__2 Cr__1 (2, 2, 0) 12.0
    Cr__1 Cr__2 (-2, -2, 0) 12.0
    >>> filtered_hamiltonian = hamiltonian.filtered(max_distance=1)
    >>> for atom1, atom2, R, J in filtered_hamiltonian:
    ...     print(atom1.fullname, atom2.fullname, R, J.iso)
    ...
    Cr__1 Cr__2 (0, 0, 0) 12.0
    Cr__2 Cr__1 (0, 0, 0) 12.0
    Cr__1 Cr__1 (1, 0, 0) 12.0
    Cr__1 Cr__1 (-1, 0, 0) 12.0
    Cr__2 Cr__2 (1, 0, 0) 12.0
    Cr__2 Cr__2 (-1, 0, 0) 12.0
    >>> filtered_hamiltonian = hamiltonian.filtered(min_distance=2.1)
    >>> for atom1, atom2, R, J in filtered_hamiltonian:
    ...     print(atom1.fullname, atom2.fullname, R, J.iso)
    ...
    Cr__2 Cr__1 (2, 2, 0) 12.0
    Cr__1 Cr__2 (-2, -2, 0) 12.0
    >>> filtered_hamiltonian = hamiltonian.filtered(R_vector=[(0, 0, 0), (1, 0, 0)])
    >>> for atom1, atom2, R, J in filtered_hamiltonian:
    ...     print(atom1.fullname, atom2.fullname, R, J.iso)
    ...
    Cr__1 Cr__2 (0, 0, 0) 12.0
    Cr__2 Cr__1 (0, 0, 0) 12.0
    Cr__1 Cr__1 (1, 0, 0) 12.0
    Cr__2 Cr__2 (1, 0, 0) 12.0
    >>> filtered_hamiltonian = hamiltonian.filtered(template=[("Cr__1", "Cr__2", (0, 0, 0))])
    >>> for atom1, atom2, R, J in filtered_hamiltonian:
    ...     print(atom1.fullname, atom2.fullname, R, J.iso)
    ...
    Cr__1 Cr__2 (0, 0, 0) 12.0

.. hint::

    You can combine filtering options together. Only the bonds, which satisfy all the
    conditions are included in the filtered Hamiltonian.


Energy
======

Ferromagnetic energy of the Hamiltonian could be calculated with
:py:meth:`.SpinHamiltonian.ferromagnetic_energy` method:

.. note::

    The notation of the Hamiltonian has to be defined for this method to work.

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(0,0,0))
    >>> hamiltonian.add_atom("Cr", spin=1.5, position=(1,0,0))
    >>> hamiltonian.add_bond("Cr__1", "Cr__1", (0, 1, 0), iso=1)
    >>> hamiltonian.add_bond("Cr__1", "Cr__2", (0, 0, 0), iso=2)
    >>> hamiltonian.ferromagnetic_energy()
    Traceback (most recent call last):
    ...
    radtools.exceptions.NotationError:
    ...
    >>> hamiltonian.notation = "standard"
    >>> hamiltonian.ferromagnetic_energy()
    -13.5

Saving
======

The Hamiltonian could be saved in as a text file with 
:py:func:`.dump_spinham_txt` function:

.. doctest::

    >>> from radtools import dump_spinham_txt 
    >>> hamiltonian = SpinHamiltonian() # doctest: +SKIP
    >>> # Saves the hamiltonian into the file "hamiltonian.txt"
    >>> dump_spinham_txt(hamiltonian, "hamiltonian.txt") # doctest: +SKIP

The format of the file is inspired by the output files of the |TB2J|_ code.
Isotropic exchange is always written. DMI, full matrix and symmetric anisotropic 
exchange can be removed from the output.

It can be serialized with :py:func:`.dump_pickle` function.

.. doctest::

    >>> from radtools import dump_pickle 
    >>> hamiltonian = SpinHamiltonian() # doctest: +SKIP
    >>> # Saves the hamiltonian into the file "hamiltonian.pickle"
    >>> dump_pickle(hamiltonian, "hamiltonian") # doctest: +SKIP

.. hint::

    Note that the file extension is added automatically only for the
    :py:func:`.dump_pickle` function. It is done intentionally to
    emphasize that the Hamiltonian is saved in a python-specific binary file. 

    While the :py:func:`.dump_spinham_txt` offers you complete control over the file name. 

.. _guide_spin-hamiltonian_dictionary:

Dictionary-like behaviour
=========================

Spin Hamiltonian supports the logic of a dictionary, where the keys are the tuples of the 
form:

.. math::

    (\text{atom}_1, \text{atom}_2, (i, j, k))

And the values are the :py:class:`.ExchangeParameter` instances.
Key-value pair is called bond. It describes the exchange interaction between two atoms, 
:math:`\text{atom}_1` is always located in the unit cell :math:`(0, 0, 0)`, 
:math:`\text{atom}_2` is located in the unit cell :math:`(i, j, k)`. 
The coordinate of unit cells are always relative to the lattice vectors.

It supports the following operations:

* ``len(hamiltonian)`` - returns the number of bonds in the Hamiltonian

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> len(hamiltonian)
    0
    >>> hamiltonian.add_bond(Cr, Cr, (1, 0, 0), iso=1)
    >>> len(hamiltonian)
    1
    >>> hamiltonian.add_bond(Cr, Cr, (0, 1, 0), iso=2)
    >>> len(hamiltonian)
    2

* Iteration over the Hamiltonian returns the key and the value

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.add_bond(Cr, Cr, (1, 0, 0), iso=1)
    >>> hamiltonian.add_bond(Cr, Cr, (0, 1, 0), iso=2)
    >>> for atom1, atom2, (i, j, k), J in hamiltonian:
    ...     print(atom1.fullname, atom2.fullname, (i, j, k), J.iso)
    ... 
    Cr__1 Cr__1 (1, 0, 0) 1.0
    Cr__1 Cr__1 (0, 1, 0) 2.0

.. note::

    It does not return only the key as for python dictionaries. It returns the key and the
    value. This is done to simplify the iteration over the Hamiltonian.

* ``in`` syntax returns ``True`` if the bond is present in the Hamiltonian

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.add_bond(Cr, Cr, (1, 0, 0), iso=1)
    >>> (Cr, Cr, (1, 0, 0)) in hamiltonian
    True
    >>> (Cr, Cr, (0, 1, 0)) in hamiltonian
    False

* Getting and setting the value of the bond

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian[Cr, Cr, (1, 0, 0)] = J1
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    1.0

* Deleting the bond

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian[Cr, Cr, (1, 0, 0)] = J1
    >>> (Cr, Cr, (1, 0, 0)) in hamiltonian
    True
    >>> del hamiltonian[Cr, Cr, (1, 0, 0)]
    >>> (Cr, Cr, (1, 0, 0)) in hamiltonian
    False


.. _library_spinham_notation-examples:

Notation examples
=================

For the detailed explanation of the notation see :ref:`library_spinham_notation`.

Setting
-------
The notation could be defined in three ways:

* By setting each individual property. For example:

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.double_counting = True
    >>> hamiltonian.spin_normalized = False
    >>> hamiltonian.factor = -1

* By setting the :py:attr:`.notation` property directly with the tuple of five ``bool``:

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.notation = (True, False, -2)
    >>> hamiltonian.notation
    (True, False, -2.0)
    >>> print(hamiltonian.notation_string)
    H = -2 \sum_{i,j} S_i J_{ij} S_j
    Double counting is present.
    Spin vectors are not normalized.
    >>> hamiltonian.double_counting
    True
    >>> hamiltonian.spin_normalized
    False
    >>> hamiltonian.factor
    -2.0


* By setting the :py:attr:`.notation` property directly with the string:

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.notation = 'standard'
    >>> hamiltonian.notation
    (True, False, -1.0)
    >>> print(hamiltonian.notation_string)
    H = - \sum_{i,j} S_i J_{ij} S_j
    Double counting is present.
    Spin vectors are not normalized.
    >>> hamiltonian.notation = 'tb2j'
    >>> hamiltonian.notation
    (True, True, -1.0)
    >>> print(hamiltonian.notation_string)
    H = - \sum_{i,j} e_i J_{ij} e_j
    Double counting is present.
    Spin vectors are normalized to 1.
    >>> hamiltonian.notation = 'vampire'
    >>> hamiltonian.notation
    (False, True, -1.0)
    >>> print(hamiltonian.notation_string)
    H = - \sum_{i>=j} e_i J_{ij} e_j
    No double counting.
    Spin vectors are normalized to 1.
    >>> hamiltonian.notation = 'spinw'
    >>> hamiltonian.notation
    (True, False, 1.0)
    >>> print(hamiltonian.notation_string)
    H = \sum_{i,j} S_i J_{ij} S_j
    Double counting is present.
    Spin vectors are not normalized.

Changing
--------

Once the full notation or any individual property is set, 
the following redefinition of the notation or corresponding property changes exchange
parameters. For example:

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.add_bond(Cr, Cr, (1, 0, 0), iso=1)
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    1.0
    >>> hamiltonian.notation = (True, False, -1)
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    1.0
    >>> # Normalize spins
    >>> hamiltonian.notation = (True, True, -1)
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    2.25
    >>> # Remove minus sign in the Hamiltonian definition
    >>> hamiltonian.notation = (True, True, 1)
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    -2.25

The rule for the notation change is simple: Conserve the value of energy.

Resetting
---------

If you want to reset the notation once it is set, but keep the parameters intact use 
:py:meth:`.SpinHamiltonian.set_interpretation` method:

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.double_counting = True
    >>> hamiltonian.spin_normalized = False
    >>> hamiltonian.factor = -2
    >>> hamiltonian.notation
    (True, False, -2.0)
    >>> print(hamiltonian.notation_string)
    H = -2 \sum_{i,j} S_i J_{ij} S_j
    Double counting is present.
    Spin vectors are not normalized.
    >>> hamiltonian.set_interpretation(double_counting=False, spin_normalized=True)
    >>> hamiltonian.notation
    (False, True, -2.0)
    >>> print(hamiltonian.notation_string)
    H = -2 \sum_{i>=j} e_i J_{ij} e_j
    No double counting.
    Spin vectors are normalized to 1.

Double counting
---------------

.. doctest::

    >>> hamiltonian = SpinHamiltonian()
    >>> hamiltonian.add_bond(Cr, Cr, (1, 0, 0), iso=1)
    >>> hamiltonian.notation = "standard"
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso, hamiltonian.ferromagnetic_energy()
    (1.0, -4.5)
    >>> hamiltonian.double_counting, hamiltonian.ferromagnetic_energy()
    (True, -4.5)
    >>> hamiltonian.double_counting = False
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    2.0
    >>> hamiltonian.double_counting, hamiltonian.ferromagnetic_energy()
    (False, -4.5)
    >>> hamiltonian.double_counting = True

Spin normalization
------------------

    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    1.0
    >>> hamiltonian.spin_normalized, hamiltonian.ferromagnetic_energy()
    (False, -4.5)
    >>> # Spin of Cr atom is 1.5
    >>> hamiltonian.spin_normalized = True
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    2.25
    >>> hamiltonian.spin_normalized, hamiltonian.ferromagnetic_energy()
    (True, -4.5)
    >>> hamiltonian.spin_normalized = False

Factor
----------

    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    1.0
    >>> hamiltonian.factor, hamiltonian.ferromagnetic_energy()
    (-1.0, -4.5)
    >>> hamiltonian.factor = 1
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    -1.0
    >>> hamiltonian.factor, hamiltonian.ferromagnetic_energy()
    (1.0, -4.5)
    >>> hamiltonian.factor = -2
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    0.5
    >>> hamiltonian.factor, hamiltonian.ferromagnetic_energy()
    (-2.0, -4.5)
    >>> hamiltonian.factor = -0.5
    >>> hamiltonian[Cr, Cr, (1, 0, 0)].iso
    2.0
    >>> hamiltonian.factor, hamiltonian.ferromagnetic_energy()
    (-0.5, -4.5)

Crystal of the spin Hamiltonian
===============================

:py:class:`.SpinHamiltonian` is a child of :py:class:`.Crystal` and inherits all its
properties and methods. Any property, which is related to the crystal is expected to 
be called directly from the :py:class:`.SpinHamiltonian` instance. For example:

.. doctest::

    >>> hamiltonian = SpinHamiltonian(lattice=tet_lattice)
    >>> hamiltonian.variation
    'TET'
    >>> hamiltonian.a1
    array([3.14159265, 0.        , 0.        ])
    >>> hamiltonian.cell
    array([[3.14159265, 0.        , 0.        ],
           [0.        , 3.14159265, 0.        ],
           [0.        , 0.        , 4.71238898]])

The crystal of the Spin Hamiltonian can be access through the 
:py:attr:`.SpinHamiltonian.crystal` attribute. It returns an independent instance of 
the :py:class:`.Crystal` class:

.. doctest::

    >>> hamiltonian = SpinHamiltonian(lattice=tet_lattice)
    >>> crystal = hamiltonian.crystal
    >>> crystal.cell
    array([[3.14159265, 0.        , 0.        ],
           [0.        , 3.14159265, 0.        ],
           [0.        , 0.        , 4.71238898]])
    >>> hamiltonian.cell
    array([[3.14159265, 0.        , 0.        ],
           [0.        , 3.14159265, 0.        ],
           [0.        , 0.        , 4.71238898]])
    >>> crystal.cell = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    >>> crystal.cell
    array([[1, 0, 0],
           [0, 1, 0],
           [0, 0, 1]])
    >>> hamiltonian.cell
    array([[3.14159265, 0.        , 0.        ],
           [0.        , 3.14159265, 0.        ],
           [0.        , 0.        , 4.71238898]])

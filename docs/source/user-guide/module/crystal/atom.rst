.. _guide_crystal_atom:

****
Atom
****

.. currentmodule:: radtools

For the full reference see :ref:`api_atom`.

:py:class:`.Atom` class describe an atom. :py:class:`.Atom` is hashable and can be
used as a dictionary key. The hash is calculated from the atom name and index.

Import
======

.. doctest::

    >>> # Exact import
    >>> from radtools.crystal.atom import Atom
    >>> # Explicit import
    >>> from radtools.crystal import Atom
    >>> # Recommended import
    >>> from radtools import Atom

Creation
========

Creation of an Atom object is straightforward:

.. doctest::

    >>> atom = Atom()
    >>> atom.name
    'X'
    >>> atom.position
    array([0., 0., 0.])
    >>> atom.type
    'X'


For the full list of constructor parameters see
:py:class:`.Atom` documentation.

.. doctest::

    >>> atom1 = Atom(name='X', index=1)
    >>> atom2 = Atom(name='X', index=1)
    >>> atom3 = Atom(name='X', index=2)
    >>> atom1 == atom2
    True
    >>> atom1 == atom3
    False
    >>> atom1 != atom3
    True

Comparison
==========

Two atoms are considered to be equal if they have the same index and name:

.. doctest::

    >>> atom1 = Atom(name='Fe', index=1)
    >>> atom2 = Atom(name='Fe', index=2)
    >>> atom3 = Atom(name='Cr', index=1, position=(1, 1, 0), spin=0.5)
    >>> atom4 = Atom(name='Cr', index=1, position=(1, 0, 0), spin=1.5)
    >>> atom1 == atom2
    False
    >>> atom1 == atom3
    False
    >>> atom1 != atom3
    True
    >>> atom3 == atom4
    True

All other properties of the atom are ignored when comparing atoms.

Usually index is automatically generated when a set of atoms appears in some context.
For example, when atoms are added to the :py:class:`.Crystal` object, the index is
silently assigned to each atom.

Hash of the atom is calculated from the atom name and index.

.. note::

    If the index of the atom is not defined, then you can still compare it to other atoms
    with a different name:

    .. doctest::

        >>> atom1 = Atom(name='Fe')
        >>> atom2 = Atom(name='Cr')
        >>> atom1 == atom2
        False
        >>> atom1 != atom2
        True

    However, if for the pair of atoms the index is not defined in at least one of them,
    then the comparison fails:

    .. doctest::

        >>> atom1 = Atom(name='Fe')
        >>> atom2 = Atom(name='Fe', index=1)
        >>> atom1 == atom2
        Traceback (most recent call last):
        ...
        ValueError: Index is not defined for the atom Fe.

Fullname
========
For the simplicity :py:meth:`.Atom.fullname` property is defined. It returns
a string which is a combination of the atom name and index. It is simply an
atom`s name with two underscores and index appended to it:

.. doctest::

    >>> atom1 = Atom(name='Fe', index=1)
    >>> atom2 = Atom(name='Fe', index=2)
    >>> atom3 = Atom(name='Cr', index=3)
    >>> atom1.fullname
    'Fe__1'
    >>> atom2.fullname
    'Fe__2'
    >>> atom3.fullname
    'Cr__3'

Fullname is defined even if the atom does not have an index:

.. doctest::

    >>> atom = Atom(name='Fe')
    >>> atom.fullname
    'Fe'

Type
====

Atom type is derived from its name:

.. doctest::

    >>> atom.name = 'Cr1'
    >>> atom.name
    'Cr1'
    >>> atom.type
    'Cr'

The atom type is automatically determined from the atom name and can not
be changed directly.

Position
========

The position of the atom can be changed by setting the position attribute:

.. doctest::

    >>> atom = Atom(name="Cr")
    >>> # It has the default value
    >>> atom.position
    array([0., 0., 0.])
    >>> atom.position = [1, 2, 3]
    >>> atom.position
    array([1., 2., 3.])

.. note::

    When :py:class:`.Atom` is used by itself, the position is not considered to be
    in relative or absolute coordinates. The interpretation of atom`s position
    depends on the context. For example, when atom is used in :py:class:`.Crystal`
    object, the position is usually considered to be in relative coordinates.

Spin
====

Spin of the atom is describe by three properties:

* :py:attr:`.spin` - total spin of the atom
* :py:attr:`.spin_vector` - spin vector of the atom
* :py:attr:`.spin_direction` - spin direction of the atom

These three properties are interconnected. Setting each one of them changes the
other two with one exception: :py:attr:`.spin` does not affect :py:attr:`.spin_direction`
and vice versa.

.. doctest::

    >>> # It has the default value
    >>> atom.spin_direction
    array([0., 0., 1.])
    >>> atom.spin = 3
    >>> atom.spin_vector, atom.spin_direction
    (array([0., 0., 3.]), array([0., 0., 1.]))
    >>> atom.spin_vector = [0, 5, 0]
    >>> atom.spin, atom.spin_direction
    (5.0, array([0., 1., 0.]))
    >>> atom.spin_direction = [1, 0, 0]
    >>> atom.spin, atom.spin_vector
    (5.0, array([5., 0., 0.]))

Magnetic moment
===============

Magnetic moment of the atom can be set by assigning a value to the
:py:attr:`.magmom` attribute:

.. doctest::

    >>> atom.magmom = [0, 0, 1]
    >>> atom.magmom
    array([0., 0., 1.])

The units of magnetic moment depend on you interpretation. In RAD-tools usually
Bohr magneton is used.

.. note::

    Magnetic moment is independent from the spin of the atom.
    This behaviour may change in the future.

Charge
======

Electrical charge of the atom can be set by assigning a value to the
:py:attr:`.charge` attribute:

.. doctest::

    >>> atom.charge = 1
    >>> atom.charge
    1.0

The units of magnetic moment depend on you interpretation. In RAD-tools usually
charge of an electron is used.

String representation
=====================
``__str__`` and ``__format__`` methods are defined for the :py:class:`.Atom` class.

.. doctest::

    >>> atom = Atom(name='Fe', index=1, position=(1, 1, 0), spin=0.5)
    >>> print(atom)
    Fe
    >>> print(f"{atom:>10}")
            Fe
    >>> print(f"{atom:#^10}")
    ####Fe####

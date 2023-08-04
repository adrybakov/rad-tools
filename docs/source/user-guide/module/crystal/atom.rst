.. _guide_crystal_atom:

****
Atom
****

For the full reference see :ref:`api_atom`.

Described by the class :py:class:`.Atom`. Atom is hashable and can be used as a
dictionary key. The hash is calculated from the atom name and index.
Two atoms are considered equal if they have the same name and index:

.. doctest::

    >>> import radtools as rad
    >>> atom1 = rad.Atom(name='X', index=1)
    >>> atom2 = rad.Atom(name='X', index=1)
    >>> atom3 = rad.Atom(name='X', index=2)
    >>> atom1 == atom2
    True
    >>> atom1 == atom3
    False
    >>> atom1 != atom3
    True

Creation and basic attributes
=============================

Creation of an Atom object is straightforward:

.. doctest::

    >>> import radtools as rad
    >>> atom = rad.Atom()
    >>> atom.name
    'X'
    >>> atom.position
    array([0, 0, 0])
    >>> atom.type
    'X'

All parameters, that can be passed to the constructor are listed
in the api reference: :py:class:`.Atom`.

The position of the atom can be changed by setting the position attribute:

.. doctest::

    >>> atom.position = [1, 2, 3]
    >>> atom.position
    array([1., 2., 3.])

.. note::

    When Atom class is used by itself, the position is not considered to be
    in relative or absolute coordinates. The interpretation of atom`s position
    depends on the context. For example, when Atom is used in :py:class:`.Crystal`
    object, the position is usually considered to be in relative coordinates.

The atom type is automatically determined from the atom name and can not
be changed directly. The atom name can be changed by setting the name 
attribute:

.. doctest::

    >>> atom.name = 'Cr1'
    >>> atom.name
    'Cr1'
    >>> atom.type
    'Cr'

As you can see the atom type is not necessarily the same as the atom name.

Atom index is an additional identifier that suppose to form a unique
combination with the atom name. In most cases it is automatically generated 
(i.e. when an Atom is added to the :py:class:`.Crystal` object). However,
it can be set manually:

.. doctest::

    >>> atom.index = 1
    >>> atom.index
    1

Properties
==========

Atom object has a few physical properties. Some of them have default values
and some of them has to be set manually. The properties which have to be set:

* :py:attr:`.charge`
* :py:attr:`.magmom`
* :py:attr:`.spin`
* :py:attr:`.spin_vector`

Properties which have default values:

* :py:attr:`.position` : [0, 0, 0]
* :py:attr:`.spin_direction` : [0, 0, 1]

All properties can be set by assigning a value to the corresponding attribute:

.. doctest::

    >>> atom.charge = 1
    >>> atom.magmom = [0, 0, 1]
    >>> atom.spin = 3
    >>> atom.spin_vector = [0, 0, 3]
    >>> atom.position = [1, 2, 3]
    >>> atom.spin_direction = [0, 0, 1]

If the property is set or have a default value, it can be accessed by
calling the corresponding attribute:

.. doctest::

    >>> atom.charge
    1.0
    >>> atom.magmom
    array([0., 0., 1.])
    >>> atom.spin
    3.0
    >>> atom.spin_vector
    array([0., 0., 3.])
    >>> atom.position
    array([1., 2., 3.])
    >>> atom.spin_direction
    array([0., 0., 1.])

:py:attr:`.spin`, :py:attr:`.spin_direction` and :py:attr:`.spin_vector` are 
interconnected. :py:attr:`.spin` and :py:attr:`.spin_direction` can be set
separately, while :py:attr:`.spin_vector` is calculated automatically:

.. doctest::

    >>> atom.spin = 3
    >>> atom.spin_direction = [0, 0, 1]
    >>> atom.spin_vector
    array([0., 0., 3.])

Correspondingly, :py:attr:`.spin` and :py:attr:`.spin_direction` can be calculated
from :py:attr:`.spin_vector`:

.. doctest::

    >>> atom.spin_vector = [0, 5, 0]
    >>> atom.spin
    5.0
    >>> atom.spin_direction
    array([0., 1., 0.])

:py:attr:`.magmom` and :py:attr:`.charge` are independent properties.












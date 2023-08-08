.. _guide_crystal_crystal:

*******
Crystal
*******

.. versionadded:: 0.7

For the full reference see :ref:`api_crystal`

.. currentmodule:: radtools

Crystal is a child of the the :py:class:`.Lattice` and utilize  
:py:class:`.Atom` for the storage of atoms. We recommend you to read 
:ref:`guide_crystal_lattice` and :ref:`guide_crystal_atom` first.

Crystal behaves like a list of atoms. See :ref:`guide_crystal_list-like` for the details.


Import 
======

    >>> # Exact import
    >>> from radtools.crystal.crystal import Crystal
    >>> # Explicit import
    >>> from radtools.crystal import Crystal
    >>> # Recommended import
    >>> from radtools import Crystal

For the examples in this page we need additional import and some predefined variables:

.. doctest::

    >>> from radtools import Lattice, Atom

Creation of the Crystal
=======================

To create a crystal you would typically need to define the lattice and the set of atoms:

* ``lattice`` - instance of the :py:class:`.Lattice` class.
* ``atoms`` - list of the :py:class:`.Atom` instances, which defines the set of atoms in the unit cell.

.. doctest::

    >>> lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0)
    >>> atoms = [Atom('Si', (0.0, 0.0, 0.0))]
    >>> crystal = Crystal(lattice, atoms)

Parameters for lattice creation can be passed to the :py:class:`.Crystal` constructor as 
keyword arguments:

.. doctest::

    >>> crystal = Crystal(a=3.0, b=3.0, c=3.0, alpha=90.0, beta=90.0, gamma=90.0)

By default the :py:class:`.Crystal` is created with the cubic lattice (:math:`a = 1`) and and no atoms.

.. doctest::

    >>> from radtools import Crystal
    >>> crystal = Crystal()
    >>> crystal.cell
    array([[1., 0., 0.],
           [0., 1., 0.],
           [0., 0., 1.]])
    >>> crystal.type()
    'CUB'
    >>> crystal.atoms
    []


Adding atoms
============

Adding atoms to the crystal is straightforward: 

.. doctest::

    >>> crystal = Crystal()
    >>> crystal.add_atom('Cr', position=(0.0, 0.0, 0.0))
    >>> crystal.add_atom('Cr', position=(0.5, 0.5, 0.5))
    >>> for atom in crystal.atoms:
    ...    print(atom.fullname, atom.position)
    ...
    Cr__1 [0. 0. 0.]
    Cr__2 [0.5 0.5 0.5]

By default atom coordinates are interpreted as relative. 
In order to add atom with absolute coordinates set
the ``relative`` keyword argument to ``False``:

.. doctest::

    >>> crystal = Crystal(cell=[[3., 0., 0.], [0., 3., 0.], [0., 0., 3.]])
    >>> crystal.add_atom(Atom('Cr', (0.0, 0.0, 0.0)), relative=False)
    >>> crystal.add_atom(Atom('Cr', (1.5, 1.5, 1.5)), relative=False)
    >>> for atom in crystal.atoms:
    ...    print(atom.fullname, atom.position)
    ...
    Cr__1 [0. 0. 0.]
    Cr__2 [0.5 0.5 0.5]

.. note::   

    :py:meth:`.Crystal.get_distance` and :py:meth:`.Crystal.get_vector` methods 
    are the ones which returns absolute coordinates by default.

It is possible to add atoms by passing the :py:class:`.Atom` instance:

.. doctest::

    >>> crystal = Crystal()
    >>> Cr = Atom('Cr', position=(0.0, 0.0, 0.0))
    >>> crystal.add_atom(Cr)

Removing atoms
==============

Removing atoms is also straightforward:

.. doctest::

    >>> crystal = Crystal()
    >>> crystal.add_atom("Cr", position=(0.0, 0.0, 0.0))
    >>> atom = Atom("Cr", (0.5, 0.5, 0.5))
    >>> crystal.add_atom(atom)
    >>> crystal.remove_atom("Cr__1")
    >>> for atom in crystal.atoms:
    ...    print(atom.fullname, atom.position)
    ...
    Cr__2 [0.5 0.5 0.5]
    >>> crystal.remove_atom(atom)
    >>> crystal.atoms
    []
    >>> crystal.add_atom("Cr", position = (0.0, 0.0, 0.0))
    >>> crystal.add_atom("Cr", position = (1.0, 0.0, 0.0))
    >>> crystal.add_atom("Cr", position = (1.0, 1.0, 0.0))
    >>> crystal.add_atom("Cr", position = (1.0, 1.0, 1.0))
    >>> # Same as len(crystal.atoms)
    >>> len(crystal)
    4
    >>> # You can remove all atoms with the same name at once
    >>> crystal.remove_atom("Cr")
    >>> len(crystal)
    0

.. note::

    The :py:meth:`.Crystal.remove_atom` method can take either the name (or fullname)
    or instance of the atom as an argument.

Atom getters
============

There is a number of properties involving atom of the crystal you have access to:

* :py:attr:`.Crystal.atoms` - list of the :py:class:`.Atom` instances in the crystal.
* :py:attr:`.Crystal.get_atom` - get atom by name or by fullname.

..doctest::


    >>> crystal = Crystal()
    >>> crystal.add_atom('Cr', position=(0.0, 0.0, 0.0))
    >>> atom1 = crystal.get_atom('Cr')
    >>> atom1.name
    'Cr'
    >>> atom1.position
    array([0., 0., 0.])
    >>> atom1.index
    1
    >>> crystal.add_atom('Cr', position=(0.5, 0.5, 0.5))
    >>> atom2 = crystal.get_atom('Cr')
    Traceback (most recent call last):
    ...
    ValueError: Multiple matches found for name = Cr, index = None
    >>> atom2 = crystal.get_atom('Cr__2')
    >>> atom2.index
    2

* :py:attr:`.Crystal.get_atom_coordinates` - get coordinates of the atom by name or by fullname.

.. doctest::

    >>> crystal = Crystal(cell=[[2, 0, 0], [0, 2, 0],[0, 0, 2]])
    >>> crystal.add_atom('Cr', position = (0.0, 0.5, 0.0))
    >>> crystal.get_atom_coordinates('Cr')
    array([0. , 0.5, 0. ])
    >>> crystal.Cr.position
    array([0. , 0.5, 0. ])
    >>> crystal.get_atom_coordinates('Cr', relative=False)
    array([0., 1., 0.])

* :py:attr:`.Crystal.get_vector` - get vector from the atom1 to atom2 by name or by fullname.

.. doctest::


    >>> crystal = Crystal(cell=[[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    >>> crystal.add_atom('Cr', position = (0.0, 0.5, 0.0))
    >>> crystal.add_atom('Cr', position = (0.5, 0.0, 0.0))
    >>> # Gives result in absolute coordinates
    >>> crystal.get_vector('Cr__1', 'Cr__2')
    array([ 1., -1.,  0.])
    >>> crystal.get_vector('Cr__1', 'Cr__2', relative=True)
    array([ 0.5, -0.5,  0. ])

* :py:attr:`.Crystal.get_distance` - get distance between atom1 and atom2 by name or by fullname.

.. doctest::

    >>> crystal = Crystal(cell=[[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    >>> crystal.add_atom('Cr', position=(0.0, 0.5, 0.0))
    >>> crystal.add_atom('Cr', position=(0.5, 0.0, 0.0))
    >>> round(crystal.get_distance('Cr__1', 'Cr__2'), 4)
    1.4142


Magnetic dipole-dipole interaction energy
=========================================

If magnetic moments of crystal`s atoms are defined, the magnetic dipole-dipole interaction energy
can be calculated. 

Two functions are used for this purpose:

* :py:meth:`.Crystal.mag_dipdip_energy` - calculate the energy of the magnetic dipole-dipole interaction.

    It returns energy in meV if atom`s positions are in Angstroms 
    and magnetic moments are in Bohr magnetons.

.. doctest::

    >>> crystal = Crystal(cell=[[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    >>> crystal.add_atom('Cr', position=(0.0, 0.25, 0.0))
    >>> crystal.add_atom('Cr', position=(0.25, 0.0, 0.0))
    >>> energy = crystal.mag_dipdip_energy(1,1,1, progress_bar=False)
    Traceback (most recent call last):
    ...
    ValueError: There are no magnetic atoms in the crystal.
    >>> crystal.Cr__1.magmom = [0, 0, 1]
    >>> crystal.Cr__2.magmom = [0, 0, 1]
    >>> energy = crystal.mag_dipdip_energy(1,1,1, progress_bar=False)
    >>> round(energy, 8)
    0.07591712
    >>> crystal.Cr__1.magmom = [0,1,0]
    >>> energy = crystal.mag_dipdip_energy(1,1,1, progress_bar=False)
    >>> round(energy, 8)
    0.0
    >>> crystal.Cr__2.magmom = [1,0,0]
    >>> energy = crystal.mag_dipdip_energy(1,1,1, progress_bar=False)
    >>> round(energy, 8)
    0.11387568

* :py:meth:`.Crystal.converge_mag_dipdip_energy` - converge the energy of the magnetic dipole-dipole interaction.

.. _guide_crystal_list-like:

List-like behaviour
===================

:py:class:`.Crystal` class supports the logic of a list of atoms. 
The following list-like methods are implemented:

* ``len(crystal)`` - returns the number of atoms in the crystal

.. doctest::

    >>> crystal = Crystal()
    >>> crystal.add_atom('Cr', position=(0.0, 0.0, 0.0))
    >>> crystal.add_atom('Cr', position=(0.5, 0.5, 0.5))
    >>> len(crystal)
    2

* Iteration over the atoms

.. doctest::

    >>> crystal = Crystal()
    >>> crystal.add_atom('Cr', position=(0.0, 0.0, 0.0))
    >>> crystal.add_atom('Cr', position=(0.5, 0.5, 0.5))
    >>> for atom in crystal:
    ...    print(atom.fullname)
    ...
    Cr__1
    Cr__2

* ``in`` syntax returns ``True`` if the atom is present in the Crystal

.. doctest::

    >>> crystal = Crystal()
    >>> atom = Atom('Cr', (0.0, 0.0, 0.0))
    >>> crystal.add_atom(atom)
    >>> atom in crystal
    True
    >>> "Cr" in crystal
    True
    >>> 'Cr__1' in crystal
    True
    >>> 'Cr__3' in crystal
    False
    >>> "Br" in crystal
    False

Attribute and item access
=========================

For the access to the atoms objects two additional interfaces are implemented. 
They both rely on the :py:meth:`.Crystal.get_atom` method.

* Access via attribute

It is identical to passing the name to 
the :py:meth:`.Crystal.get_atom` method with ``return_all=False``.

.. doctest::

    >>> crystal = Crystal()
    >>> crystal.add_atom('Cr', position=(0.0, 0.0, 0.0))
    >>> crystal.add_atom('Cr', position=(0.5, 0.5, 0.5))
    >>> atom = crystal.Cr__1
    >>> atom.fullname
    'Cr__1'
    >>> atom = crystal.Cr__2
    >>> atom.fullname
    'Cr__2'
    >>> crystal.Cr
    Traceback (most recent call last):
    ...
    AttributeError: 'Crystal' object has no attribute 'Cr'

* Access via item

It is identical to passing the name to 
the :py:meth:`.Crystal.get_atom` method with ``return_all=True``.

.. doctest::

    >>> crystal = Crystal()
    >>> crystal.add_atom('Cr', position=(0.0, 0.0, 0.0))
    >>> crystal.add_atom('Cr', position=(0.5, 0.5, 0.5))
    >>> atoms = crystal["Cr__1"]
    >>> atoms[0].fullname
    'Cr__1'
    >>> atoms = crystal["Cr__2"]
    >>> atoms[0].fullname
    'Cr__2'
    >>> atoms = crystal["Cr"]
    >>> for atom in atoms:
    ...    print(atom.fullname)
    ...
    Cr__1
    Cr__2

.. note::
    Second method returns a list of atoms, while the first one returns a single atom.






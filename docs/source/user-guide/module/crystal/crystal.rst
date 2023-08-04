.. _guide_crystal_crystal:

*******
Crystal
*******

.. versionadded:: 0.7

For the full reference see :ref:`api_crystal`

.. currentmodule:: radtools

Crystal is a child of the the :py:class:`.Lattice` and utilize  
:py:class:`.Atom` for the storage of atoms.


Creation of the Crystal
=======================

To create a crystal you would typically need to define the lattice and the set of atoms:

* ``lattice`` - instance of the :py:class:`.Lattice` class.
* ``atoms`` - list of the :py:class:`.Atom` instances, which defines the set of atoms in the unit cell.

.. doctest::

    >>> from radtools import Crystal, Lattice, Atom
    >>> lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0)
    >>> atoms = [Atom('Si', (0.0, 0.0, 0.0))]
    >>> crystal = Crystal(lattice, atoms)

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

Identification of the Bravais lattice type
==========================================

Identification of the Bravais lattice type is in the same way as 
for the :ref:`guide_crystal_lattice`. The only difference is that
you may want to identify the primitive unit cell of the crystal first. By default
the current unit cell is considered as primitive.

.. warning::

    Primitive cell identification is not implemented yet. 
    Any cell is interpreted as primitive.

Adding and removing atoms
=========================

Adding atoms to the crystal is straightforward: 

.. doctest::

    >>> from radtools import Crystal
    >>> crystal = Crystal()
    >>> crystal.add_atom('Cr', position=(0.0, 0.0, 0.0))
    >>> crystal.add_atom('Cr', position=(0.5, 0.5, 0.5))
    >>> for atom in crystal.atoms:
    ...    print(atom.fullname, atom.position)
    ...
    Cr__1 [0. 0. 0.]
    Cr__2 [0.5 0.5 0.5]

By default atom coordinates are interpreted as relative. 
In order to add atom in the fractional coordinates
the ``relative`` keyword argument should be set to ``False``:

.. doctest::

    >>> from radtools import Crystal, Lattice, Atom
    >>> lattice = Lattice(3.0, 3.0, 3.0, 90.0, 90.0, 90.0)
    >>> crystal = Crystal(lattice)
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

Removing atoms is also straightforward:

.. doctest::

    >>> from radtools import Crystal, Lattice, Atom
    >>> crystal = Crystal()
    >>> crystal.add_atom(Atom("Cr", (0.0, 0.0, 0.0)))
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

    The :py:meth:`.Crystal.remove_atom` method can take either the string literal
    or the atom instance as an argument. String literal is interpreted as the
    :py:attr:`.Atom.name` if only one atom with that name is present in the crystal.
    Otherwise it is interpreted as the :py:attr:`.Atom.fullname`.

Atom getters
============

There is a number of properties involving atom of the crystal you have access to:

* :py:attr:`.Crystal.atoms` - list of the :py:class:`.Atom` instances in the crystal.
* :py:attr:`.Crystal.get_atom` - get atom by name or by fullname.

..doctest::

    >>> from radtools import Crystal, Atom
    >>> crystal = Crystal()
    >>> crystal.add_atom(Atom('Cr', (0.0, 0.0, 0.0)))
    >>> atom1 = crystal.get_atom('Cr')
    >>> atom1.name
    'Cr'
    >>> atom1.position
    array([0., 0., 0.])
    >>> atom1.index
    1
    >>> crystal.add_atom(Atom('Cr', (0.5, 0.5, 0.5)))
    >>> atom2 = crystal.get_atom('Cr')
    Traceback (most recent call last):
    ...
    ValueError: Multiple matches found for name = Cr, index = None
    >>> atom2 = crystal.get_atom('Cr__2')
    >>> atom2.index
    2

* :py:attr:`.Crystal.get_atom_coordinates` - get coordinates of the atom by name or by fullname.

.. doctest::

    >>> from radtools import Crystal, Lattice
    >>> lattice = Lattice([[2,0,0],[0,2,0],[0,0,2]]) 
    >>> crystal = Crystal(lattice)
    >>> crystal.add_atom('Cr', position = (0.0, 0.5, 0.0))
    >>> crystal.get_atom_coordinates('Cr')
    array([0. , 0.5, 0. ])
    >>> crystal.Cr.position
    array([0. , 0.5, 0. ])
    >>> crystal.get_atom_coordinates('Cr', relative=False)
    array([0., 1., 0.])

* :py:attr:`.Crystal.get_vector` - get vector from the atom1 to atom2 by name or by fullname.

.. doctest::

    >>> from radtools import Crystal
    >>> crystal = Crystal(cell=[[2,0,0],[0,2,0],[0,0,2]])
    >>> crystal.add_atom('Cr', position = (0.0, 0.5, 0.0))
    >>> crystal.add_atom('Cr', position = (0.5, 0.0, 0.0))
    >>> # Gives result in absolute coordinates
    >>> crystal.get_vector('Cr__1', 'Cr__2')
    array([ 1., -1.,  0.])
    >>> crystal.get_vector('Cr__1', 'Cr__2', relative=True)
    array([ 0.5, -0.5,  0. ])

* :py:attr:`.Crystal.get_distance` - get distance between atom1 and atom2 by name or by fullname.

.. doctest::

    >>> from radtools import Crystal
    >>> crystal = Crystal(cell=[[2,0,0],[0,2,0],[0,0,2]])
    >>> crystal.add_atom('Cr', position=(0.0, 0.5, 0.0))
    >>> crystal.add_atom('Cr', position=(0.5, 0.0, 0.0))
    >>> print(f"{crystal.get_distance('Cr__1', 'Cr__2'):.4f}")
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

    >>> from radtools import Crystal
    >>> crystal = Crystal(cell=[[2,0,0],[0,2,0],[0,0,2]])
    >>> crystal.add_atom('Cr', position=(0.0, 0.25, 0.0))
    >>> crystal.add_atom('Cr', position=(0.25, 0.0, 0.0))
    >>> energy = crystal.mag_dipdip_energy(1,1,1, progress_bar=False)
    Traceback (most recent call last):
    ...
    ValueError: There are no magnetic atoms in the crystal.
    >>> crystal.Cr__1.magmom = [0,0,1]
    >>> crystal.Cr__2.magmom = [0,0,1]
    >>> energy = crystal.mag_dipdip_energy(1,1,1, progress_bar=False)
    >>> print(f"{energy:.8f}")
    0.07591712
    >>> crystal.Cr__1.magmom = [0,1,0]
    >>> energy = crystal.mag_dipdip_energy(1,1,1, progress_bar=False)
    >>> print(f"{energy:.8f}")
    0.00000000
    >>> crystal.Cr__2.magmom = [1,0,0]
    >>> energy = crystal.mag_dipdip_energy(1,1,1, progress_bar=False)
    >>> print(f"{energy:.8f}")
    0.11387568

* :py:meth:`.Crystal.converge_mag_dipdip_energy` - converge the energy of the magnetic dipole-dipole interaction.

Iterations and stuff
====================

:py:class:`.Crystal` class is iterable over the atoms of the crystal.

.. doctest::

    >>> from radtools import Crystal
    >>> crystal = Crystal()
    >>> crystal.add_atom('Cr', position=(0.0, 0.0, 0.0))
    >>> crystal.add_atom('Cr', position=(0.5, 0.5, 0.5))
    >>> for atom in crystal:
    ...    print(atom.fullname)
    ...
    Cr__1
    Cr__2

You can check if the atom is in the crystal:

.. doctest::

    >>> from radtools import Crystal, Atom
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

You can get an atom object by its name or fullname 
(it is identical to passing the name to 
the :py:meth:`.Crystal.get_atom` method with ``return_all=False``):

.. doctest::

    >>> from radtools import Crystal
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

Similar method
(it is identical to passing the name to 
the :py:meth:`.Crystal.get_atom` method with ``return_all=True``):

.. doctest::

    >>> from radtools import Crystal
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






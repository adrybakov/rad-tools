from math import sqrt

import numpy as np

from radtools.crystal.atom import Atom
from radtools.crystal.lattice import Lattice
from radtools.crystal.bravais_lattice import bravais_lattice_from_cell
from radtools.routines import absolute_to_relative


class Crystal:
    r"""
    Crystal class.

    Iterable over atoms. All attributes of the :py:class:`.Lattice`
    are accessible directly from the crystal or from the lattice attribute,
    i.e. the following lines are equivalent:

    .. doctest::

        >>> import radtools as rad
        >>> cub = rad.lattice_example("CUB")
        >>> crystal = rad.Crystal(cub)
        >>> crystal.lattice.pearson_symbol
        'cP'
        >>> crystal.pearson_symbol
        'cP'

    For the full description of the lattice attributes and methods
    see :ref:`rad-tools_lattice`.

    Parameters
    ----------
    lattice : :py:class:`.Lattice`, optional
        Lattice of the crystal. If not provided,
        then orthonormal lattice is used ("CUB with :math:`a = 1`).
    atoms : list, optional
        List of :py:class:`Atom` objects.

    Attributes
    ----------
    lattice : :py:class:`.Lattice`
    atoms : list
        List of atoms of the crystal.
    """

    def __init__(self, lattice: Lattice = None, atoms=None) -> None:
        self._lattice = None
        self.atoms = []

        if lattice is None:
            lattice = Lattice([1, 0, 0], [0, 1, 0], [0, 0, 1])
        self.lattice = lattice
        if atoms is not None:
            for a in atoms:
                self.add_atom(a)

    @property
    def cell(self):
        r"""
        Cell of the lattice.

        Main difference of the cell attribute of the crystal from the cell
        attribute of the lattice is that change of the crystal`s cell
        leads to the computation of the atom coordinate.

        Examples
        --------

        .. doctest::

            >>> import radtools as rad
            >>> c = rad.Crystal()
            >>> c.add_atom(rad.Atom("Cr", (0.5, 0.5, 0.5)))
            >>> c.cell
            array([[1, 0, 0],
                   [0, 1, 0],
                   [0, 0, 1]])
            >>> c.atoms[0].position
            array([0.5, 0.5, 0.5])
            >>> c.lattice.cell = [[2,0,0],[0,2,0],[0,0,2]]
            >>> c.cell
            array([[2, 0, 0],
                   [0, 2, 0],
                   [0, 0, 2]])
            >>> c.atoms[0].position
            array([0.5, 0.5, 0.5])
            >>> c.lattice.cell = [[1,0,0],[0,1,0],[0,0,1]]
            >>> c.cell = [[2,0,0],[0,2,0],[0,0,2]]
            >>> c.cell
            array([[2, 0, 0],
                   [0, 2, 0],
                   [0, 0, 2]])
            >>> c.atoms[0].position
            array([1., 1., 1.])
        """
        return self.lattice.cell

    @cell.setter
    def cell(self, new_cell):
        old_cell = self.cell
        self.lattice.cell = new_cell
        for atom in self:
            relative = absolute_to_relative(old_cell, atom.position)
            atom.position = relative @ self.cell

    def __iter__(self):
        return CrystalIterator(self)

    def __contains__(self, atom):
        return atom in self.atoms

    def __getitem__(self, index) -> Atom:
        return self.atoms[index]

    def __getattr__(self, name):
        # Fix copy/deepcopy RecursionError
        if name in ["__setstate__"]:
            raise AttributeError(name)
        return getattr(self.lattice, name)

    @property
    def lattice(self):
        r"""
        Lattice of the crystal.

        See :ref:`rad-tools_lattice` for details.
        """
        return self._lattice

    @lattice.setter
    def lattice(self, new_lattice: Lattice):
        if not isinstance(new_lattice, Lattice):
            raise TypeError(
                "New lattice is not a lattice. " + f"Received {type(new_lattice)}."
            )
        self._lattice = new_lattice

    def add_atom(self, new_atom: Atom):
        if not isinstance(new_atom, Atom):
            raise TypeError("New atom is not an Atom. " + f"Received {type(new_atom)}.")
        if new_atom not in self.atoms:
            try:
                i = new_atom.index
            except ValueError:
                new_atom.index = len(self.atoms) + 1
            self.atoms.append(new_atom)

    def remove_atom(self, atom: Atom):
        r"""
        Remove atom from the crystal.

        Parameters
        ----------
        atom : :py:class:`.Atom`
            Atom object.
        """

        self.atoms.remove(atom)

    def get_atom(self, literal, index=None):
        r"""
        Return atom of the crystal.

        Return all atoms with the same literal. If index is provided,
        then return all atoms with literal and index.

        Notes
        -----
        ``index`` is supposed to be a unique value,
        however it uniqueness is not strictly checked,
        pay attention in custom cases.

        Parameters
        ----------
        literal : str
            Name of the atom. In general not unique.
        index : any, optional
            Index of the atom.

        Returns
        -------
        atom : :py:class:`.Atom` or list or None
            If there is no atom in the crystal then return ``None``.
            If only one atom is found, then :py:class:`.Atom` object is returned.
            If several atoms are found then list of :py:class:`.Atom` objects is returned.
        """

        atoms = []

        for atom in self:
            if atom.literal == literal:
                if index is None:
                    atoms.append(atom)
                elif atom.index == index:
                    atoms.append(atom)
        if len(atoms) == 0:
            return None
        elif len(atoms) == 1:
            return atoms[0]
        return atoms

    def get_atom_coordinates(self, atom: Atom, R=(0, 0, 0)):
        r"""
        Getter for the atom coordinates.

        Parameters
        ----------
        atom : :py:class:`.Atom`
            Atom object.
        R : (3,) |array_like|_, default (0, 0, 0)
            Radius vector of the unit cell for atom2 (i,j,k).

        Returns
        -------
        coordinates : 1 x 3 array
            Coordinates of atom in the cell R in absolute coordinates.
        """

        if atom not in self.atoms:
            raise ValueError(f"There is no {atom} in the crystal.")

        return np.array(R) @ self.lattice.cell + atom.position

    def get_vector(self, atom1, atom2, R=(0, 0, 0)):
        r"""
        Getter for vector between the atom1 and atom2.

        Parameters
        ----------
        atom1 : :py:class:`.Atom`
            Atom object in (0, 0, 0) unit cell.
        atom2 : :py:class:`.Atom`
            Atom object in R unit cell.
        R : (3,) |array_like|_, default (0, 0, 0)
            Radius vector of the unit cell for atom2 (i,j,k).

        Returns
        -------
        v : (3,) :numpy:`ndarray`
            Vector from atom1 in (0,0,0) cell to atom2 in R cell.
        """

        atom1 = self.get_atom_coordinates(atom1)
        atom2 = self.get_atom_coordinates(atom2, R)

        return atom2 - atom1

    def get_distance(self, atom1, atom2, R=(0, 0, 0)):
        r"""
        Getter for distance between the atom1 and atom2.

        Parameters
        ----------
        atom1 : :py:class:`.Atom`
            Atom object in (0, 0, 0) unit cell.
        atom2 : :py:class:`.Atom`
            Atom object in R unit cell.
        R : (3,) |array_like|_, default (0, 0, 0)
            Radius vector of the unit cell for atom2 (i,j,k).

        Returns
        -------
        distance : floats
            Distance between atom1 in (0,0,0) cell and atom2 in R cell.
        """

        return sqrt(np.sum(self.get_vector(atom1, atom2, R) ** 2))

    def find_primitive_cell(self):
        r"""
        Detect primitive cell.
        """
        pass

    def identify(self, find_primitive=True):
        r"""
        Identify Bravais lattice type.

        Parameters
        ----------
        find_primitive : bool, default True
            Whether to find primitive cell before identification.
        """

        # Define primitive cell
        if find_primitive:
            self.find_primitive_cell()

        self.lattice = bravais_lattice_from_cell(self.lattice.cell)


class CrystalIterator:
    def __init__(self, crystal: Crystal) -> None:
        self._list = crystal.atoms
        self._index = 0

    def __next__(self) -> Atom:
        if self._index < len(self._list):
            result = self._list[self._index]
            self._index += 1
            return result
        raise StopIteration

    def __iter__(self):
        return self


if __name__ == "__main__":
    l = Lattice([[2, 0, 0], [0, 3, 0], [0, 0, 1]])
    c = Crystal(l)
    c.identify()
    c.plot("primitive")
    c.plot("brillouin_kpath")
    c.show()

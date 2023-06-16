from math import floor, log10, sqrt
from typing import Union

import numpy as np

from radtools.crystal.atom import Atom
from radtools.crystal.bravais_lattice import bravais_lattice_from_cell
from radtools.crystal.lattice import Lattice
from radtools.crystal.properties import dipole_dipole_energy, dipole_dipole_interaction
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

        Returns
        -------
        cell : (3, 3) :numpy:`ndarray`
            Cell of the crystal`s lattice.


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

    def __contains__(self, atom: Union[Atom, str]):
        if isinstance(atom, str):
            atom = self.get_atom(atom, return_all=True)
            if isinstance(atom, list):
                atom = atom[0]
        return atom in self.atoms

    def __len__(self):
        return self.atoms.__len__()

    def __getitem__(self, index) -> Atom:
        return self.atoms[index]

    def __getattr__(self, name):
        # Fix copy/deepcopy RecursionError
        if name in ["__setstate__"]:
            raise AttributeError(name)
        try:
            index = None
            if "_" in name:
                index = int(name.split("_")[1])
                name = name.split("_")[0]
            atom = self.get_atom(name=name, index=index)
            return atom
        except ValueError:
            return getattr(self.lattice, name)

    @property
    def lattice(self):
        r"""
        Lattice of the crystal.

        See :ref:`rad-tools_lattice` for details.

        Returns
        -------
        lattice : :py:class:`.Lattice`
            Lattice of the crystal.
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
        r"""
        Add atom to the crystall.

        If ``new_atom```s literal and index are the same as of some atom of the crystal,
        then ``new_atom`` is not added.

        If index of ``new_atom`` is not defined, it is set.

        Parameters
        ----------
        new_atoms : :py:class:`.Atom`
            New atom.

        Raises
        ------
        TypeError
            If ``new_atom`` is not an :py:class:`.Atom`.
        ValueError
            If the atom is already present in the crystal.
        """
        if not isinstance(new_atom, Atom):
            raise TypeError("New atom is not an Atom. " + f"Received {type(new_atom)}.")
        try:
            i = new_atom.index
        except ValueError:
            new_atom.index = len(self.atoms) + 1
        if new_atom not in self.atoms:
            self.atoms.append(new_atom)
        else:
            raise ValueError("Atom is already in the crystal.")

    def remove_atom(self, atom: Union[Atom, str]):
        r"""
        Remove atom from the crystal.

        If type(``atom``) == ``str``, then all atoms with the name ``atom`` are removed.

        Parameters
        ----------
        atom : :py:class:`.Atom` or str
            :py:class`.Atom` object or atom`s name.
            If name, then it has to be unique among atoms of the crystal.
        """

        if isinstance(atom, str):
            atoms = self.get_atom(atom, return_all=True)
        else:
            atoms = [atom]
        for atom in atoms:
            self.atoms.remove(atom)

    def get_atom(self, name, index=None, return_all=False):
        r"""
        Return atom object of the crystal.

        Notes
        -----
        :py:attr:`.index` in combination with :py:attr:`.name` is supposed to be a unique value,
        however it uniqueness is not strictly checked,
        pay attention in custom cases.

        Parameters
        ----------
        name : str
            Name of the atom. In general not unique.
        index : any, optional
            Index of the atom.
        return_all : bool, default False
            Whether to return the list of non-unique matches or raise an ``ValueError``.

        Returns
        -------
        atom : :py:class:`.Atom` or list
            If only one atom is found, then :py:class:`.Atom` object is returned.
            If several atoms are found and ``return_all`` is ``True``,
            then list of :py:class:`.Atom` objects is returned.

        Raises
        ------
        ValueError
            If no match is found or the match is not unique and ``return_all`` is False.
        """

        atoms = []

        for atom in self:
            if atom.name == name:
                if index is None or atom.index == index:
                    atoms.append(atom)
        if len(atoms) == 0:
            raise ValueError(f"No match found for name = {name}, index = {index}")
        elif len(atoms) == 1:
            return atoms[0]
        elif not return_all:
            raise ValueError(
                f"Multiple matches found for name = {name}, index = {index}"
            )
        return atoms

    def get_atom_coordinates(self, atom: Union[Atom, str], R=(0, 0, 0)):
        r"""
        Getter for the atom coordinates.

        Parameters
        ----------
        atom : :py:class:`.Atom` or str
            :py:class`.Atom` object or atom`s name.
            If name, then it has to be unique among atoms of the crystal.
        R : (3,) |array_like|_, default (0, 0, 0)
            Radius vector of the unit cell for atom2 (i,j,k).

        Returns
        -------
        coordinates : 1 x 3 array
            Coordinates of atom in the cell R in absolute coordinates.
        """

        if isinstance(atom, str):
            atom = self.get_atom(atom)

        if atom not in self.atoms:
            raise ValueError(f"There is no {atom} in the crystal.")

        return np.array(R) @ self.lattice.cell + atom.position

    def get_vector(self, atom1: Union[Atom, str], atom2: Union[Atom, str], R=(0, 0, 0)):
        r"""
        Getter for vector between the atom1 and atom2.

        Parameters
        ----------
        atom1 : :py:class:`.Atom` or str
            :py:class`.Atom` object or atom`s name in (0, 0, 0) unit cell.
            If name, then it has to be unique among atoms of the crystal.
        atom2 : :py:class:`.Atom` or str
            :py:class`.Atom` object or atom`s name in ``R`` unit cell.
            If name, then it has to be unique among atoms of the crystal.
        R : (3,) |array_like|_, default (0, 0, 0)
            Radius vector of the unit cell for atom2 (i,j,k).

        Returns
        -------
        v : (3,) :numpy:`ndarray`
            Vector from atom1 in (0,0,0) cell to atom2 in R cell.
        """

        coord1 = self.get_atom_coordinates(atom1)
        coord2 = self.get_atom_coordinates(atom2, R)

        return coord2 - coord1

    def get_distance(
        self, atom1: Union[Atom, str], atom2: Union[Atom, str], R=(0, 0, 0)
    ):
        r"""
        Getter for distance between the atom1 and atom2.

        Parameters
        ----------
        atom1 : :py:class:`.Atom` or str
            :py:class`.Atom` object or atom`s name in (0, 0, 0) unit cell.
            If name, then it has to be unique among atoms of the crystal.
        atom2 : :py:class:`.Atom` or str
            :py:class`.Atom` object or atom`s name in ``R`` unit cell.
            If name, then it has to be unique among atoms of the crystal.
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

    def mag_dipdip_energy(self, na, nb, nc, progress_bar=True):
        r"""
        Computes magnetic dipole-dipole energy.

        .. math::

            E = -\frac{\mu_0}{4\pi}\sum_{i > j}\left(3(\vec{m_i} \cdot \vec{r_{ij}})(\vec{m_j} \cdot \vec{r_{ij}}) - (\vec{m_i}\cdot\vec{m_j})\right)

        Parameters
        ----------
        na : int
            Translations along :py:attr:`.Crystal.a1`.
        nb : int
            Translations along :py:attr:`.Crystal.a2`.
        nc : int
            Translations along :py:attr:`.Crystal.a3`.
        progress_bar : bool, default True
            Whether to show progressbar.

        Returns
        -------
        energy : float
            Dipole-dipole energy in meV.

        See Also
        --------
        dipole_dipole_energy
        """

        translations = (
            np.transpose(np.indices((na, nb, nc)), (1, 2, 3, 0)).reshape(
                (na * nb * nc, 3)
            )
            @ self.cell
        )
        n_t = len(translations)

        magnetic_atoms = []
        for atom in self:
            try:
                tmp = atom.magmom
                magnetic_atoms.append(atom)
            except ValueError:
                pass

        magnetic_centres = np.zeros((n_t * len(magnetic_atoms), 2, 3), dtype=float)
        for a_i, atom in enumerate(magnetic_atoms):
            magnetic_centres[a_i * n_t : (a_i + 1) * n_t, 0] = np.tile(
                atom.magmom, (n_t, 1)
            )
            magnetic_centres[a_i * n_t : (a_i + 1) * n_t, 1] = (
                translations + atom.position
            )

        return dipole_dipole_energy(magnetic_centres, progress_bar=progress_bar)

    def converge_mag_dipdip_energy(
        self,
        start=(10, 10, 10),
        step=(10, 10, 10),
        eps=10e-3,
        progress_bar=True,
        verbose=False,
    ):
        r"""
        Converge magnetic dipole-dipole energy.

        Parameters
        ----------
        start : tuple of 3 int, default (10, 10, 10)
            (na, nb, nc) starting values.
        step : tuple of 3 int, default (10, 10, 10)
            (na, nb, nc) step values. If 0, then no step is applied.
        eps : float, default 1e-3
            Convergence parameter.
        progress_bar : bool, default False
            Whether to show progressbar.
        verbose : bool, default False
            Whether to print information about each step.

        Returns
        -------
        energies : list
            List of (na,nb,nc) and energies for each step.

            .. code-block:: python

                energies = [((na, nb, nc), energy), ...]

        See Also
        --------
        dipole_dipole_energy
        dipole_dipole_interaction
        """

        na, nb, nc = start
        da, db, dc = step

        translations = (
            np.transpose(np.indices((na, nb, nc)), (1, 2, 3, 0)).reshape(
                (na * nb * nc, 3)
            )
            @ self.cell
        )

        magnetic_atoms = []
        for atom in self:
            try:
                tmp = atom.magmom
                magnetic_atoms.append(atom)
            except ValueError:
                pass
        n_atoms = len(magnetic_atoms)

        n_t = len(translations)
        mc1 = np.zeros((n_t * n_atoms, 2, 3), dtype=float)
        for a_i, atom in enumerate(magnetic_atoms):
            mc1[a_i * n_t : (a_i + 1) * n_t, 0] = np.tile(atom.magmom, (n_t, 1))
            mc1[a_i * n_t : (a_i + 1) * n_t, 1] = translations + atom.position

        energy1 = dipole_dipole_energy(mc1, progress_bar=progress_bar, normalize=False)
        energies = [(start, energy1 / len(mc1))]
        if verbose:
            n_digits = abs(floor(log10(abs(eps)))) + 2
            print(
                f"{'Size':^{6+len(str(na))+len(str(nb))+len(str(nc))}} {'Energy':^{n_digits+4}} {'Difference':^{n_digits+4}}"
            )
            print(
                f"({na}, {nb}, {nc}) "
                + f"{energy1 /len(mc1):<{n_digits+4}.{n_digits}f}"
            )

        if da == 0 and db == 0 and dc == 0:
            return energies

        difference = 10 * eps
        while difference > eps:
            translations = []
            if dc != 0:
                for k in range(nc, nc + dc):
                    for j in range(0, nb + db):
                        for i in range(0, na + da):
                            translations.append((i, j, k))
            for k in range(0, nc):
                if db != 0:
                    for j in range(nb, nb + db):
                        for i in range(0, na + da):
                            translations.append((i, j, k))
                for j in range(0, nb):
                    for i in range(na, na + da):
                        translations.append((i, j, k))
            translations = np.array(translations) @ self.cell

            n_t = len(translations)
            mc2 = np.zeros((n_t * n_atoms, 2, 3), dtype=float)
            for a_i, atom in enumerate(magnetic_atoms):
                mc2[a_i * n_t : (a_i + 1) * n_t, 0] = np.tile(atom.magmom, (n_t, 1))
                mc2[a_i * n_t : (a_i + 1) * n_t, 1] = translations + atom.position
            energy2 = dipole_dipole_energy(
                mc2, progress_bar=progress_bar, normalize=False
            )
            energy_inter = dipole_dipole_interaction(
                mc1, mc2, progress_bar=progress_bar, normalize=False
            )
            new_energy = energy1 + energy_inter + energy2

            energy1 = new_energy
            nc += dc
            nb += db
            na += da
            energies.append(((na, nb, nc), energy1 / (len(mc2) + len(mc1))))
            difference = abs(energies[-2][1] - energies[-1][1])
            if verbose:
                print(
                    f"({na}, {nb}, {nc}) "
                    + f"{energy1 / (len(mc2) +len(mc1)):<{n_digits+4}.{n_digits}f} "
                    + f"{difference:<{n_digits+4}.{n_digits}f}"
                )

            mc1 = np.concatenate((mc1, mc2))

        return energies


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
    crystal = Crystal()
    crystal.cell = [
        [3.513012, 0.000000000, 0.000000000],
        [0.000000000, 4.752699, 0.000000000],
        [0.000000000, 0.000000000, 23.5428674],
    ]

    Cr1 = Atom(
        "Cr1",
        np.array((0.2500000000, 0.7500000000, 0.1616620796)) @ crystal.cell,
    )
    Cr2 = Atom(
        "Cr2",
        np.array((0.7500000000, 0.2500000000, 0.0737774367)) @ crystal.cell,
    )

    crystal.add_atom(Cr1)
    crystal.add_atom(Cr2)

    crystal.Cr1.magmom = [0, 0, 3]
    crystal.Cr2.magmom = [0, 0, 3]

    z10 = crystal.mag_dipdip_energy(10, 10, 1)
    z15 = crystal.mag_dipdip_energy(15, 15, 1)
    z20 = crystal.mag_dipdip_energy(20, 20, 1)
    print(z10, z15, z20, sep="\n")
    energies = crystal.converge_mag_dipdip_energy(
        (10, 10, 1), (5, 5, 0), eps=0.001, verbose=1, progress_bar=False
    )
    for e in energies:
        print(f"{e[0]} ({e[1]})<------")

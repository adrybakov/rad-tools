from math import floor, log10
from typing import Union

import numpy as np

import radtools.crystal.cell as Cell
from radtools.crystal.atom import Atom
from radtools.crystal.lattice import Lattice
from radtools.crystal.properties import dipole_dipole_energy, dipole_dipole_interaction
from radtools.geometry import absolute_to_relative


class Crystal(Lattice):
    r"""
    Crystal class.

    It is a child of Lattice class.

    Iterable over atoms. All attributes of the :py:class:`.Lattice`
    are accessible directly from the crystal:

    .. doctest::

        >>> import radtools as rad
        >>> cub = rad.lattice_example("CUB")
        >>> crystal = rad.Crystal(cub)
        >>> crystal.pearson_symbol
        'cP'

    For the full description of the lattice attributes and methods
    see :ref:`guide_crystal_lattice`.

    Parameters
    ----------
    lattice : :py:class:`.Lattice`, optional
        Lattice of the crystal. If not provided,
        then orthonormal lattice is used ("CUB with :math:`a = 1`).
    atoms : list, optional
        List of :py:class:`Atom` objects.
    relative : bool, default True
        Whether ``atoms`` positions are in relative coordinates.
    standardize : bool, default True
        Whether to standardize the lattice.
    **kwargs
        Keyword arguments for :py:class:`.Lattice` initialization.

    Attributes
    ----------
    atoms : list
        List of atoms of the crystal.
    """

    def __init__(
        self,
        lattice: Lattice = None,
        atoms=None,
        relative=True,
        standardize=True,
        **kwargs,
    ) -> None:
        self.atoms = []
        if lattice is None:
            if len(kwargs) == 0:
                kwargs["cell"] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        else:
            kwargs = {}
            kwargs["cell"] = lattice.cell

        super().__init__(standardize=standardize, **kwargs)

        if atoms is not None:
            for a in atoms:
                self.add_atom(a, relative=relative)

    def __iter__(self):
        return CrystalIterator(self)

    def __contains__(self, atom: Union[Atom, str]):
        if isinstance(atom, str):
            try:
                atom = self.get_atom(atom, return_all=True)
                return True
            except ValueError:
                return False
        return atom in self.atoms

    def __len__(self):
        return self.atoms.__len__()

    def __getitem__(self, name) -> Atom:
        return self.get_atom(name, return_all=True)

    def __getattr__(self, name):
        # Fix copy/deepcopy RecursionError
        if name in ["__setstate__"]:
            raise AttributeError(name)
        try:
            atom = self.get_atom(name=name)
            return atom
        except ValueError:
            raise AttributeError(f"'Crystal' object has no attribute '{name}'")

    @property
    def lattice(self):
        r"""
        Lattice of the crystal.

        It returns an independent instance of the :py:class:`.Lattice` class.
        You can use it to play with the crystal`s lattice independently,
        but it will not affect the crystal itself.

        See :ref:`guide_crystal_lattice` for details.

        Returns
        -------
        lattice : :py:class:`.Lattice`
            Lattice of the crystal.

        Notes
        -----
        New :py:class:`.Lattice` object is created each time you call this property.
        It is created with ``standardize=False`` parameter.
        """
        return Lattice(self.cell, standardize=False)

    def add_atom(self, new_atom: Union[Atom, str] = None, relative=True, **kwargs):
        r"""
        Add atom to the crystal.

        If name and index of the ``new_atom`` are the same as of some atom of the crystal,
        then ``new_atom`` is not added.

        If index of ``new_atom`` is not defined, it is set.

        Parameters
        ----------
        new_atoms : :py:class:`.Atom` or str, optional
            New atom. All kwargs are ignored if ``new_atom`` is not ``None`` and of type :py:class:`.Atom`.
            If ``str``, then pair ``name : new_atom`` is added to ``kwargs``.
        relative : bool, default True
            Whether ``new_atom`` position is in relative coordinates.
        **kwargs
            Keyword arguments for :py:class:`.Atom` initialization.

        Raises
        ------
        TypeError
            If ``new_atom`` is not an :py:class:`.Atom`.
        ValueError
            If the atom is already present in the crystal.
        """

        if isinstance(new_atom, str):
            kwargs["name"] = new_atom
            new_atom = None

        if new_atom is None:
            new_atom = Atom(**kwargs)

        if not isinstance(new_atom, Atom):
            raise TypeError("New atom is not an Atom. " + f"Received {type(new_atom)}.")
        try:
            i = new_atom.index
        except ValueError:
            new_atom.index = len(self.atoms) + 1

        if not relative:
            new_atom.position = absolute_to_relative(self.cell, new_atom.position)

        if new_atom not in self.atoms:
            self.atoms.append(new_atom)
        else:
            raise ValueError("Atom is already in the crystal.")

    def remove_atom(self, atom: Union[Atom, str], index=None):
        r"""
        Remove atom from the crystal.

        If type(``atom``) == ``str``, then all atoms with the name ``atom`` are removed
        if ``index`` == ``None`` or only atom with the name ``atom`` and index ``index`` is removed.

        It type(``atom``) == :py:class:`.Atom`, ``index`` is ignored.

        Parameters
        ----------
        atom : :py:class:`.Atom` or str
            :py:class`.Atom` object or atom`s name.
            If name, then it has to be unique among atoms of the crystal.
        index : optional
            Index of the atom.

        Raises
        ------
        ValueError
            If no match is found.
        """

        if isinstance(atom, str):
            atoms = self.get_atom(atom, index=index, return_all=True)
        else:
            atoms = [atom]
        for atom in atoms:
            self.atoms.remove(atom)

    # Modification of position has to be avoided here.
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
            Name of the atom. In general not unique. If ``name`` contains "__", then it is split
            into ``name`` and ``index``.
        index : int, optional
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
            If no match is found or the match is not unique and ``return_all`` is ``False``.
        """

        atoms = []
        if "__" in name and index is None:
            name, index = name.split("__")
            index = int(index)

        for atom in self.atoms:
            if atom.name == name:
                if index is None or atom.index == index:
                    atoms.append(atom)

        if len(atoms) == 0:
            raise ValueError(f"No match found for name = {name}, index = {index}")
        elif len(atoms) == 1:
            if return_all:
                return atoms
            return atoms[0]
        elif not return_all:
            raise ValueError(
                f"Multiple matches found for name = {name}, index = {index}"
            )
        return atoms

    def get_atom_coordinates(
        self, atom: Union[Atom, str], R=(0, 0, 0), index=None, relative=True
    ):
        r"""
        Getter for the atom coordinates.

        Parameters
        ----------
        atom : :py:class:`.Atom` or str
            :py:class`.Atom` object or atom`s name.
            If name, then it has to be unique among atoms of the crystal.
        index : int, optional
            Index of the atom.
        R : (3,) |array_like|_, default (0, 0, 0)
            Radius vector of the unit cell for atom2 (i,j,k).
        relative : bool, default True
            Whether to return relative coordinates.

        Returns
        -------
        coordinates : 1 x 3 array
            Coordinates of atom in the cell R in absolute coordinates.
        """

        if isinstance(atom, str):
            atom = self.get_atom(atom, index=index)
        elif atom not in self.atoms:
            raise ValueError(f"There is no {atom} in the crystal.")

        rel_coordinates = np.array(R + atom.position)

        if relative:
            return rel_coordinates

        return rel_coordinates @ self.cell

    def get_vector(
        self,
        atom1: Union[Atom, str],
        atom2: Union[Atom, str],
        R=(0, 0, 0),
        index1=None,
        index2=None,
        relative=False,
    ):
        r"""
        Getter for vector from atom1 to atom2.

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
        relative : bool, default False
            Whether to return the vector relative coordinates.
        index1 : int, optional
            Index of the atom1.
        index2 : int, optional
            Index of the atom2.

        Returns
        -------
        v : (3,) :numpy:`ndarray`
            Vector from atom1 in (0,0,0) cell to atom2 in R cell.
        """

        coord1 = self.get_atom_coordinates(atom1, index=index1, relative=relative)
        coord2 = self.get_atom_coordinates(atom2, R, index=index2, relative=relative)

        return coord2 - coord1

    def get_distance(
        self,
        atom1: Union[Atom, str],
        atom2: Union[Atom, str],
        R=(0, 0, 0),
        index1=None,
        index2=None,
        relative=False,
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
        relative : bool, default False
            Whether to use relative coordinates. (Strange, but if you wish)
        index1 : int, optional
            Index of the atom1.
        index2 : int, optional
            Index of the atom2.

        Returns
        -------
        distance : floats
            Distance between atom1 in (0,0,0) cell and atom2 in R cell.
        """

        return np.linalg.norm(
            self.get_vector(
                atom1, atom2, R, index1=index1, index2=index2, relative=relative
            )
        )

    def find_primitive_cell(self):
        r"""
        Detect primitive cell.

        Before the detection of the primitive cell the corresponding bravais lattice type may not
        be correct, since it is determined with the current cell, which is not necessary primitive one.
        """
        self.cell, self.atoms = Cell.primitive(self.cell, self.atoms)

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
        for atom in self.atoms:
            try:
                tmp = atom.magmom
                magnetic_atoms.append(atom)
            except ValueError:
                pass

        if len(magnetic_atoms) == 0:
            raise ValueError("There are no magnetic atoms in the crystal.")

        magnetic_centres = np.zeros((n_t * len(magnetic_atoms), 2, 3), dtype=float)
        for a_i, atom in enumerate(magnetic_atoms):
            magnetic_centres[a_i * n_t : (a_i + 1) * n_t, 0] = np.tile(
                atom.magmom, (n_t, 1)
            )
            magnetic_centres[a_i * n_t : (a_i + 1) * n_t, 1] = (
                translations + atom.position @ self.cell
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
        for atom in self.atoms:
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
            mc1[a_i * n_t : (a_i + 1) * n_t, 1] = (
                translations + atom.position @ self.cell
            )

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
                mc2[a_i * n_t : (a_i + 1) * n_t, 1] = (
                    translations + atom.position @ self.cell
                )
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

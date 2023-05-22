from rad_tools.crystal.lattice import Lattice
from rad_tools.crystal.bravais_lattice import lattice_example
from rad_tools.crystal.atom import Atom


class Crystal:
    r"""
    Crystal class.

    Iterable over atoms. All attributes of the :py:class:`.Lattice`
    are accessible directly from the crystal or from the lattice attribute,
    i.e. the following lines are equivalent:

    .. doctest::

        >>> import rad_tools as rad
        >>> crystal = rad.Crystal()
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
        then cubic lattice is used (:math:`a = \pi`).
    atoms : list, optional
        List of :py:class:`Atom` objects.

    Attributes
    ----------
    lattice : :py:class:`.Lattice`
        Lattice of the crystal.
    atoms : list
        List of atoms of the crystal.
    """

    def __init__(self, lattice: Lattice = None, atoms=None) -> None:
        self._lattice = None
        self.atoms = []

        if lattice is None:
            lattice = lattice_example("CUB")
        self.lattice = lattice
        if atoms is not None:
            for a in atoms:
                self.add_atom(a)

    def __iter__(self):
        return CrystalIterator(self)

    def __contains__(self, atom):
        return atom in self.atoms

    def __getitem__(self, index) -> Atom:
        return self.atoms[index]

    def __getattr__(self, name):
        return getattr(self.lattice, name)

    @property
    def lattice(self):
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
        self.atoms.append(new_atom)


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

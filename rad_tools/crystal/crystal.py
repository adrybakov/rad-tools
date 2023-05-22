from rad_tools.crystal.lattice import Lattice
from rad_tools.crystal.atom import Atom


class Crystal:
    def __init__(self, lattice: Lattice, atoms=None) -> None:
        self._lattice = None
        self.atoms = []

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

from rad_tools.crystal.lattice import Lattice


class Crystal:
    def __init__(self, lattice) -> None:
        self._lattice = None

        self.lattice = lattice

    @property
    def lattice(self):
        if self._lattice is not None:
            return self._lattice

    @lattice.setter
    def lattice(self, new_lattice):
        if not isinstance(new_lattice, Lattice):
            raise TypeError(
                "New lattice is not a lattice. " + f"Received {type(new_lattice)}."
            )

r"""Atom class"""

import numpy as np


class Atom:
    r"""
    Atom.

    Parameters
    ----------
    literal : str, default X
        Name of the atom.
    position : (3,) array_like, default [0,0,0]
        Position of the atom in relative coordinates.
    spin : (3,) array_like, default None
        Classical spin vector of the atom.
    magmom : (3,) array_like, default None
        Magnetic moment of the atom.

    Attributes
    ----------
    literal : str
    position : (3,) array_like
        Position of the atom in relative coordinates.
    spin : (3,) array_like
    magmom : (3,) array_like

    """

    def __init__(
        self,
        literal="X",
        position=(
            0,
            0,
            0,
        ),
        spin=None,
        magmom=None,
    ) -> None:
        self.literal = literal
        self.position = np.array(position)
        self._spin = None
        self._magmom = None
        if spin is not None:
            self.spin = spin
        if magmom is not None:
            self.magmom = magmom

    @property
    def spin(self):
        r"""
        Classical spin vector of the atom.

        .. code-block:: python
            S_vec = [S_x, S_y, S_z]
            S = np.linalg.norm(np.array(S))
        """

        if self._spin is None:
            raise ValueError(f"Spin is not defined for the atom {self.literal}.")
        return self._spin

    @spin.setter
    def spin(self, new_spin):
        try:
            new_spin = np.array(new_spin)
        except:
            raise ValueError(f"New spin value is not array-like, new_spin = {new_spin}")
        if new_spin.shape != (3,):
            raise ValueError(
                f"New spin has to be a 3 x 1 vector, shape: {new_spin.shape}"
            )
        self._spin = new_spin

    @property
    def magmom(self):
        r"""
        Magnetic moment of the atom.

        .. code-block:: python
            magmom = [S_x, S_y, S_z]
        """

        if self._magmom is None:
            raise ValueError(
                f"Magnetic moment is not defined for the atom {self.literal}."
            )
        return self._magmom

    @magmom.setter
    def magmom(self, new_magmom):
        try:
            new_magmom = np.array(new_magmom)
        except:
            raise ValueError(
                f"New magnetic moment value is not array-like, new_magmom = {new_magmom}"
            )
        if new_magmom.shape != (3,):
            raise ValueError(
                f"New magnetic moment has to be a 3 x 1 vector, shape: {new_magmom.shape}"
            )
        self._magmom = new_magmom

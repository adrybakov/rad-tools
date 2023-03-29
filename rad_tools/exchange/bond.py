r"""
Exchange bond.
"""

from math import sqrt

import numpy as np


class Bond:
    r"""
    Exchange bond.

    If ``matrix`` is specified then ``iso``, ``aniso`` and ``dmi`` are
    ignored and derived from ``matrix``. If ``matrix`` is not specified
    then it is derived from ``iso``, ``aniso`` and ``dmi``.
    If nothing is provided exchange matrix is set to zero matrix.

    Parameters
    ----------
    iso : int or float, default None
        Value of isotropic exchange parameter.
    aniso : 3 x 3 array, None
        3 x 3 matrix of symmetric anisotropic exchange.
    dmi : 1 x 3 array, None
        Dzyaroshinsky-Moria interaction vector :math:`(D_x, D_y, D_z)`.
    matrix : 3 x 3 array, None
        Exchange matrix.

    Attributes
    ----------
    matrix : 3 x 3 array of floats
        Exchange matrix.

        .. code-block:: python

            [[J_xx, J_xy, J_xz],
             [J_yx, J_yy, J_yz],
             [J_zx, J_zy, J_zz]]

    symm_matrix : 3 x 3 array of floats
    asymm_matrix : 3 x 3 array of floats
    iso : float
    iso_matrix : 3 x 3 array of floats
    aniso : 3 x 3 array of floats
    aniso_diagonal : 1 x 3 array of floats
    aniso_diagonal_matrix : 3 x 3 array of floats
    dmi : 1 x 3 array of floats
    dmi_matrix : 3 x 3 array of floats
    dmi_module : float
    """

    def __init__(self, matrix=None, iso=None, aniso=None, dmi=None) -> None:
        self._matrix = np.zeros((3, 3), dtype=float)
        self._iso = 0.0
        self._aniso = np.zeros((3, 3), dtype=float)
        self._dmi = np.zeros(3, dtype=float)

        if matrix is None:
            self.iso = iso
            self.dmi = dmi
            self.aniso = aniso
        else:
            self.matrix = matrix

    @property
    def matrix(self):
        return self._matrix

    @matrix.setter
    def matrix(self, new_matrix):
        if new_matrix is not None:
            new_matrix = np.array(new_matrix, dtype=float)
            if new_matrix.shape != (3, 3):
                raise ValueError("Matrix shape has to be equal to (3, 3)")
        else:
            new_matrix = np.zeros((3, 3), dtype=float)
        self._matrix = new_matrix

    @property
    def symm_matrix(self):
        r"""
        Symmetric part of exchange matrix.

        .. math::

            J_symm = \dfrac{J + J^T}{2}
        """

        return (self.matrix + self.matrix.T) / 2

    @property
    def asymm_matrix(self):
        """
        Asymmetric part of exchange matrix.

        .. math::

            J_asymm = \dfrac{J - J^T}{2}
        """

        return (self.matrix - self.matrix.T) / 2

    @property
    def iso(self):
        r"""
        Value of isotropic exchange parameter.

        Derived from the exchange matrix (:math:`\mathbf{J}`) as

        .. math::
            J_{iso} = \dfrac{1}{3}Tr(\mathbf{J})
        """

        return np.trace(self.symm_matrix) / 3

    @iso.setter
    def iso(self, new_iso):
        if new_iso is not None:
            # correction is correct (see matrix update)
            new_iso = new_iso - self.iso
        else:
            new_iso = -self.iso
        self.matrix = self.matrix + (new_iso) * np.identity(3, dtype=float)

    @property
    def iso_matrix(self):
        r"""
        Isotropic part of the exchange matrix.

        Matrix form:

        .. code-block:: python

            [[J, 0, 0],
             [0, J, 0],
             [0, 0, J]]
        """

        return self.iso * np.identity(3, dtype=float)

    @property
    def aniso(self):
        r"""
        3 x 3 matrix of symmetric anisotropic exchange.

        .. code-block:: python

            [[J_xx, J_xy, J_xz],
             [J_xy, J_yy, J_yz],
             [J_xz, J_yz, J_zz]]

        Derived from the exchange matrix (:math:`\mathbf{J}`) as

        .. math::
            \mathbf{J}_{aniso} = \mathbf{J}_{symm} - \dfrac{1}{3}Tr(\mathbf{J})
            \cdot \mathbf{I}

        where :math:`\mathbf{I}` is a :math:`3\times3` identity matrix.
        """

        return self.symm_matrix - self.iso * np.identity(3, dtype=float)

    @aniso.setter
    def aniso(self, new_aniso):
        if new_aniso is not None:
            new_aniso = np.array(new_aniso, dtype=float)
            if new_aniso.shape != (3, 3):
                raise ValueError("Aniso matrix shape have " + "to be equal to (3, 3)")
            # correction is correct (see matrix update)
            new_aniso = new_aniso - self.aniso
        else:
            new_aniso = -self.aniso
        self.matrix = self.matrix + new_aniso

    @property
    def aniso_diagonal(self):
        r"""
        Diagonal part of the symmetric anisotropic exchange.

        .. code-block:: python

            [J_xx, J_yy, J_zz]
        """

        return np.diag(self.aniso)

    @property
    def aniso_diagonal_matrix(self):
        r"""
        Diagonal part of the symmetric anisotropic exchange matrix.

        .. code-block:: python

            [[J_xx, 0, 0],
             [0, J_yy, 0],
             [0, 0, J_zz]]
        """

        return np.diag(np.diag(self.aniso))

    @property
    def dmi(self):
        r"""
        Dzyaroshinsky-Moria interaction vector (Dx, Dy, Dz).

        .. code-block:: python

            [D_x, D_y, D_z]
        """

        return np.array(
            [self.asymm_matrix[1][2], self.asymm_matrix[2][0], self.asymm_matrix[0][1]],
            dtype=float,
        )

    @dmi.setter
    def dmi(self, new_dmi):
        if new_dmi is not None:
            new_dmi = np.array(new_dmi, dtype=float)
            if new_dmi.shape != (3,):
                raise ValueError(f"DMI have to be a 3 component vector. {new_dmi}")
            new_dmi = new_dmi - self.dmi
        else:
            new_dmi = -self.dmi
        dmi_matrix = np.array(
            [
                [0, new_dmi[2], -new_dmi[1]],
                [-new_dmi[2], 0, new_dmi[0]],
                [new_dmi[1], -new_dmi[0], 0],
            ],
            dtype=float,
        )
        self.matrix = self.matrix + dmi_matrix

    @property
    def dmi_matrix(self):
        r"""
        Asymmetric part of the exchange matrix.

        .. code-block:: python

            [[0, D_z, -D_y],
             [-D_z, 0, D_x],
             [D_y, -D_x, 0]]
        """

        return self.asymm_matrix

    @property
    def dmi_module(self):
        r"""
        Length of the DMI vector in the units of exchange interaction.
        """

        return sqrt(np.sum(self.dmi**2))

    def round(self, decimals=4):
        r"""Round exchange values.

        Parameters
        ----------
        accuracy : int, default 4
            number of decimals after the point.
        """

        self.matrix = np.round(self.matrix, decimals=decimals)

    def __add__(self, other):
        if isinstance(other, Bond):
            matrix = self.matrix + other.matrix
            return Bond(matrix=matrix)
        else:
            raise TypeError

    def __sub__(self, other):
        if isinstance(other, Bond):
            matrix = self.matrix - other.matrix
            return Bond(matrix=matrix)
        else:
            raise TypeError

    def __mul__(self, number):
        matrix = self.matrix * number
        return Bond(matrix=matrix)

    def __rmul__(self, number):
        return self.__mul__(number)

    def __truediv__(self, number):
        matrix = self.matrix / number
        return Bond(matrix=matrix)

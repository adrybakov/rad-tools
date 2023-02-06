import numpy as np

class Bond:
    r"""
    Exchange bond.

    If ``matrix`` is specified then ``iso``, ``aniso`` and ``dmi`` will 
    be ignored and derived from ``matrix``. If ``matrix`` is not specified 
    then it will be derived from ``iso``, ``aniso`` and ``dmi``.

    Parameters
    ----------
    iso : int or float, default None
        Value of isotropic exchange parameter in meV.
    aniso : 3 x 3 array, None
        3 x 3 matrix of symmetric anisotropic exchange in meV.
    dmi : 3 x 1 array, None
        Dzyaroshinsky-Moria interaction vector :math:`(D_x, D_y, D_z)` in meV.
    matrix : 3 x 3 array, None
        Exchange matrix in meV. If ``matrix`` is specified then ``iso`` ,
        ``aniso`` and ``dmi`` will be ignored and derived from ``matrix`` .
        If ``matrix`` is not specified then it will be derived from
        ``iso`` , ``aniso`` and ``dmi`` .

    Attributes
    ----------
    matrix : 3 x 3 array of floats
        Exchange matrix in meV. 

        Matrix form 

        .. code-block:: python

            [[J_xx, J_xy, J_xz],
             [J_yx, J_yy, J_yz],
             [J_zx, J_zy, J_zz]]

    iso : float
    aniso : 3 x 3 array of floats
    dmi : 3 x 1 array of floats
    symm_matrix : 3 x 3 array of floats
    asymm_matrix : 3 x 3 array of floats
    dmi_module : float
    dmi_vs_iso : float
    """

    def __init__(self,
                 iso=None,
                 aniso=None,
                 dmi=None,
                 matrix=None) -> None:
        self._matrix = np.zeros((3, 3), dtype=float)
        self._iso = 0.
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
        if new_matrix is None:
            new_matrix = np.zeros((3, 3), dtype=float)
        new_matrix = np.array(new_matrix)
        if new_matrix.shape != (3, 3):
            raise ValueError("Matrix shape have to be equal to (3, 3)")
        self._matrix = new_matrix

    @property
    def symm_matrix(self):
        r"""
        Symmetric part of exchange matrix.
        """

        return (self.matrix + self.matrix.T) / 2

    @property
    def asymm_matrix(self):
        """
        Asymmetric part of exchange matrix.
        """

        return (self.matrix - self.matrix.T) / 2

    @property
    def iso(self):
        r"""
        Value of isotropic exchange parameter in meV. 
        If ``iso`` is not specified then it will be 0.

        Matrix form: 

        .. code-block:: python

            [[J, 0, 0],
             [0, J, 0],
             [0, 0, J]]

        Derived from the exchange matrix (:math:`\mathbf{J}`) as 

        .. math::
            J_{iso} = \dfrac{1}{3}Tr(\mathbf{J})
        """
        return np.trace(self.symm_matrix) / 3

    @iso.setter
    def iso(self, new_iso):
        if new_iso is not None:
            new_iso = new_iso - self.iso
        else:
            new_iso = -self.iso
        self.matrix = (self.matrix +
                       (new_iso) *
                       np.identity(3, dtype=float))

    @property
    def aniso(self):
        r"""
        3 x 3 matrix of symmetric anisotropic exchange in meV. 

        Matrix form: 

        .. code-block:: python

            [[J_xx, J_xy, J_xz],
             [J_xy, J_yy, J_yz],
             [J_xz, J_yz, J_zz]]

        Derived from the exchange matrix (:math:`\mathbf{J}`) as 

        .. math::
            \mathbf{J}_{aniso} = \mathbf{J}_{symm} - \dfrac{1}{3}Tr(\mathbf{J}) \cdot \mathbf{I}

        where :math:`\mathbf{I}` is an :math:`3\times3` identity matrix.
        """
        return self.symm_matrix - self.iso * np.identity(3, dtype=float)

    @aniso.setter
    def aniso(self, new_aniso):
        if new_aniso is not None:
            new_aniso = np.array(new_aniso, dtype=float)
            if new_aniso.shape != (3, 3):
                raise ValueError("Aniso matrix shape have " +
                                 "to be equal to (3, 3)")
            new_aniso = new_aniso - self.aniso
        else:
            new_aniso = -self.aniso
        self.matrix = self.matrix + new_aniso

    @property
    def dmi(self):
        r"""
        Dzyaroshinsky-Moria interaction vector (Dx, Dy, Dz) in meV. 

        Vector form: 

        .. code-block:: python

            [D_x, D_y, D_z]

        Matrix form: 

        .. code-block:: python

            [[0, D_z, -D_y],
             [-D_z, 0, D_x],
             [D_y, -D_x, 0]]

        Derived from antisymmetric part of exchange matrix (:math:`\mathbf{J}`).
        """
        return np.array([self.asymm_matrix[1][2],
                         self.asymm_matrix[2][0],
                         self.asymm_matrix[0][1]],
                        dtype=float)

    @dmi.setter
    def dmi(self, new_dmi):
        if new_dmi is not None:
            new_dmi = np.array(new_dmi, dtype=float)
            if new_dmi.shape != (3,):
                raise ValueError(
                    f"DMI have to be a 3 component vector. {new_dmi}")
            new_dmi = new_dmi - self.dmi
        else:
            new_dmi = -self.dmi
        dmi_matrix = np.array([[0, new_dmi[2], -new_dmi[1]],
                               [-new_dmi[2], 0, new_dmi[0]],
                               [new_dmi[1], -new_dmi[0], 0]],
                              dtype=float)
        self.matrix = self.matrix + dmi_matrix

    @property
    def dmi_module(self):
        r"""
        Length of the DMI vector in th e units of exchange interaction.
        """

        return sqrt(np.sum(self.dmi**2))

    @property
    def dmi_vs_iso(self):
        r"""
        Relative strength of DMI.

        .. math::

            \dfrac{\vert\vec{D}\vert}{\vert J_{iso}\vert}
        """
        return abs(self.dmi_module/self.iso)

    def round(self, accuracy=4):
        r"""Round exchange values.

        Parameters
        ----------
        accuracy : int, default 4
            number of decimals after the point.
        """
        self.matrix = np.round(self.matrix, decimals=accuracy)

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


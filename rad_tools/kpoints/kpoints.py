from math import sqrt

import numpy as np

from rad_tools.kpoints.high_symmetry_point import HighSymmetryPoints


class KPoints(HighSymmetryPoints):
    r"""
    K-points
    """

    def __init__(self, kpoints=None, path=None) -> None:
        super().__init__(kpoints, path)

    @property
    def labels(self):
        r"""
        Labels of the high symmetry k points, ready to plot.

        Returns
        -------
        label_list : list of str
            List of high symmetry k points labels.
        """
        label_list = []
        for i in range(0, len(self.path)):
            if i != 0:
                label_list[-1] += self._PLOT_LITERALS[self.path[i][0]]
            else:
                label_list.append(self._PLOT_LITERALS[self.path[i][0]])
            for j in range(1, len(self.path[i])):
                label_list.append(self._PLOT_LITERALS[self.path[i][j]])
            if i != len(self.path) - 1:
                label_list[-1] += "$\\vert$"
        return label_list

    def labels_coordinates(self, b1=None, b2=None, b3=None):
        r"""
        Coordinates of the high symmetry k points.

        Returns
        -------
        coordinate_list : list of float
            List of high symmetry k points coordinates in reciprocal space
            if ``b1``, ``b2``, ``b3`` are specified, in relative otherwise.
        """

        coordinates_list = self.coordinates(n=1, b1=b1, b2=b2, b3=b3, linear=True)
        return np.concatenate(
            ([coordinates_list[0]], coordinates_list[1:-1:2], [coordinates_list[-1]])
        )

    def coordinates(self, n=100, b1=None, b2=None, b3=None, linear=False):
        r"""
        Coordinates of all k points in the path.

        Return coordinates in reciprocal space units if ``b1``, ``b2``, ``b3``
        are specified, otherwise return in relative to the reciprocal
        lattice vectors coordinates.

        Parameters
        ----------
        n : int
            Number of k points between two high symmetry points + 1.
        b1 : 1 x 3 array
            First reciprocal lattice vector.
        b2 : 1 x 3 array
            Second reciprocal lattice vector.
        b3 : 1 x 3 array
            Third reciprocal lattice vector.
        linear : bool, default False
            Return linear coordinates for the plot.

        Returns
        -------
        k_points : N x 3 array
            Array of N total k points in the path.

            :math:`N = i \cdot n`,
            where i - number of intervals between high symmetry k points.
        """

        point_list = []
        for i in range(0, len(self.path)):
            for j in range(1, len(self.path[i])):
                first_point = self.kpoints[self.path[i][j - 1]]
                diff = self.kpoints[self.path[i][j]] - self.kpoints[self.path[i][j - 1]]

                # Move from relative to absolute coordinates
                if b1 is not None and b2 is not None and b3 is not None:
                    first_point = (
                        first_point[0] * b1 + first_point[1] * b2 + first_point[2] * b3
                    )
                    diff = diff[0] * b1 + diff[1] * b2 + diff[2] * b3
                # Move from (kx, ky, kz) to |k|
                if linear:
                    diff = sqrt(np.sum(np.array(diff) ** 2))
                    if i == 0 and j == 1:
                        first_point = 0
                    else:
                        first_point = point_list[-1]

                for step in range(0, n + 1):
                    point_list.append(first_point + step / n * diff)

        return np.array(point_list)

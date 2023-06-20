r"""
General 3D lattice.
"""

__all__ = ["Kpoints"]


class Kpoints:
    r"""
    K-point path.

    Parameters
    ----------
    path : list of list of str
        K points path. Each subpath is a list of the high symmetry points.
    points : dict
        Dictionary of the high symmetry points.
        Coordinates are given in absolute coordinates in reciprocal space.
    labels : dict
        Dictionary of the high symmetry points labels. Has to have the same
        keys as ``points``.
    n : int
        Number of points between each pair of the high symmetry points.
    """

    def __init__(self, path, points, labels, n=100):
        self._path = path
        self._hs_points = [point for point in points]
        self_hs_coordinates = points
        self._labels = labels
        self._n = n

    @property
    def n(self):
        r"""
        Amount of points between each pair of the high symmetry points.

        Returns
        -------
        n : int
        """

        return self._n

    @n.setter
    def n(self, new_n):
        if not isinstance(new_n, int):
            raise ValueError(f"n has to be integer. Given: {new_n}")
        self._n = new_n

    @property
    def labels(self):
        r"""
        Labels of high symmetry points, ready to be plotted.

        For example for point "Gamma" it returns r"$\Gamma$".

        If there are two high symmetry points following one another in the path,
        it returns "X|Y" where X and Y are the labels of the two high symmetry points.

        Returns
        -------
        labels : list of str
            Labels, ready to be plotted. Same length as :py:attr:`.coordinates`.
        """

        labels = []
        for s_i, subpath in enumerate(self._path):
            if s_i != 0:
                labels[-1] += "|" + self._labels[subpath[0]]
            else:
                labels.append(self._labels[subpath[0]])
            for name in subpath[1:]:
                labels.append(self._labels[name])

        return labels

    @property
    def coordinates(self):
        r"""
        Flatten coordinates of the high symmetry points, ready to be plotted.

        Returns
        -------
        coordinates : list of float
            Coordinates, ready to be plotted. Same length as :py:attr:`.labels`.
        """

        coordinates = []
        for s_i, subpath in enumerate(self._path):
            if s_i == 0:
                coordinates.append()  # TODO
            for name in subpath[1:]:
                coordinates.append()  # TODO

        return coordinates

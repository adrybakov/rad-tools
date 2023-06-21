r"""
General 3D lattice.
"""

from typing import Iterable
import numpy as np

__all__ = ["Kpoints"]


class Kpoints:
    r"""
    K-point path.

    Parameters
    ----------
    points : dict
        Dictionary of the high symmetry points.
        Coordinates are given in absolute coordinates in reciprocal space.
    labels : dict
        Dictionary of the high symmetry points labels. Has to have the same
        keys as ``points``.
    path : str, optional
        K points path.
    n : int
        Number of points between each pair of the high symmetry points
        (high symmetry points excluded).
    """

    def __init__(self, points, labels, path=None, n=100):
        self._path = None
        self._hs_points = [point for point in points]
        self._hs_coordinates = dict(
            [(point, np.array(points[point])) for point in points]
        )
        self._labels = labels
        self._n = n

        if path is None:
            path = self._hs_points[0]
            for i in self._hs_points[1:]:
                path += "-" + i
        self.path = path

    @property
    def path(self):
        r"""
        K points path.

        Returns
        -------
        path : list of list of str
            K points path. Each subpath is a list of the high symmetry points.
        """

        return self._path

    @path.setter
    def path(self, new_path):
        if isinstance(new_path, str):
            tmp_path = new_path.split("|")
            new_path = []
            for i in range(len(tmp_path)):
                subpath = tmp_path[i].split("-")
                # Each subpath has to contain at least two points.
                if len(subpath) != 1:
                    new_path.append(subpath)
        elif isinstance(new_path, Iterable):
            tmp_path = new_path
            new_path = []
            for subpath in tmp_path:
                if isinstance(subpath, str) and "-" in subpath:
                    subpath = subpath.split("-")
                    # Each subpath has to contain at least two points.
                    if len(subpath) != 1:
                        new_path.append(subpath)
                elif (
                    not isinstance(subpath, str)
                    and isinstance(subpath, Iterable)
                    and len(subpath) != 1
                ):
                    new_path.append(subpath)
                else:
                    new_path = [tmp_path]
                    break
        # Check if all points are defined.
        for subpath in new_path:
            for point in subpath:
                if point not in self._hs_points:
                    message = f"Point '{point}' is not defined. Defined points are:"
                    for defined_point in self._hs_points:
                        message += f"\n  {defined_point} : {self._hs_coordinates[defined_point]}"
                    raise ValueError(message)
        self._path = new_path

    @property
    def path_string(self):
        r"""
        K points path as a string.

        Returns
        -------
        path : str
        """

        result = ""
        for s_i, subpath in enumerate(self.path):
            for i, name in enumerate(subpath):
                if i != 0:
                    result += "-"
                result += name
            if s_i != len(self.path) - 1:
                result += "|"

        return result

    @property
    def n(self):
        r"""
        Amount of points between each pair of the high symmetry points
        (high symmetry points excluded).

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
        for s_i, subpath in enumerate(self.path):
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
        coordinates : :numpy:`ndarray`
            Coordinates, ready to be plotted. Same length as :py:attr:`.labels`.
        """

        coordinates = []
        for s_i, subpath in enumerate(self.path):
            if s_i == 0:
                coordinates.append(0)
            for i, name in enumerate(subpath[1:]):
                coordinates.append(
                    np.linalg.norm(
                        self._hs_coordinates[name] - self._hs_coordinates[subpath[i]]
                    )
                    + coordinates[-1]
                )

        return np.array(coordinates)

    @property
    def points(self):
        r"""
        Coordinates of all points with n points between each pair of the high
        symmetry points (high symmetry points excluded).

        Returns
        -------
        points : (N, 3) :numpy:`ndarray`
            Coordinates of all points.
        """

        points = None
        for subpath in self.path:
            for i in range(len(subpath) - 1):
                name = subpath[i]
                next_name = subpath[i + 1]
                new_points = np.linspace(
                    self._hs_coordinates[name],
                    self._hs_coordinates[next_name],
                    self._n + 2,
                )
                if points is None:
                    points = new_points
                else:
                    points = np.concatenate((points, new_points))
        return points

    # It can not just call for points and flatten them, because it has to treat "|" as a special case.
    @property
    def flatten_points(self):
        r"""
        Flatten coordinates of all points with n points between each pair of the high
        symmetry points (high symmetry points excluded). Used to plot band structure, dispersion, etc.

        Returns
        -------
        flatten_points : (N, 3) :numpy:`ndarray`
            Flatten coordinates of all points.
        """

        flatten_points = None
        for s_i, subpath in enumerate(self.path):
            for i in range(len(subpath) - 1):
                name = subpath[i]
                next_name = subpath[i + 1]
                points = (
                    np.linspace(
                        self._hs_coordinates[name],
                        self._hs_coordinates[next_name],
                        self._n + 2,
                    )
                    - self._hs_coordinates[name]
                )
                delta = np.linalg.norm(points, axis=1)
                if s_i == 0 and i == 0:
                    flatten_points = delta
                else:
                    delta += flatten_points[-1]
                    flatten_points = np.concatenate((flatten_points, delta))
        return flatten_points

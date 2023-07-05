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
    b1 : (3,) array_like
        First reciprocal lattice vector :math:`\mathbf{b}_1`.
    b2 : (3,) array_like
        Second reciprocal lattice vector :math:`\mathbf{b}_2`.
    b3 : (3,) array_like
        Third reciprocal lattice vector :math:`\mathbf{b}_3`.
    points : dict, optional
        Dictionary of the high symmetry points.
        Coordinates are given in relative coordinates in reciprocal space.
    labels : dict, optional
        Dictionary of the high symmetry points labels. Has to have the same
        keys as ``points``.
    path : str, optional
        K points path.
    n : int
        Number of points between each pair of the high symmetry points
        (high symmetry points excluded).
    """

    def __init__(self, b1, b2, b3, points=None, labels=None, path=None, n=100) -> None:
        self.b1 = np.array(b1)
        self.b2 = np.array(b2)
        self.b3 = np.array(b3)
        self._path = None
        if points is None:
            self.hs_points = []
        else:
            self.hs_points = [point for point in points]
            self.hs_coordinates = dict(
                [(point, np.array(points[point])) for point in points]
            )
        self._n = n

        if labels is None:
            self._labels = []
        else:
            if len(labels) != len(self.hs_points):
                raise ValueError(
                    f"Amount of labels ({len(labels)}) does not match amount of points ({len(self.hs_points)})."
                )
            self._labels = [
                (point, labels[i]) for i, point in enumerate(self.hs_points)
            ]

        if path is None:
            path = ""
            for i in self.hs_points:
                if i != 0:
                    path += "-"
                path += i
        self.path = path

    def add_hs_point(self, name, coordinates, plot_label, relative=True):
        r"""
        Add high symmetry point.

        Parameters
        ----------
        name : str
            Name of the high symmetry point.
        coordinates : (3,) array_like
            Coordinates of the high symmetry point.
        plot_label : str
            Label of the high symmetry point, ready to be plotted.
        relative : bool, optional
            Whether to interpret coordinates as relative or absolute.
        """

        if name in self.hs_points:
            raise ValueError(f"Point '{name}' already defined.")
        if relative:
            cell = np.array([self.b1, self.b2, self.b3])
        else:
            cell = np.eye(3)
        self.hs_points.append(name)
        self.hs_coordinates[name] = np.array(coordinates) @ cell
        self._labels.append(plot_label)

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
                if point not in self.hs_points:
                    message = f"Point '{point}' is not defined. Defined points are:"
                    for defined_point in self.hs_points:
                        message += f"\n  {defined_point} : {self.hs_coordinates[defined_point]}"
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

    def coordinates(self, relative=False):
        r"""
        Flatten coordinates of the high symmetry points, ready to be plotted.

        Parameters
        ----------
        relative : bool, optional
            Whether to use relative coordinates instead of the absolute ones.
        Returns
        -------
        coordinates : :numpy:`ndarray`
            Coordinates, ready to be plotted. Same length as :py:attr:`.labels`.
        """

        if relative:
            cell = np.eye(3)
        else:
            cell = np.array([self.b1, self.b2, self.b3])

        coordinates = []
        for s_i, subpath in enumerate(self.path):
            if s_i == 0:
                coordinates.append(0)
            for i, name in enumerate(subpath[1:]):
                coordinates.append(
                    np.linalg.norm(
                        self.hs_coordinates[name] @ cell
                        - self.hs_coordinates[subpath[i]] @ cell
                    )
                    + coordinates[-1]
                )

        return np.array(coordinates)

    def points(self, relative=False):
        r"""
        Coordinates of all points with n points between each pair of the high
        symmetry points (high symmetry points excluded).

        Parameters
        ----------
        relative : bool, optional
            Whether to use relative coordinates instead of the absolute ones.

        Returns
        -------
        points : (N, 3) :numpy:`ndarray`
            Coordinates of all points.
        """

        if relative:
            cell = np.eye(3)
        else:
            cell = np.array([self.b1, self.b2, self.b3])

        points = None
        for subpath in self.path:
            for i in range(len(subpath) - 1):
                name = subpath[i]
                next_name = subpath[i + 1]
                new_points = np.linspace(
                    self.hs_coordinates[name] @ cell,
                    self.hs_coordinates[next_name] @ cell,
                    self._n + 2,
                )
                if points is None:
                    points = new_points
                else:
                    points = np.concatenate((points, new_points))
        return points

    # It can not just call for points and flatten them, because it has to treat "|" as a special case.
    def flatten_points(self, relative=False):
        r"""
        Flatten coordinates of all points with n points between each pair of the high
        symmetry points (high symmetry points excluded). Used to plot band structure, dispersion, etc.

        Parameters
        ----------
        relative : bool, optional
            Whether to use relative coordinates instead of the absolute ones.

        Returns
        -------
        flatten_points : (N, 3) :numpy:`ndarray`
            Flatten coordinates of all points.
        """

        if relative:
            cell = np.eye(3)
        else:
            cell = np.array([self.b1, self.b2, self.b3])

        flatten_points = None
        for s_i, subpath in enumerate(self.path):
            for i in range(len(subpath) - 1):
                name = subpath[i]
                next_name = subpath[i + 1]
                points = (
                    np.linspace(
                        self.hs_coordinates[name] @ cell,
                        self.hs_coordinates[next_name] @ cell,
                        self._n + 2,
                    )
                    - self.hs_coordinates[name] @ cell
                )
                delta = np.linalg.norm(points, axis=1)
                if s_i == 0 and i == 0:
                    flatten_points = delta
                else:
                    delta += flatten_points[-1]
                    flatten_points = np.concatenate((flatten_points, delta))
        return flatten_points

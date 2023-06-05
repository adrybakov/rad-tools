from math import cos, sin, tan

import numpy as np


class HighSymmetryPoints:
    r"""
    Generator of the high symmetry k-points.

    Parameters
    ----------
    points : dict
        Dictionary of the high symmetry k-points and coordinate
        (fractions of the reciprocal vectors).

        .. code-block:: python

            points = {"name" : [xb1, xb2, xb3]}
    """

    def __init__(self, points=None, path=None) -> None:
        self._points = points
        self._path = path

    @property
    def points(self):
        r"""
        Dictionary of K-points.

        .. code-block:: python

            points = {Point1 : [k1, k2, k3], ...}

        where ``k1``, ``k2``, ``k3`` are
        the relative coordinates of the ``Point1``.
        """

        if self._points is None:
            return {}
        else:
            return self._points

    def add_kpoint(self, name, coordinates, plot_name=None):
        r"""
        Add one kpoint to the set of High symmetry points.

        Parameters
        ----------
        name : str
            Name of the points, which is used by the code as a key.
        coordinates : 1x3 array
            Relative coordinates of the kpoint.
        plot_name : str
            name of the points to be display on the plot.
            Equal to the ``name`` by default.
        """

        coordinates = np.array(coordinates)
        if coordinates.shape != (3,):
            raise ValueError
        if self._points is None:
            self._points = {}
        self._points[name] = coordinates

        if plot_name is None and name not in self._PLOT_LITERALS:
            self._PLOT_LITERALS[name] = name
        elif plot_name is not None:
            self._PLOT_LITERALS[name] = plot_name

    @property
    def path(self):
        r"""
        Path of K-points.

        Always is a list of lists.

        .. code-block:: python

            path = [subpath1, subpath2, ...]
            subpath1 = [hs_kpoint1, hs_kpoint2, ...]
        """

        if self._path is None:
            return []
        else:
            return self._path

    @path.setter
    def path(self, new_path):
        if new_path is None:
            self._path = []

        if isinstance(new_path, str):
            if "|-" in new_path or "-|" in new_path:
                raise ValueError("Check the format of the new path")
            subpaths = new_path.split("|")
            if subpaths[-1] == "":
                subpaths = subpaths[0:-1]
            if subpaths[0] == "":
                subpaths = subpaths[1:]

            for i in range(0, len(subpaths)):
                subpaths[i] = subpaths[i].split("-")
                for j in range(0, len(subpaths[i])):
                    if subpaths[i][j] not in self.points:
                        raise ValueError(
                            f"Provided  high symmetry "
                            + f"k-point {subpaths[i][j]} is not known."
                        )
            self._path = subpaths

        if isinstance(new_path, list):
            if len(new_path) == 0:
                self._path = []
            elif isinstance(new_path[0], list):
                for i in range(0, len(new_path)):
                    for j in range(0, len(new_path[i])):
                        if new_path[i][j] not in self.points:
                            raise ValueError(
                                f"Provided  high symmetry "
                                + f"k-point {new_path[i][j]} is not known."
                            )
                if len(new_path) == 1 and len(new_path[0]) == 0:
                    self._path = []
                else:
                    self._path = new_path
            elif isinstance(new_path[0], str):
                for i in range(0, len(new_path)):
                    for j in range(0, len(new_path[i])):
                        if new_path[i][j] not in self.points:
                            raise ValueError(
                                f"Provided  high symmetry "
                                + f"k-point {new_path[i][j]} is not known."
                            )
                self._path = [new_path]

            else:
                raise ValueError(f"Path format is not supported: {new_path}")

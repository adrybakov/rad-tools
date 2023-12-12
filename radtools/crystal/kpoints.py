# RAD-tools - program for spin Hamiltonian and magnons.
# Copyright (C) 2022-2023  Andrey Rybakov
#
# e-mail: anry@uv.es, web: adrybakov.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

r"""
General 3D lattice.
"""

from typing import Iterable
from collections import Counter
import numpy as np
from scipy.interpolate import griddata
from scipy.spatial.transform import Rotation
from radtools.geometry import absolute_to_relative

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
    coordinates : list, optional
        Coordinates are given in relative coordinates in reciprocal space.
    names: list, optional
        Names of the high symmetry points. Used for programming, not for plotting.
    labels : list, optional
        Dictionary of the high symmetry points labels for plotting.
        Has to have the same length as ``coordinates``.
    path : str, optional
        K points path.
    n : int
        Number of points between each pair of the high symmetry points
        (high symmetry points excluded).

    Attributes
    ----------
    b1 : (3,) :numpy:`ndarray`
        First reciprocal lattice vector :math:`\mathbf{b}_1`.
    b2 : (3,) :numpy:`ndarray`
        Second reciprocal lattice vector :math:`\mathbf{b}_2`.
    b3 : (3,) :numpy:`ndarray`
        Third reciprocal lattice vector :math:`\mathbf{b}_3`.
    hs_names : list
        Names of the high symmetry points. Used for programming, not for plotting.
    hs_coordinates : dict
        Dictionary of the high symmetry points coordinates.

        .. code-block:: python

            {"name": [k_a, k_b, k_c], ... }

    hs_labels : dict
        Dictionary of the high symmetry points labels for plotting.

        .. code-block:: python

            {"name": "label", ... }
    """

    def __init__(
        self, b1, b2, b3, coordinates=None, names=None, labels=None, path=None, n=100
    ) -> None:
        self.b1 = np.array(b1)
        self.b2 = np.array(b2)
        self.b3 = np.array(b3)

        if coordinates is None:
            coordinates = []

        # Fill names and labels with defaults
        if names is None:
            names = [f"K{i+1}" for i in range(len(coordinates))]
            if labels is None:
                labels = [f"K$_{i+1}$" for i in range(len(coordinates))]
        if labels is None:
            labels = [name for name in names]
        else:
            if len(labels) != len(coordinates):
                raise ValueError(
                    f"Amount of labels ({len(labels)}) does not match amount of points ({len(coordinates)})."
                )

        # Define high symmetry points attributes
        self.hs_coordinates = dict(
            [(names[i], np.array(coordinates[i])) for i in range(len(coordinates))]
        )
        self.hs_labels = dict([(names[i], labels[i]) for i in range(len(coordinates))])
        self.hs_names = names

        self._n = n

        self._path = None
        if path is None:
            if len(self.hs_names) > 0:
                path = [self.hs_names[0]]
            else:
                path = []
            for point in self.hs_names[1:]:
                path.append(f"-{point}")
            path = "".join(path)
        self.path = path

    def add_hs_point(self, name, coordinates, label, relative=True):
        r"""
        Add high symmetry point.

        Parameters
        ----------
        name : str
            Name of the high symmetry point.
        coordinates : (3,) array_like
            Coordinates of the high symmetry point.
        label : str
            Label of the high symmetry point, ready to be plotted.
        relative : bool, optional
            Whether to interpret coordinates as relative or absolute.
        """

        if name in self.hs_names:
            raise ValueError(f"Point '{name}' already defined.")

        if not relative:
            coordinates = absolute_to_relative(
                np.array([self.b1, self.b2, self.b3]), coordinates
            )

        self.hs_names.append(name)
        self.hs_coordinates[name] = np.array(coordinates)
        self.hs_labels[name] = label

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
                if point not in self.hs_names:
                    message = f"Point '{point}' is not defined. Defined points are:"
                    for defined_name in self.hs_names:
                        message += (
                            f"\n  {defined_name} : {self.hs_coordinates[defined_name]}"
                        )
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
            raise ValueError(
                f"n has to be integer. Given: {new_n}, type = {type(new_n)}"
            )
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
                labels[-1] += "|" + self.hs_labels[subpath[0]]
            else:
                labels.append(self.hs_labels[subpath[0]])
            for name in subpath[1:]:
                labels.append(self.hs_labels[name])

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

# function checking if the k points of the list k_points_subgrid are inside the brillouin zone
def check_inside_brillouin_zone(
    k_points_subgrid,
    reciprocal_vectors_2d,
    brillouin_primitive_vectors,
    plane_2d,
):
    matrix_crystal_to_cartesian = np.zeros((3, 3))
    matrix_cartesian_to_crystal = np.zeros((3, 3))
    k_points_subgrid_tmp = np.zeros((len(k_points_subgrid[:, 0]), 3))
    matrix_crystal_to_cartesian[:, 0] = brillouin_primitive_vectors[0, :]
    matrix_crystal_to_cartesian[:, 1] = brillouin_primitive_vectors[1, :]
    matrix_crystal_to_cartesian[:, 2] = brillouin_primitive_vectors[2, :]
    matrix_cartesian_to_crystal = np.linalg.inv(np.matrix(matrix_crystal_to_cartesian))
    # writing k points in crystal coordinates
    for r in range(len(k_points_subgrid[:, 0])):
        for s in range(0, 3):
            for t in range(0, 3):
                k_points_subgrid_tmp[r, s] = (
                    k_points_subgrid_tmp[r, s]
                    + matrix_cartesian_to_crystal[s, t] * k_points_subgrid[r, t]
                )
        count = 0
        # if the crystal coordinate of a k point is beyond +-1, it means that the k point is outside the brillouin zone
        for l in range(0, 3):
            # considering only those primitive vectors defining the 2D plane chosen
            if plane_2d[l] != 0:
                if abs(k_points_subgrid_tmp[r, l]) >= 1:
                    for s in range(0, 3):
                        # translating the respective cartesian coordinates
                        k_points_subgrid[r, s] = (
                            k_points_subgrid[r, s]
                            - int(k_points_subgrid_tmp[r, l])
                            * reciprocal_vectors_2d[count][s]
                        )
            # progressing on the 2D plane primitive vectors
            count = count + 1
    # returning the properly translated k points
    return k_points_subgrid

def symmetry_analysis(
    k_origin,
    k_origin_weight,
    k_points_subgrid,
    symmetries,
    threshold_k_grid,
):
    r"""

    #TODO Short description

    The symmetry analysis is applied on a subset of k points ``k_points_subgrid``,
    the origin of the subsystem is considered as well ``k_origin``
    in our case the origin is the point to which the refinement procedure is applied
    from the symmetry analysis is clear the distribution of the origin weight
    ``k_origin_weight`` between the different k points of the subset

    Parameters
    ----------
    k_origin : (3,) |array-like|_
        #TODO description
    k_origin_weight : float
        #TODO description
    k_points_subgrid : (N, 3) |array-like|_
        #TODO description
    symmetries : list of lists
        #TODO better description
        A symmetry is a list of 3 elements
        (the versor is the axis of rotation, while the modulus is the angle).
    threshold_k_grid : float
        Threshold to recognize a symmetry.

    Returns
    -------
    k_points_subgrid_weight_tmp: (N,) :numpy:`ndarray`
        Distribution of the origin k point weight between the different k points
    """

    k_points_subgrid = np.array(k_points_subgrid)
    k_origin = np.array(k_origin)

    # It is used 10 times, so it is better to calculate it once
    N = k_points_subgrid.shape[0]  # former len(k_points_subgrid[:, 0])

    k_points_subgrid_weight_tmp = (k_origin_weight / 4) * np.ones(N, dtype=float)

    # If there are no symmetry operations,
    # than the weight on each subgrid k point is exactly 1/4 of the original weight
    if symmetries[0][0] == symmetries[0][1] == symmetries[0][2] == 0:
        return k_points_subgrid_weight_tmp

    # Defining a mpiatrix to save the degeneracies, after each symmetry operation
    check_degeneracy = np.zeros((N, N), dtype=bool)
    # Saving number of degeneracies detected
    no_degeneracy = True
    # Each k point is compared to each other k point after the symmetry operation has been applied
    for i in range(0, N - 1):
        for j in range(i + 1, N):
            # For each pair considering all the symmetry operations
            for rotvec in symmetries:
                k_point_transformed = (
                    Rotation.from_rotvec(rotvec).as_matrix()
                    @ (k_points_subgrid[j] - k_origin)
                    + k_origin
                )
                if np.allclose(
                    k_points_subgrid[i], k_point_transformed, atol=threshold_k_grid
                ):
                    check_degeneracy[i][j] = True
                    no_degeneracy = False
                    break
    # If no degeneracy is detected, the common result is given
    if no_degeneracy:
        return k_points_subgrid_weight_tmp
    else:
        # If degeneracy is detected, of the degenerate points only one is chosen
        list = dict([(str(i), [i]) for i in range(N)])
        # Assuming all k points as not-degenerate, and creating a dictionary
        # where to each eigenspace is associated an eigenvector(k_point)
        for l in range(0, N - 1):
            for j in range(l + 1, N):
                if len(list) != 1:
                    # Reading the degeneracy matrix, the eigenspaces are properly populated,
                    # in particular degenerate k points are put in the same eigenspace
                    if check_degeneracy[l][j] == True:
                        for key, values in list.items():
                            for value in values:
                                if value == l:
                                    positionl = key
                                if value == j:
                                    positionj = key
                        list_new = {}
                        # New list is created with the proper eigenspaces, but at the same
                        # time is compared with the old one, to check if there are old
                        # degeneracies to take into account
                        if positionl != positionj:
                            count = 0
                            for key, values in list.items():
                                if key == positionl:
                                    list_new[str(count)] = values + list[positionj]
                                    count = count + 1
                                else:
                                    if key != positionj:
                                        list_new[str(count)] = values
                                        count = count + 1
                        # The new list is save to be used as an old one
                        list = {}
                        for key, values in list_new.items():
                            list[str(key)] = list_new[str(key)]
        # The chosen k points
        for key, values in list_new.items():
            list_new[str(key)] = sorted(set(list_new[str(key)]))
        if len(list_new) != 0:
            k_weight_tmp = k_origin_weight / len(list_new)
        else:
            k_weight_tmp = 0.0
        for key, values in list_new.items():
            if len(values) == 1:
                k_points_subgrid_weight_tmp[values] = k_weight_tmp
            else:
                count = 0
                # Chosing of the degenerate points only the first one
                for value in values:
                    if count == 0:
                        k_points_subgrid_weight_tmp[value] = k_weight_tmp
                        count = count + 1
                    else:
                        k_points_subgrid_weight_tmp[value] = 0
        # Returning the matrix with the respective weights
        return k_points_subgrid_weight_tmp


def local_refinment(
    refined_grid,
    reciprocal_vectors_2d,
    not_refined_grid,
    refinment_spacing,
    refinment_iteration,
    symmetry,
    threshold_k_grid,
    brillouin_primitive_vectors,
    plane_2d,
):
    r"""
    Starting from a set of k points (not_refined_grid) with certain weight, to each k point iteratively (refinment_iteration) associates a subset of k points,
    applying to each step a symmetry analysis (symmetry), and consequently distributing the initial weight of the starting k point to the subset k points
    The refinment subgrid is built on the plane pointed out as an input (reciprocal_vectors_2d) [generalization to 3D is quite straightforward]

    Parameters
    ----------
    refined_grid empty list []
    reciprocal_vectors_2d : (2,3) |array|
    not_refined_grid: (:,4) |array| the first 3 columns are the coordinates, while the 4th column is the weight
    refinment_spacing: |double| initial subgrid dimension, half of the preceding subgrid dimension is considered at each iteration
    refinment_iteration: |int| number of refinment iterations considered
    symmetry: |list of lists| a symmetry is a list of 3 elements (the versor is the axis of rotation, while the modulus is the angle)
    threshold_k_grid: |double| therhold to recognize a symmetry
    """
    if not_refined_grid.ndim == 1:
        length_not_refined_grid = int(len(not_refined_grid) / 4)
    else:
        length_not_refined_grid = int(len(not_refined_grid))
    if refinment_iteration == 0:
        if length_not_refined_grid > 1:
            for s in range(length_not_refined_grid):
                if not_refined_grid[s, 3] != 0:
                    refined_grid.extend(not_refined_grid[s, :])
        else:
            if not_refined_grid[3] != 0:
                refined_grid.extend(not_refined_grid)
    else:
        iter = 0
        while iter != length_not_refined_grid:
            if length_not_refined_grid == 1:
                k_tmp = not_refined_grid[:3]
                k_weight_tmp = not_refined_grid[3]
                iter = length_not_refined_grid
            else:
                k_tmp = not_refined_grid[iter, :3]
                k_weight_tmp = not_refined_grid[iter, 3]
                iter = iter + 1
            if k_weight_tmp != 0:
                k_points_subgrid = np.zeros((4, 4))
                k_points_subgrid[0, :3] = (
                    k_tmp + reciprocal_vectors_2d[0, :] * refinment_spacing
                )
                k_points_subgrid[1, :3] = (
                    k_tmp - reciprocal_vectors_2d[0, :] * refinment_spacing
                )
                k_points_subgrid[2, :3] = (
                    k_tmp + reciprocal_vectors_2d[1, :] * refinment_spacing
                )
                k_points_subgrid[3, :3] = (
                    k_tmp - reciprocal_vectors_2d[1, :] * refinment_spacing
                )
                k_points_subgrid[:, 3] = symmetry_analysis(
                    k_tmp,
                    k_weight_tmp,
                    k_points_subgrid[:, :3],
                    symmetry,
                    threshold_k_grid,
                )
                k_points_subgrid[:, :3] = check_inside_brillouin_zone(
                    k_points_subgrid[:, :3],
                    reciprocal_vectors_2d,
                    brillouin_primitive_vectors,
                    plane_2d,
                )
                for i in range(0, 4):
                    new_refinment_spacing = refinment_spacing / 2
                    if k_points_subgrid[i, 3] != 0:
                        local_refinment(
                            refined_grid,
                            reciprocal_vectors_2d,
                            k_points_subgrid[i, :],
                            new_refinment_spacing,
                            refinment_iteration - 1,
                            symmetry,
                            threshold_k_grid,
                            brillouin_primitive_vectors,
                            plane_2d,
                        )

def dynamical_refinment(
    k_point_to_refine,
    reciprocal_vectors_2d,
    refinment_spacing,
    refinment_iteration,
    symmetry,
    threshold_k_grid,
    brillouin_primitive_vectors,
    plane_2d,
):
    r"""
    In case to a point is already associated a subset of k points, but we are interested in a new refinment starting from the old one

    Parameters
    ----------
    k_point_to_refine : (:,4) |array| subgrid of k points and respective weights
    reciprocal_vectors_2d : (2,3) |array|
    refinment_spacing: |double| initial subgrid dimension, half of the preceding subgrid dimension is considered at each iteration
    refinment_iteration: |int| number of refinment iterations considered
    symmetry: |list of lists| a symmetry is a list of 3 elements (the versor is the axis of rotation, while the modulus is the angle)
    threshold_k_grid: |double| therhold to recognize a symmetry
    """
    refined_grid_tmp = []
    local_refinment(
        refined_grid_tmp,
        reciprocal_vectors_2d,
        k_point_to_refine,
        refinment_spacing,
        refinment_iteration,
        symmetry,
        threshold_k_grid,
        brillouin_primitive_vectors,
        plane_2d,
    )
    refined_grid_tmp = np.reshape(refined_grid_tmp, (int(len(refined_grid_tmp) / 4), 4))
    for s in range(len(k_point_to_refine)):
        k_point_to_refine = np.delete(k_point_to_refine, 0, axis=0)
    k_point_to_refine = np.append(k_point_to_refine, refined_grid_tmp)
    k_point_to_refine = np.reshape(
        k_point_to_refine, (int(len(k_point_to_refine) / 4), 4)
    )
    return k_point_to_refine


def k_points_grid_2d_refinment_and_symmetry(
    brillouin_primitive_vectors,
    plane_2d,
    grid_spacing,
    shift_in_plane,
    shift_in_space,
    symmetry,
    refinment_spacing,
    refinment_iteration,
    threshold_k_grid,
):
    r"""
    2d k points grid
    Parameters
    ----------
    plane_2d : (1,3) |array_like|_  ex. [0,1,1] the plane is the b2xb3 in reciprocal space where the bi are the primitive vectors
    grid_spacing: |double| spacing of the largest k grid
    shift_in_plane : (1,2) |array_like|_  shift of the reciprocal plane_2d with respcet to the reciprocal(crystal) coordinates
    shift_in_space : (1,3) |array_like|_  shift of the reciprocal plane_2d with respcet to the cartesian coordinates
    symmetry : list [symmetry1,symmetry2...] each symmetry is a list of three elements: the versor is associated to the axis of rotation, while the modulus is the angle of rotation
    refinment_spacing : |double| first refinment iteration dimension, the next refinment itarations are multiple of this
    refinment_iteration : |int| number of refinment iterations
    threshold_k_grid: |double| threshold for the symmetry checking
    Return
    -------------
    normalized_chosen_reciprocal_plane: (2,3) |matrix| versor of the chosen plane
    k0,k1 (1,1) |int,int| the number of k points along the 2 direction for the largest k grid
    total_number_k_points |int|, total number of k points, taking into account also the refinment k points
    k_points_grid_2d k0xk1 grid |k_point class| the largest k grid, where each point is associated to a refinment grid
    """
    count = 0
    chosen_reciprocal_plane = np.zeros((2, 3))
    for i in range(0, 3):
        if plane_2d[i] != 0 and count < 2:
            if i == 0:
                chosen_reciprocal_plane[count] = brillouin_primitive_vectors[0]
            elif i == 1:
                chosen_reciprocal_plane[count] = brillouin_primitive_vectors[1]
            else:
                chosen_reciprocal_plane[count] = brillouin_primitive_vectors[2]
            count = count + 1
    # np.dot(v1,v2) same as v1 @ v2
    k0 = int(
        np.dot(chosen_reciprocal_plane[0], chosen_reciprocal_plane[0]) / grid_spacing
    )
    k1 = int(
        np.dot(chosen_reciprocal_plane[1], chosen_reciprocal_plane[1]) / grid_spacing
    )
    weight = 1 / (k0 * k1)
    normalized_chosen_reciprocal_plane = np.zeros((2, 3))
    # Try to use numpy here
    for i in range(0, len(chosen_reciprocal_plane[:, 0])):
        normalized_chosen_reciprocal_plane[i] = chosen_reciprocal_plane[i] / np.dot(
            chosen_reciprocal_plane[i], chosen_reciprocal_plane[i]
        )
    k_points_grid_2d = np.zeros((k0, k1), dtype=object)
    k_points_grid_2d_tmp = np.zeros(4)
   ## with open("k_points.txt", "w") as file:
    for i in range(0, k0):
        for j in range(0, k1):
            #print(
            #    "k grid: {}% and {}%".format(int(i / k0 * 100), int(j / k1 * 100)),
            #    end="\r",
            #    flush=True,
            #)
            k_points_grid_2d_tmp[:3] = ((i + shift_in_plane[0]) / k0) * (
                shift_in_space + chosen_reciprocal_plane[0]
            ) + ((j + shift_in_plane[1]) / k1) * (
                shift_in_space + chosen_reciprocal_plane[1]
            )
            k_points_grid_2d_tmp[3] = weight
            refined_grid_tmp = []
            local_refinment(
                refined_grid_tmp,
                normalized_chosen_reciprocal_plane,
                k_points_grid_2d_tmp,
                refinment_spacing,
                refinment_iteration,
                symmetry,
                threshold_k_grid,
                brillouin_primitive_vectors,
                plane_2d,
            )
            refined_grid_tmp = np.reshape(
                refined_grid_tmp, (int(len(refined_grid_tmp) / 4), 4)
            )
            k_points_grid_2d[i][j] = np.zeros((len(refined_grid_tmp), 4))
            k_points_grid_2d[i][j] = refined_grid_tmp
  ##              file.write(f"{k_points_grid_2d[i][j]}  \n")
    ##the k points are saved in the file k_points.txt
    ##example of dynamical refinment
    ##print(k_points_grid_2d[0][0])
    ##k_points_grid_2d[0][0]=dynamical_refinment(
    ##    k_points_grid_2d[0][0],
    ##    chosen_reciprocal_plane,
    ##    refinment_spacing,
    ##    refinment_iteration,
    ##    symmetry,threshold_k_grid
    ##)
    ##print(k_points_grid_2d[0][0])
    return (normalized_chosen_reciprocal_plane, k0, k1, k_points_grid_2d)

def mapping_to_square_grid(
        k_points_list_with_weights,
        normalized_chosen_reciprocal_plane,
        epsilon
):
    r"""
    mapping of the 2d k point list (generated using known primitive vectors) into a square k point list, i.e. naming the position of the single k points
    1) first the origin is defined (as the leftest and lowest point)
    2) to each point is associated a k vector with respect to the origin
    2b) it is considered the pair of vectors producing the larger angle
    2) each k vector is projected along the two found vectors
    3) the bording points are selected and correspondingly enumerated 
    1) the points (1),(2), (3) and (4) are repeted
    ACTHUNG!!! be carefull in the choice of epsilon, otherwise you risk not to find any point
    there is the possibility that on a border there would be less than two points, this is evaluated a posteriori
    Parameters
    ----------
    k_points_list_with_weights : (2,4) |array_like|_  ex. the list contain the 3 xyz coordinates and a weight
    """

    ###defining a recursive function to define the borders
    def bording_definition(origin_indices,k_points_list,k_points_border,normalized_chosen_reciprocal_plane,epsilon):
        r"""
        parameters
        origin_indices: (1,1) |array|_ indices associated with the origin
        k_points_list: (n,3) |array|_ list of k points
        brillouin_primitive_vectors: (2,3) |array|_ brillouin primitive vectors defining a 2D plane in the Brillouin zone
        
        (iterative) return, until k=n
        k_points_list: (n-k,3) |array|_ list of k points not part of the border
        k_points_border: (k,5) |array|_ list of k points part of the border, to each k point are associated the 3 kxkykz coordinates and the two indices, defined with respect to the origin indices
        """
        if len(k_points_list)==0:
            return 0
        
        ###assuming at the beginning as origin (0,0,0)
        ###finding the real origin (the lowest and leftest point)
        number_k_points=len(k_points_list)
        k_vectors=np.zeros((number_k_points,3))
        k_vectors_projections=np.zeros((number_k_points,2))
        ###the k vectors have a weight in the 3rd position
        for i in range(0,number_k_points):
            k_vectors[i]=k_points_list[i][:3]
            for j in range(0,2):
                k_vectors_projections[i][j]=np.dot(k_vectors[i],brillouin_primitive_vectors[j])
        ###print(k_vectors_projections)
        min_value=k_vectors_projections[0][:]
        position_min_value=0
        for i in range(0,number_k_points):
            if min_value[0] > k_vectors_projections[i][0] and min_value[1] > k_vectors_projections[i][1]:
                position_min_value=i
                min_value=k_vectors_projections[i][:]
        origin=k_points_list[position_min_value][:]
        
        ###finding the maximal angle pair (the parallelogram containing all the points)
        ###recentering the vectors with respect to the real origin
        for i in range(0,number_k_points):
            if i!=position_min_value:
                k_vectors[i]=k_points_list[i][:3]-origin[:3]
                norm=np.dot(k_vectors[i],k_vectors[i])
                ###normalizing the vectors
                if norm !=0:
                    k_vectors[i]=k_vectors[i]/norm
        min_sin_angle=1
        position_min_sin_angle=[0,0]
        for i in range(0,number_k_points):
            for j in range(i+1,number_k_points):
                if i!=position_min_value:
                    sin_angle=np.dot(k_vectors[i],k_vectors[j])
                    if sin_angle<min_sin_angle:
                        min_sin_angle=sin_angle
                        position_min_sin_angle=[i,j]
        ###saving the two found axis in a variable (a1,a2)
        new_axis_vectors=np.zeros((2,3))
        for i in range(0,2):
            new_axis_vectors[i]=k_vectors[position_min_sin_angle[i]]

        epsilon_tmp=epsilon
        ###finding the k points around the two found axis (the epsilon used here is the same of the k point grid building); the control flag is used to be sure the found k points are at leat two (the two of the axis..)
        control_flag=0
        while control_flag==0:
            ####defining every point as a vector with respect to the found origin and differentiating between points along a1 and a2 and not
            k_points_not_along_a=[]
            k_points_along_a=[[],[]]
            for i in range(0,number_k_points):
                if i!=position_min_value:
                    ####here we need the not normilzed vectors
                    k_vectors[i]=k_points_list[i][:3]-origin[:3]
                    ###flag to check that the projection is not null on both the primitive reciprocal vectors
                    count_not_ptojected=0
                    count_projected=0
                    for j in range(0,2):
                        k_vectors_projections[i][j]=np.dot(k_vectors[i],new_axis_vectors[j])
                        ####grouping all elements with one of the two projections equal to zero (those at the border)
                        if abs(k_vectors_projections[i][j])<=epsilon_tmp and count_projected!=1:
                            k_points_along_a[1-j].append([i,k_vectors_projections[i][1-j]])
                            count_projected=count_projected+1
                        else:
                            count_not_ptojected=count_not_ptojected+1
                            if count_not_ptojected == 2 and count_projected==0:
                                k_points_not_along_a.append(k_points_list[i])
            number_k_points_along_a=np.zeros(2)
            for i in range(0,2):
                number_k_points_along_a[i]=int(len(k_points_along_a[i]))
            if number_k_points_along_a[0]<1 or number_k_points_along_a[1]<1:
                ###going back in the refinment procedure
                epsilon_tmp=epsilon_tmp*2
            else:
                control_flag=1

        ##reshaping a0 and a1 lists, and trasforming them to arrays
        matrix_k_points_along_a={}
        for i in range(0,2):
            matrix_k_points_along_a[i]=np.zeros((int(number_k_points_along_a[i]),2))
        for i in range(0,2):
            for j in range(int(number_k_points_along_a[i])):
                for r in range(0,2):
                    matrix_k_points_along_a[i][j][r]=k_points_along_a[i][j][r]
        
        ###addding the origin to the border_list
        k_points_border.append([k_points_list[position_min_value],origin_indices])
        ###ordering to border points, and giving them proper indices
        for i in range(0,2):
            matrix_k_points_along_a[i]=np.reshape(sorted(matrix_k_points_along_a[i],key=lambda x: x[1]),(int(number_k_points_along_a[i]),2))
            for j in range(int(number_k_points_along_a[i])):
                position=int(matrix_k_points_along_a[i][j][0])
                if i == 0:
                    indices=[origin_indices[i],j+1]
                else:
                    indices=[j+1,origin_indices[i]]
                k_points_border.append([k_points_list[position],indices])
        
        #print(k_points_border)
        #increasing the index of the origin
        new_origin_indices=[origin_indices[0]+1,origin_indices[1]+1]
        bording_definition(new_origin_indices, k_points_not_along_a, k_points_border, normalized_chosen_reciprocal_plane,epsilon)
    
    print("START OF THE COUNTOR-CUT")
    ###applying bording definition in an iterative way
    k_points_border=[]
    INFO=bording_definition([0,0],k_points_list_with_weights,k_points_border,normalized_chosen_reciprocal_plane,epsilon)
    print("END OF THE COUNTOR-CUT")
    k_points_border_matrix=np.zeros((len(k_points_border),2),dtype=list)
    for i in range(len(k_points_border)):
        k_points_border_matrix[i,0]=k_points_border[i][0]
        k_points_border_matrix[i,1]=k_points_border[i][1]
   
    ####now the list has to be converted into a matrix
    ####the rows of the matrix can have different number of columns
    ####building a square matrix, putting the missing spots to None
    def create_square_matrix(list_of_k_points):
        print(list_of_k_points)
        print(list_of_k_points[:,1])
        #Determine the matrix size based on the maximum indices
        max_pair=max(list_of_k_points[:,1])
        maxi=max_pair[0]
        maxj=max_pair[1]
        print(maxi,maxj)
        #Initialize an empty square matrix with lists containing four values (3 coordinatex kxkykz and weight)
        square_matrix = np.((maxi,maxj),dtype=list)
                                 
        [None,None,None,None])
        print(square_matrix)
        
        #Fill the matrix with specified values at specified indices
        for pair in range(0,len(list_of_k_points[:,1])):
            pairi=list_of_k_points[pair,1][0]
            pairj=list_of_k_points[pair,1][1]
            #saving coordinates and weight
            square_matrix[pairi][pairj] = list_of_k_points[pair][0]  
        return square_matrix
    
    square_grid=create_square_matrix(k_points_border_matrix)
    print(square_grid)
    #####the missing points in the square grid are filled with interpolated values
    ##def fill_matrix_with_interpolation(matrix):
    ##    # Find the indices of elements that are not arrays of 4 None elements
    ##    known_indices = np.argwhere(np.array([element is not None and len(element) == 4 and all(e is not None for e in element) for row in matrix for element in row]))
##
    ##    # Extract the coordinates and values of known positions
    ##    known_coords = np.argwhere(np.array(matrix) != None)
    ##    known_values = np.array([matrix[i][j] for i, j in known_indices])
##
    ##    # Find the indices of None values in the matrix
    ##    none_indices = np.argwhere(np.array(matrix) == None)
##
    ##    # Extract the coordinates of None values
    ##    none_coords = tuple(map(tuple, none_indices))
##
    ##    # Perform interpolation for each element of the array separately
    ##    interpolated_values = [griddata(known_coords, known_values[:, i], none_coords, method='linear') for i in range(4)]
##
    ##    # Fill the matrix with interpolated values for each element of the array
    ##    for i in range(4):
    ##        matrix[tuple(np.transpose(none_indices))] = [interpolated_values[j][i] for j in range(len(none_coords))]
##
    ##    return matrix

    #square_grid=fill_matrix_with_interpolation(square_grid)
    return square_grid    

if __name__ == "__main__":
    import timeit
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    ##from termcolor import cprint
    from radtools.io.internal import load_template
    from radtools.io.tb2j import load_tb2j_model
    from radtools.magnons.dispersion import MagnonDispersion
    from radtools.decorate.stats import logo
    from radtools.spinham.constants import TXT_FLAGS
    from radtools.decorate.array import print_2d_array
    from radtools.decorate.axes import plot_hlines

    brillouin_primitive_vectors=np.zeros((3,3))
    brillouin_primitive_vectors[0]=[1,0,0]
    brillouin_primitive_vectors[1]=[-1/2,3,0]
    brillouin_primitive_vectors[2]=[0,0,1]
    plane_2d=[1,1,0]
    #threshold for understanding if there is a degeneracy
    threshold_k_grid=0.0000001
    shift_in_plane=[0,0]
    shift_in_space=[0,0,0]
    symmetry=[[0,0,0]]
    grid_spacing=0.1
    refinment_iteration=0
    refinment_spacing=0.01
    epsilon=refinment_spacing
    
    normalized_chosen_reciprocal_plane,k0,k1,refined_grid_2d=k_points_grid_2d_refinment_and_symmetry(brillouin_primitive_vectors,plane_2d,
        grid_spacing,
        shift_in_plane,
        shift_in_space,
        symmetry,
        refinment_spacing,
        refinment_iteration,
        threshold_k_grid)
    list_of_k_points=[]
    list_tmp=[]
    for i in range(0,k0):
        for j in range(0,k1):
            list_tmp=refined_grid_2d[i][j]
            for r in range(0,len(list_tmp)):
                list_of_k_points.append(list_tmp[r])

    list_of_k_points=mapping_to_square_grid(list_of_k_points,normalized_chosen_reciprocal_plane,epsilon)
    print(list_of_k_points)
    print(k0,k1)

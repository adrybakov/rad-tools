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

import numpy as np

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

def symmetry_transformation(
    k_origin, 
    k_point, 
    axis
    ):
    r"""
        symmetry transformation of a k vector, from k_origin to k_point, axis is the rotation vector: 
        i.e. its modulus is equal to the rotation angle, while its verso is the axis of rotation
      
        Parameters
        ----------
        k_origin: (,3) |double|
        k_point: (,3) |double|
        axis: (3) |array|
        Return
        ----------
        k_point rotated
       """
    rotation = np.sqrt(np.dot(axis, axis))
    # print(k_point-k_origin,axis)
    axis = axis / rotation
    k_vector = np.zeros(3)
    for i in range(0, 3):
        k_vector[i] = k_point[i] - k_origin[i]
    module = np.sqrt(np.dot(k_vector, k_vector))
    unit_k_vector = k_vector / module
    W = np.zeros((3, 3))
    I = np.eye(3)
    R = np.zeros((3, 3))
    W[0][1] = -axis[2]
    W[0][2] = axis[1]
    W[1][2] = -axis[0]
    W[1][0] = axis[2]
    W[2][0] = -axis[1]
    W[2][1] = axis[0]
    # numpy does it for your: matrix * scalar
    # Try to get rid of this function and use numpy
    def mul_mat_by_scalar(matrix, scalar):
        return [[scalar * j for j in i] for i in matrix]
    # numpy does it for your: matrix * matrix
    # Try to get rid of this function and use numpy
    def mul_mat_by_mat(matrixa, matrixb):
        product = np.zeros((3, 3))
        for k in range(0, 3):
            for l in range(0, 3):
                for r in range(0, 3):
                    product[k][l] = (
                        product[k][l] + matrixa[k][r] * matrixb[r][l]
                    )
        return product
    R = (
        I
        + mul_mat_by_scalar(W, np.sin(rotation))
        + mul_mat_by_scalar(
            mul_mat_by_mat(W, W), 2 * pow(np.sin(rotation / 2), 2)
        )
    )
    k_point_transformed = np.zeros(3)
    # Try to rewrite with numpy
    for i in range(0, 3):
        for j in range(0, 3):
            k_point_transformed[i] = (
                k_point_transformed[i] + R[i][j] * unit_k_vector[j]
            )
    k_point_transformed = k_point_transformed * module + k_origin
    return k_point_transformed

def symmetry_analysis(
    k_origin,
    k_origin_weight,
    k_points_subgrid,
    symmetry,
    threshold_k_grid,
):
    r"""
    the symmetry analisys is applied on a subset of k points (k_points_subgrid), the origin of the subsystem is considered as well (k_origin)
    in our case the origin is the point to which the refinment procedure is applied
    from the symmetry annalysis is clear the distribution of the origin weight (k_origin_weight) between the different k points of the subset

    Parameters
    ----------
    k_origin: (,3) |double|
    k_origin_weight: |double|
    k_points_subgrid: (:,3) |array|
    symmetry: |list of lists| a symmetry is a list of 3 elements (the versor is the axis of rotation, while the modulus is the angle)
    threshold_k_grid: |double| therhold to recognize a symmetry

    Return
    ----------
    k_points_subgrid_weight_tmp: (:,1) |array| distribution of the origin k point weight between the different k points
    """
    k_points_subgrid_weight_tmp = np.zeros(len(k_points_subgrid[:, 0]))
    # Use numpy
    for i in range(0, len(k_points_subgrid[:, 0])):
        k_points_subgrid_weight_tmp[i] = k_origin_weight / 4
    # if there are no symmetry operations, the weight on each subgrid k point is exactly 1/4 of the original weight
    if symmetry[0][0] == symmetry[0][1] == symmetry[0][2] == 0:
        return k_points_subgrid_weight_tmp
    
    # defining a matrix to save the degeneracies, after each symmetry operation
    check_degeneracy = np.zeros(
        (len(k_points_subgrid[:, 0]), len(k_points_subgrid[:, 0])),
        dtype=bool,
    )
    #saving number of degeneracies detected
    degeneracy = 0
    #each k point is compared to each other k point after the symmetry operation has been applied
    for i in range(0, len(k_points_subgrid[:, 0]) - 1):
        for j in range(i + 1, len(k_points_subgrid[:, 0])):
            r = 0
            #for each pair considering all the symmetry operations 
            while r < len(symmetry):
                k_point_transformed = symmetry_transformation(
                    k_origin, k_points_subgrid[j, :], symmetry[r]
                )
                flag = False
                #degeneracy_tmp is counting how many coordinates are equal
                degeneracy_tmp = 0
                #checking if all three coordinates of the two points are "equal"
                for s in range(len(k_point_transformed)):
                    flag = np.isclose(
                        k_points_subgrid[i, s],
                        k_point_transformed[s],
                        atol=threshold_k_grid,
                    )
                    if flag == True:
                        degeneracy_tmp = degeneracy_tmp + 1
                #if all three coordinates are equal, then the two points are saved as degenerate
                if degeneracy_tmp == len(k_point_transformed):
                    check_degeneracy[i][j] = True
                    degeneracy = degeneracy + 1
                    r = len(symmetry)
                else:
                    r = r + 1
    #if no degeneracy is detected, the common result is given
    if degeneracy == 0:
        for i in range(0, len(k_points_subgrid[:, 0])):
            k_points_subgrid_weight_tmp[i] = k_origin_weight / 4
        return k_points_subgrid_weight_tmp
    else:
        #if degeneracy is detected, of the degenerate points only one is chosen
        list = {}
        for i in range(0, len(k_points_subgrid[:, 0])):
            list[str(i)] = [i]
        #assuming all k points as not-degenerate, and creating a dictionary where to each eigenspace is associated an eigenvector(k_point)
        for l in range(0, len(k_points_subgrid[:, 0]) - 1):
            for j in range(l + 1, len(k_points_subgrid[:, 0])):
                if len(list) != 1:
                    #reading the degeneracy matrix, the eigenspaces are properly populated, in particular degenerate k points are put in the same eigenspace
                    if check_degeneracy[l][j] == True:
                        for key, values in list.items():
                            for value in values:
                                if value == l:
                                    positionl = key
                                if value == j:
                                    positionj = key
                        list_new = {}
                        #new list is created with the proper eigenspaces, but at the same time is compared with the old one, 
                        #to check if there are old degeneracies to take into account
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
                        #the new list is save to be used as an old one
                        list = {}
                        for key, values in list_new.items():
                            list[str(key)] = list_new[str(key)]
        #the chosen k points 
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
                #chosing of the degenerate points only the first one
                for value in values:
                    if count == 0:
                        k_points_subgrid_weight_tmp[value] = k_weight_tmp
                        count = count + 1
                    else:
                        k_points_subgrid_weight_tmp[value] = 0
        #returning the matrix with the respective weights
        return k_points_subgrid_weight_tmp

#function checking if the k points of the list k_points_subgrid are inside the brillouin zone
def check_inside_brillouin_zone(
    k_points_subgrid,
    reciprocal_vectors_2d,
    brillouin_primitive_vectors,
    plane_2d,
):
    matrix_crystal_to_cartesian=np.zeros((3,3))
    matrix_cartesian_to_crystal=np.zeros((3,3))
    k_points_subgrid_tmp=np.zeros((len(k_points_subgrid[:,0]),3))
    matrix_crystal_to_cartesian[:,0]=brillouin_primitive_vectors[0,:]
    matrix_crystal_to_cartesian[:,1]=brillouin_primitive_vectors[1,:]
    matrix_crystal_to_cartesian[:,2]=brillouin_primitive_vectors[2,:]
    matrix_cartesian_to_crystal=np.linalg.inv(np.matrix(matrix_crystal_to_cartesian))
    #writing k points in crystal coordinates
    for r in range(len(k_points_subgrid[:,0])):
        for s in range(0,3):
            for t in range(0,3):
                k_points_subgrid_tmp[r,s]= k_points_subgrid_tmp[r,s]+matrix_cartesian_to_crystal[s,t]*k_points_subgrid[r,t]
        count=0
        #if the crystal coordinate of a k point is beyond +-1, it means that the k point is outside the brillouin zone
        for l in range(0,3):
            #considering only those primitive vectors defining the 2D plane chosen 
            if plane_2d[l]!=0:
                if abs(k_points_subgrid_tmp[r,l])>=1:
                    for s in range(0,3):
                        #translating the respective cartesian coordinates
                        k_points_subgrid[r,s]=k_points_subgrid[r,s]-int(k_points_subgrid_tmp[r,l])*reciprocal_vectors_2d[count][s]
            #progressing on the 2D plane primitive vectors
            count=count+1
    #returning the properly translated k points
    return k_points_subgrid[:,:3]

def local_refinment(
    refined_grid,
    reciprocal_vectors_2d,
    not_refined_grid,
    refinment_spacing,
    refinment_iteration,
    symmetry,
    threshold_k_grid,
    brillouin_primitive_vectors,
    plane_2d
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
                k_points_subgrid[:,:3]=check_inside_brillouin_zone(k_points_subgrid[:,:3],reciprocal_vectors_2d,brillouin_primitive_vectors,plane_2d)
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
                            plane_2d
                        )


def dynamical_refinment(
    k_point_to_refine,
    reciprocal_vectors_2d,
    refinment_spacing,
    refinment_iteration,
    symmetry,
    threshold_k_grid,
    brillouin_primitive_vectors,
    plane_2d
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
        plane_2d
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
    #np.dot(v1,v2) same as v1 @ v2
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
    with open("k_points.txt", "w") as file:
        for i in range(0, k0):
            for j in range(0, k1):
                print(
                    "k grid: {}% and {}%".format(int(i / k0 * 100), int(j / k1 * 100)),
                    end="\r",
                    flush=True,
                )
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
                    plane_2d
                )
                refined_grid_tmp = np.reshape(
                    refined_grid_tmp, (int(len(refined_grid_tmp) / 4), 4)
                )
                k_points_grid_2d[i][j] = np.zeros((len(refined_grid_tmp), 4))
                k_points_grid_2d[i][j] = refined_grid_tmp
                file.write(f"{i} {j} {k_points_grid_2d[i][j]}  \n")
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

##if __name__ == "__main__":
##    import timeit
##    import os
##    import matplotlib.pyplot as plt
##    import numpy as np
##    from termcolor import cprint
##    from radtools.io.internal import load_template
##    from radtools.io.tb2j import load_tb2j_model
##    from radtools.magnons.dispersion import MagnonDispersion
##    from radtools.decorate.stats import logo
##    from radtools.spinham.constants import TXT_FLAGS
##    from radtools.decorate.array import print_2d_array
##    from radtools.decorate.axes import plot_hlines
##
##    brillouin_primitive_vectors=np.zeros((3,3))
##    brillouin_primitive_vectors[0]=[1,0,0]
##    brillouin_primitive_vectors[1]=[-1/2,3,0]
##    brillouin_primitive_vectors[2]=[0,0,1]
##    plane_2d=[1,1,0]
##    #threshold for understanding if there is a degeneracy
##    threshold_k_grid=0.0000001
##    shift_in_plane=[0,0]
##    shift_in_space=[0,0,0]
##    symmetry=[[0,0,np.pi]]
##    grid_spacing=0.01
##    refinment_iteration=3
##    refinment_spacing=0.005
##    normalized_chosen_reciprocal_plane,k0,k1,refined_grid_2d=k_points_grid_2d_refinment_and_symmetry(brillouin_primitive_vectors,plane_2d,
##        grid_spacing,
##        shift_in_plane,
##        shift_in_space,
##        symmetry,
##        refinment_spacing,
##        refinment_iteration,
##        threshold_k_grid)
##    print(k0,k1)
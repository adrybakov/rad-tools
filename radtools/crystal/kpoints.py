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
from scipy.spatial import Delaunay
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

# function checking if the k points in the list k_list are inside the Brillouin zone
def check_inside_brillouin_zone(
    k_point_list,
    brillouin_primitive_vectors_3d
):
    r"""
        check if any point in a list is inside the first brillouin zone, in case shifting it back to the first brillouin zone
        Parameters
        ----------
        k_point_list: (N, 3) :|array_like| kx,ky,kz (k points are given in cartesian coordinates)
        brillouin_primitive_vectors_3d: (3,3) :|array_like| columns: kx,ky,kz, rows: b1,b2,b3
        Returns
        -------
        transformed_k_point_list: (N, 3) :|array_like| kx,ky,kz (k points are given in cartesian coordinates)
        """

    matrix_transformation_crystal_to_cartesian = np.zeros((3,3))
    matrix_transformation_cartesian_to_crystal = np.zeros((3,3))
    
    number_k_point_elements = k_point_list.shape[0]
    transformed_k_point_list = np.zeros((number_k_point_elements,3),dtype=float)

    matrix_transformation_crystal_to_cartesian = np.matrix(brillouin_primitive_vectors_3d).transpose()

    matrix_transformation_cartesian_to_crystal=np.linalg.inv(matrix_transformation_crystal_to_cartesian)

    # writing k points in crystal coordinates
    for i in range(number_k_point_elements):
        transformed_k_point_list[i,:] = matrix_transformation_cartesian_to_crystal @ k_point_list[i,:]

    # checking if any point of the list is outside the first brillouin zone
    for r in range(3):
        modules=np.absolute(transformed_k_point_list[:,r])
        for i in range(number_k_point_elements):
            if modules[i] >=1:
                for d in range(3):
                    k_point_list[i,d]-=int(transformed_k_point_list[i,d])*brillouin_primitive_vectors_3d[r,d]
    return k_point_list

def check_inside_closed_shape(
        subset_triangles,
        k_point_list,
        all_k_point_list
):
    r"""
    Given a set of triangles with one common vertex and adjacent sides (pair by pair), 
    the function checks if the k points in the list are outside any of the triangles
    if any of the k points is outside it is erased from the k point list
    Parameters
    ----------
    triangles: (N,3) |array_like| the triangles with a common origin and adjacent sides 
                each triangle has three indices, which point to one of the k points
                from the all k point list

    all_k_point_list (N,3) |array_like| (kx,ky,kz)
    k_point_list: (k,3) |array_like| the k points that need to be checked (kx,ky,kz)
    Returns
    ----------
    k_point_list: (k-s,3) |array_like| the k points that resulted to be inside one of the subset of triangles
    """

    number_elements=k_point_list.shape[0]
    number_triangles=subset_triangles.shape[0]

    # for each triangle the vecotrial products between the different vertices, and between the vertices and the to-be-checked k points
    # are used to calculate the barycentric weights, which give information about the to-be-checked k point being inside or outside the triangle itself
    vectorial_products=np.zeros(6)
    weights=np.zeros(3)
    elements_to_erase=[]
    for j in range(number_triangles):
        indices=np.unique(np.where(all_k_point_list[:][0]==subset_triangles[j])[0],axis=0)
        vectorial_products[0]=np.linalg.norm(np.cross(all_k_point_list[indices[0]],all_k_point_list[indices[1]]))
        vectorial_products[1]=np.linalg.norm(np.cross(all_k_point_list[indices[1]],all_k_point_list[indices[2]]))
        vectorial_products[2]=np.linalg.norm(np.cross(all_k_point_list[indices[2]],all_k_point_list[indices[0]]))
        normalization=vectorial_products[0]+vectorial_products[1]+vectorial_products[2]
        for i in range(number_elements):
            vectorial_products[3]=np.linalg.norm(np.cross(k_point_list[i],all_k_point_list[indices[1]]-all_k_point_list[indices[2]]))
            vectorial_products[4]=np.linalg.norm(np.cross(k_point_list[i],all_k_point_list[indices[2]]-all_k_point_list[indices[0]]))
            vectorial_products[5]=np.linalg.norm(np.cross(k_point_list[i],all_k_point_list[indices[0]]-all_k_point_list[indices[1]]))
            weights[1]=(vectorial_products[1]+vectorial_products[3])/normalization
            weights[2]=(vectorial_products[2]+vectorial_products[4])/normalization
            weights[3]=(vectorial_products[0]+vectorial_products[5])/normalization
        if np.min(weights) < 0 or np.max(weights)> 1:
            elements_to_erase.append(i)
    
    # eliminating all k points not contained in any of the triangles
    np.delete(k_point_list,elements_to_erase,axis=0)

    return k_point_list


# function applying symmetry analysis to a list of k points, the symmetry operations considered are the point group ones
# the fixed point is given; the analysis aim is to redistribute the fixed point weight between the list of k points
# the found symmetry orbits are counted in the redistribution of the fixed point weight a number of times equal to the number of k points inside the orbit itself
def symmetry_analysis(
    k_origin,
    k_origin_weight,
    k_list,
    symmetries,
    threshold,
):
    r"""
    Given a list of k points "k_list" and an origin "k_origin" with a certain weight "k_origin_weight", different symmetry operations are applied to the list of k points, keeping the origin fixed.
    The origin weight is propelry distributed between the points:
        the orbits of the symmetry operations acquire a certain percentage of the origin weight, which takes into account how many orbits are found 
        of the orbits one k point (from the k list) is chosen as representative, and the weight of the orbit is assigned only to it
    Parameters
    ----------
    k_origin:(3) |array_like| fixed point coordinates (kx,ky,kz), fixed point considered in the point-group symmetry operations
    k_origin_weight : |float|  weight to distribute between the k points of the list 
    k_list : (N, 3) |array_like| list of k points considered in the symmetry analysis (kx,ky,kz)
    symmetries : list of lists, a symmetry is a list of 3 elements (the versor is the axis of rotation, while the modulus is the angle).
    threshold : |float| threshold to recognize a symmetry.
    Returns
    -------
    new_k_list: (N,4) :|array_like| list of k points coordinates with respcetive weights from the symmetry analysis
    """

    number_elements = k_list.shape[0]  

    new_k_list = np.c_[ k_list, (k_origin_weight / number_elements) * np.ones(number_elements, dtype=float)]

    # if there are no symmetry operations, than the weight on each k point is exactly 1/N of the original weight
    if symmetries is None or (symmetries[0][0] == symmetries[0][1] == symmetries[0][2] == 0):
        return new_k_list

    # defining a matrix to save the degeneracies, after each symmetry operation
    check_degeneracy = np.zeros((number_elements, number_elements), dtype=bool)
    # saving flag if no degeneracy is detected
    no_degeneracy = True
    
    k_eigenspaces=np.zeros(number_elements)
    for i in range(number_elements):
        k_eigenspaces[i]=i

    # each k point is compared to each other k point after each of the symmetry operations has been applied
    for i in range(number_elements - 1):
        for j in range(i + 1, number_elements):
            # for each pair considering all the symmetry operations
            for rotvec in symmetries:
                k_point_transformed = (
                    Rotation.from_rotvec(rotvec).as_matrix()
                    @ (k_list[j] - k_origin)
                    + k_origin
                )
                # if any degeneracy is detected, the check on all the symmetry operations given is stopped
                if np.allclose(
                    k_list[i], k_point_transformed, atol=threshold
                ):
                    check_degeneracy[i][j] = True
                    no_degeneracy = False
                    break

    # if no degeneracy is detected, the no-symmetry-operations result is given
    if no_degeneracy:
        return new_k_list
    else:
        # if degeneracy is detected, of the degenerate points only one is chosen as representative of the respective orbit (the first to appear in the k list)
        # assuming all k points as not-degenerate, and creating a dictionary
        # where to each eigenspace (orbit) is associated a not null k_point
        for i in range(0,number_elements-1):
            for j in range(i+1,number_elements):
                if check_degeneracy[i][j]==True:
                    min_value=min(i,j)
                    k_eigenspaces[k_eigenspaces==min_value]=min_value
        
        # counting number of different eigenspaces
        number_eigenspaces=np.unique(k_eigenspaces)
        weights=k_origin_weight/number_eigenspaces
        # ordering eigenspaces values
        ordered_k_eigenspaces=np.zeros(number_elements)
        counting=0
        for i in range(number_elements):
            if counting<number_eigenspaces:
                ordered_k_eigenspaces[k_eigenspaces==k_eigenspaces[i]]=counting
                counting+=1
        # considering a representative for each eigenspace
        for i in range(number_eigenspaces):
            list_positions=[ordered_k_eigenspaces==i]
            new_k_list[list_positions[0]][3]=weights
            new_k_list[list_positions[0:]][3]=0

        # returning the matrix with the respective weights
        return new_k_list

# local refinment of the k list in the plane (parallelepiped) pointed out by the brillouin primitive vectors given (it can be 2D, 3D or ..)
def local_refinment(
    new_k_list,
    old_k_list,
    refinment_spacing,
    refinment_iteration,
    symmetries,
    threshold,
    brillouin_primitive_vectors_3d,
    brillouin_primitive_vectors_2d,
    normalized_brillouin_primitive_vectors_2d

):
    r"""
    Starting from a set of k points (old_k_list) with certain weights, to each k point iteratively (refinment_iteration) a new subset of k points is associated,
    a symmetry analysis (symmetry) is performed on the subset, and consequently the initial weight of the k point is distributed in the subset of k points.
    This procedure goes on for the given number of refinment iterations.
    At each iteration for each k point four new k points are generated along the reciprocal primitive vectors of the selcted plane (2d refinment), at a distance from the original k point 
    equal to the refinment spacing (the refinment spacing is halved at each iteration).
    Obviously no refinment is applied to k points with null weight.

    if check_inside_dominion == True instead of checking that the point is inside the 2d plane of the brillouin zone, it is checked that it is inside a circle of a certain radius
    (radius), and if it is outside the circle it is downfolded back 

    Parameters
    ----------
    new_k_list: (,4) |array_like| list of values produced by the refinment procedure (expected to be empty at the beginning of the refinment procedure)
    old_k_list: (,4) |list| list of values give to the refinment procedure (kx,ky,kz,weight)
    refinment_spacing: |double| initial spacing given to generate the refinment, half of the preceding refinment spacing is considered at each refinment iteration
    refinment_iteration: |int| number of refinment iterations considered
    symmetry: |list of lists| a symmetry is a list of 3 elements (the versor is the axis of rotation, while the modulus is the angle)
    threshold: |double| threshold to recognize a symmetry
    brillouin_primitive_vectors : (2,3) |array_like| vectors definining the 2D plane in the reciprocal space, chosen for the refinment procedure
    brillouin_primitive_vectors_all: (3,3) |array_like| vectors definining the 3D reciprocal space
    (the k points of the list are expected to be in the same 2D plane)
    """
    length_old_k_list=old_k_list.shape[0]
    
    if refinment_iteration == 0:
        for i in range(length_old_k_list):
            if old_k_list[i,3]!=0:
                new_k_list=np.append(new_k_list,old_k_list[i,:])
    else:
        iter = 0
        while iter != length_old_k_list:
            # reading the k points inside the list and refine around them
            if length_old_k_list == 1:
                k_tmp = old_k_list[0,:3]
                k_tmp_weight = old_k_list[0,3]
                iter = length_old_k_list
            else:
                k_tmp = old_k_list[iter,:3]
                k_tmp_weight = old_k_list[iter,3]
                iter = iter + 1
            # refinment procedure is applied to the single k point if its weight is not null
            if k_tmp_weight != 0:
                # a subgrid of points is associated to each k point
                k_tmp_subgrid = np.zeros((4, 4))
                k_tmp_subgrid[:,:3] = k_tmp
                k_tmp_subgrid[0,:3] +=  normalized_brillouin_primitive_vectors_2d[0,:]*refinment_spacing
                k_tmp_subgrid[1,:3] -=  normalized_brillouin_primitive_vectors_2d[0,:]*refinment_spacing
                k_tmp_subgrid[2,:3] +=  normalized_brillouin_primitive_vectors_2d[1,:]*refinment_spacing
                k_tmp_subgrid[3,:3] -=  normalized_brillouin_primitive_vectors_2d[1,:]*refinment_spacing
                
                k_tmp_subgrid = symmetry_analysis(
                    k_tmp,
                    k_tmp_weight,
                    k_tmp_subgrid[:,:3],
                    symmetries,
                    threshold
                )
                
                k_tmp_subgrid=np.delete(k_tmp_subgrid, np.where(k_tmp_subgrid[:,3]==0),axis=0)

                k_tmp_subgrid[:,:3] = check_inside_brillouin_zone(
                    k_tmp_subgrid[:,:3],
                    brillouin_primitive_vectors_3d,
                )

                # applying to the points of the subgrid the refinment procedure
                new_refinment_spacing = refinment_spacing / 2
                
                local_refinment(
                    new_k_list,
                    k_tmp_subgrid,
                    new_refinment_spacing,
                    refinment_iteration - 1,
                    symmetries,
                    threshold,
                    brillouin_primitive_vectors_3d,
                    brillouin_primitive_vectors_2d,
                    normalized_brillouin_primitive_vectors_2d 
                )

def interpolation_k_points_weights(
    k_points_list_with_weights,
    brillouin_primitive_vectors_2d,
    number_elements_interpolation
    ):
    r"""
    Through an interpolation, the given list of k points defined in a 2D plane (generated by the given brillouin_primitive_vectors) is used to produce a list of ordered k points (properly weighted).
    The number of k points produced is equal to number_elements x number_elements.

    Parameters 
    ----------
    k_points_list_with_weights : (,4) |array_like| : list of k points (kx,ky,kz,w) where w is the respective weight of the k point
    brillouin_primitive_vectors: (2,3) |array_like| : 2D plane in the reciprocal space (the list of k points given has to be generated from these reciprocal vectors)
    number_elements_interpolation: (2) number of elements considered along the two reciprocal vectors given, in the interpolation procedure
    
    Return
    -------------
    new_k_points_list : (:6) |array_like| : list of k points (kx,ky,kz,w,i,j) where w is the respective weight of the k point, and (i,j) the indices ordering the 2D k point grid 
    new_k_points_grid : (n0,n1) |array_like| (4) (kx,ky,kz,w) the list of k points ( new_k_points_list_with_weights_and_ordering) is written as a grid
    n0,n1 ...
    note: the weights, being added some points in the interpolation, are renormalized over the entire set of points (it is assumed that the initial list is already normalized)
    """
    # considering first a mapping from the list of k vectors in the 3D cartesian coordinates with the respective weight, to a list of k vectors in the 2D crystal coordinates with the respective weight
    # the projection of the k vectors on the primitive reciprocal vectors (defining the 2D plane given) is therefore considered, instead of the cartesian coordinates
    
    number_k_points=k_points_list_with_weights.shape[0]
    k_vectors=k_points_list_with_weights[:,:3]
    k_vectors_projections=np.zeros((number_k_points,2))
    
    for i in range(number_k_points):
        k_vectors_projections[i,:]=brillouin_primitive_vectors_2d@k_vectors[i,:]
    
    #this pairs (the projection coordinates) are localized in a box [0,1)x[0,1), making an interpolation quite straightforward and reliable
    #this is the main motivation to consider the projections along the two directions
    xmin,xmax=0,1
    ymin,ymax=0,1
    
    #number of points considered in the interpolation (a differentiation between the two directions can be done)
    if number_elements_interpolation is None:
        ny,nx=number_k_points,number_k_points
        number_elements_interpolation=np.zeros(2,dtype=int)
        number_elements_interpolation[0]=number_k_points
        number_elements_interpolation[1]=number_k_points
    else:
        nx,ny=number_elements_interpolation[0],number_elements_interpolation[1]

    #generate a regular grid to interpolate the data
    xi=np.linspace(xmin,xmax,nx)
    yi=np.linspace(ymin,ymax,ny)
    xi,yi=np.meshgrid(xi,yi)

    #considering the projections as a pair of data with a corresponding weight z
    x=k_vectors_projections[:,0]
    y=k_vectors_projections[:,1]
    z=k_points_list_with_weights[:,3]

    #interpolation 
    zi=griddata((x,y),z,(xi,yi),method='nearest')

    #normalizing the new list to one
    total_norm=np.sum(zi)
    zi=list(map(lambda x: x/total_norm,zi))

    ##now mapping back to the brillouin zone saving the indices
    new_k_points_list=np.zeros((number_elements_interpolation[0]*number_elements_interpolation[1],6))
    new_k_points_grid=np.zeros((number_elements_interpolation[0],number_elements_interpolation[1],4))
    coordinates=np.zeros(3)
    for i in range(number_elements_interpolation[0]):
        for j in range(number_elements_interpolation[1]):
            vector=np.asanyarray([xi[i,j],yi[i,j]])
            coordinates=(brillouin_primitive_vectors_2d.T)@vector
            new_k_points_grid[i][j][:3]=coordinates
            new_k_points_grid[i][j][3]=zi[i][j]
            new_k_points_list[i*number_elements_interpolation[1]+j,:]=np.asarray(coordinates,zi[i][j],i,j)

    return(new_k_points_list,new_k_points_grid,number_elements_interpolation[0],number_elements_interpolation[1])

# Generation of a 2D k points grid through a refinment procedure
# Because of the disorder of the k points due to the refinment procedure, an interpolation scheme can been added to obtain a propelry ordered k points grid
# however this approach diminuishes the advantages to consider a symmetry refinment procedure 
# therefore a triangulation can be considered instead
def k_points_grid_generator_2D(
    brillouin_primitive_vectors_3d,
    chosen_plane,
    initial_grid_spacing,
    shift_in_plane,
    shift_in_space,
    symmetries,
    threshold,
    refinment,
    refinment_spacing,
    refinment_iteration,
    interpolation,
    triangulation
):
    r"""
    2D k points grid generator
    Parameters
    ----------
    brillouin_primitive_vectors_3d: (3,3) |array_like| (columns: kx,ky,kz, rows: b1,b2,b3)
    chosen_plane: (,3) |array_like|  ex. [0,1,1] the plane is the b2 x b3 plane in the reciprocal space where the bi are the brillouin primitive vectors
    initial_grid_spacing: |double| spacing of the k grid before the refinment procedure
    shift_in_plane : (,2) |array_like| shift of the chosen reciprocal 2D plane with respecet to the crystal coordinates
    shift_in_space : (,3) |array_like|  shift of the chosen reciprocal 2D plane with respecet to the cartesian coordinates 
    (this allow to build a 3D k points grid, considering different shifts along the third brillouin primitive vector, if there are no symmetries along this direction to consider)
    symmetries: |list of lists| a symmetry is a list of 3 elements (the versor is the axis of rotation, while the modulus is the angle)
    refinment_spacing: |double| initial spacing given to generate the refinment, half of the preceding refinment spacing is considered at each refinment iteration
    refinment_iteration: |int| number of refinment iterations considered
    threshold: |double| therhold to recognize a symmetry
    
    Return
    -------------
    k0,k1 (1,1) |int| the number of k points along the two directions generating the chosen 2D reciprocal plane
    k points grid 2D k0xk1 grid
    """

    if not refinment_spacing:
        refinment_spacing=initial_grid_spacing/2
    if not shift_in_plane:
        shift_in_plane=[0,0]
    if not shift_in_space:
        shift_in_space=[0,0,0]
    
    brillouin_primitive_vectors_2d = np.zeros((2, 3))
    normalized_brillouin_primitive_vectors_2d = np.zeros((2, 3))

    count = 0
    for i in range(3):
        if chosen_plane[i] != 0:
            brillouin_primitive_vectors_2d[count] = brillouin_primitive_vectors_3d[i]
            normalized_brillouin_primitive_vectors_2d[count] = brillouin_primitive_vectors_2d[count]/(brillouin_primitive_vectors_2d[count]@brillouin_primitive_vectors_2d[count])
            count += 1

    # initial 2D grid dimension
    k0 = int(
        (brillouin_primitive_vectors_2d[0] @ brillouin_primitive_vectors_2d[0])/ initial_grid_spacing
    )
    k1 = int(
        (brillouin_primitive_vectors_2d[1] @ brillouin_primitive_vectors_2d[1]) / initial_grid_spacing
    )
    initial_weights = 1 / (k0 * k1)
    
    if refinment == True:
        # building the initial 2D k points grid
        k_points_grid = np.zeros((k0*k1,4),dtype=float)
        count=0
        for i in range(k0):
            for j in range(k1):
                k_points_grid[count,:3] = (float(i + shift_in_plane[0]) / k0) * (shift_in_space + brillouin_primitive_vectors_2d[0]) \
                    + (float(j + shift_in_plane[1]) / k1) * (shift_in_space + brillouin_primitive_vectors_2d[1])
                k_points_grid[count,3] = initial_weights
                count=count+1

        # applying the refinment procedure to the 2d k points grid
        refined_k_points_list = []
        local_refinment(
            refined_k_points_list,
            k_points_grid,
            refinment_spacing,
            refinment_iteration,
            symmetries,
            threshold,
            brillouin_primitive_vectors_3d,
            brillouin_primitive_vectors_2d,
            normalized_brillouin_primitive_vectors_2d  
        )     
        #transforming from a list to array_like elements
        refined_k_points_list = np.reshape(refined_k_points_list,(int(len(refined_k_points_list)/4), 4))

        if interpolation == True:
            #interpolating and ordering the 2D k points grid
            new_k_points_list,new_k_points_grid,n0,n1=interpolation_k_points_weights(refined_k_points_list,brillouin_primitive_vectors_2d,None)
            return (new_k_points_list,new_k_points_grid,n0,n1)
        elif triangulation == True:
            refined_k_points_list,triangles=k_points_list_tesselation_2d(refined_k_points_list,brillouin_primitive_vectors_2d)
            return refined_k_points_list, triangles
        else:
            return(refined_k_points_list)
    else:
        k_points_grid_not_refined = np.zeros((k0,k1,4))
        k_points_list_not_refined = np.zeros((k0*k1,6))
        for i in range(k0):
            for j in range(k1):
                k_points_grid_not_refined[i,j][:3] = (float(i + shift_in_plane[0]) / k0) * (shift_in_space + brillouin_primitive_vectors_2d[0]) \
                    + (float(j + shift_in_plane[1]) / k1) * (shift_in_space + brillouin_primitive_vectors_2d[1])
                k_points_grid_not_refined[i,j][3] = initial_weights
                k_points_list_not_refined[i*k1+j,:]=[k_points_grid_not_refined[i,j][:],i,j]
                
        return (k_points_list_not_refined,k_points_grid_not_refined,k0,k1)

def k_points_list_tesselation_2d(
    k_points_list,
    brillouin_primitive_vectors_2d,
    count
):
    r"""
    the sparse k points (i,kx,ky,kz,w) are linked togheter through a triangulation procedure (Delaunay procedure)
    the so-built triangles are listed, each triangle is characterized by three indices pointing to the vertices position in the k poins list  
    the ordering of the vertices is so to assure a cloack-wise path along the triangle
    Parameters
    ---------
    k_points_list: (n,4) |array_like| (kx,ky,kz,w)
    brillouin_primitive_vectors_2d: (2,3) |array_like|
    
    Returns
    ---------
    k_points_list: (n,5) |array_like| (i,kx,ky,kz,w)
    triangles: (s,3) |array_like| 
    """
    number_elements=k_points_list.shape[0]
    k_points_list_projections=np.zeros((number_elements,3))
    new_k_points_list=np.zeros((number_elements,5))
    # choosing one point as an origin in the 2d brillouin plane
    origin=k_points_list[0][:3]
    # calculating the projections of the different k points on the primitive vectors
    for i in range(number_elements):
        k_points_list_projections[i,:2]=np.dot(k_points_list[i][:3]-origin,brillouin_primitive_vectors_2d)
        k_points_list_projections[i,2]=i
        new_k_points_list[i,1:]=k_points_list[i,:]
        new_k_points_list[i,0]=i+count

    # triangulation of the set of points
    triangles = Delaunay(k_points_list_projections)
    # these are the different triangles, where each element has a pair of number representating the ordering in the triangle itself and a number representing the respective element of the list
    number_of_triangles = len(triangles) 
    triangles_renamed = np.zeros((number_of_triangles,3),dtype=float)
    for i in range(number_of_triangles):
        # looking for the clock-wise path
        sorted(triangles[i],key=lambda x: x[1])
        if triangles[i][0,:1]==np.min(triangles[i][:,0]):
            sorted(triangles[i],key=lambda x: x[0])
        for s in range(3):
            triangles_renamed[i][s]=triangles[i][2]
    # for each triangle the ordering of the vertices is now clock-wise
    return new_k_points_list,triangles_renamed

def dynamical_refinment_tesselation_2d(
    triangles,
    positions_k_points_to_refine,
    all_k_points_list,
    refinment_iteration,
    symmetries,
    threshold,
    brillouin_primitive_vectors_3d,
    brillouin_primitive_vectors_2d,
    normalized_brillouin_primitive_vectors_2d
):
    new_k_points_list=[]
    new_triangles=[]
    new_k_points_list.append(all_k_points_list)
    number_elements_to_refine=len(positions_k_points_to_refine)
    indices_k_points_to_refine=np.unique(np.where(all_k_points_list[:][0]==positions_k_points_to_refine)[0],axis=0)
    count=len(all_k_points_list)
    for i in range(number_elements_to_refine):
        triangles_to_refine=[]
        positions_triangles_to_refine=np.unique(np.where(triangles_to_refine==positions_k_points_to_refine[i])[0])
        number_subset_triangles=len(positions_triangles_to_refine)
        for r in range(number_subset_triangles):
            triangles_to_refine.append(triangles[positions_triangles_to_refine[r]])
        
        border_k_points=[]
        # here we are counting vertices twice
        distances=np.zeros(number_subset_triangles*3)
        for j in range(number_subset_triangles):
            for s in range(3):
                distances[j*3+s]=np.linalg.norm(all_k_points_list[triangles[positions_triangles_to_refine[j]][s]]-all_k_points_list[indices_k_points_to_refine[i]])
                if distances[j*3+s]==0:
                    distances[j*3+s]=np.nan
                else:
                    border_k_points.append(triangles[positions_triangles_to_refine[j]][s])
        # remove any double counting
        border_k_points=np.unique(np.array(border_k_points)).tolist()
        # applying local refinment to the selected point
        # chose as distance the minimum distance with respect to the border
        # here some problems can emerge 
        minimum_distance=np.min(distances)/2
        refined_k_points_list=[]
        local_refinment(
            refined_k_points_list,
            all_k_points_list[indices_k_points_to_refine[i]],
            minimum_distance,
            refinment_iteration,
            symmetries,
            threshold,
            brillouin_primitive_vectors_3d,
            brillouin_primitive_vectors_2d,
            normalized_brillouin_primitive_vectors_2d  
        )
        # transforming from a list to array_like elements
        refined_k_points_list = np.reshape(refined_k_points_list,(int(len(refined_k_points_list)/4), 4))
        # check if the points are inside the border 
        refined_k_points_list=check_inside_closed_shape(
            triangles_to_refine,
            refined_k_points_list,
            all_k_points_list)
        # triangulation of the new k points
        border_k_points.append(refined_k_points_list)
        border_k_points = np.reshape(refined_k_points_list,(int(len(refined_k_points_list)/4), 4))    
        k_point_list,added_triangles=k_points_list_tesselation_2d(
            border_k_points,
            brillouin_primitive_vectors_2d,
            count
            )
        count=count+len(refined_k_points_list)
        new_k_points_list.append(k_point_list)
        new_triangles.append(added_triangles)

    new_triangles=np.unique(np.array(new_triangles)).tolist()
    return new_k_points_list,new_triangles


# Function considering one point in a list of k points and applying a refinment procedure to the selected k point (in the plane pointed out by brillouin_primitive_vectors)
def dynamical_refinment(
    old_k_list,
    position_elements_to_refine,
    brillouin_primitive_vectors,
    refinment_iteration,
    symmetries,
    threshold,
    ordering
):
    r"""
    Dynamical refinment around a k point in a 2D plane
    Note: the k points are supposed to be characterized by (kx,ky,kz,w,i,j) their cartesian coordinates, their weight and their pair of indices ordering them in the 2D plane
    
    Parameters
    -------
    old_k_list: (:6) |array_like| : list of k points (kx,ky,kz,w,i,j) where w is the respective weight of the k point, and (i,j) the indices ordering the 2D k point grid 
    elements_to_refine: (n,4) |array_like| : k points to refine (kx,ky,kz,w)
    position_elements_to_refine: (n) |array_like| (int) the position in the list of the k point to refine
    brillouin_primitive_vectors: (2,3) |array_like| : 2D plane in the reciprocal space (the list of k points given has to be generated from these reciprocal vectors)
    symmetries: |list of lists| a symmetry is a list of 3 elements (the versor is the axis of rotation, while the modulus is the angle)
    refinment_iteration: |int| number of refinment iterations considered
    threshold: |double| therhold to recognize a symmetry
    number_elements_interpolation: (2) number of elements considered along the two reciprocal vectors given, in the interpolation procedure

    Return
    -------------
    k0,k1 (1,1) |int| the number of k points along the two directions generating the chosen 2D reciprocal plane
    k points grid 2D k0xk1 grid
    """
    # applying the refinment procedure to the k point selected
    refined_grid = []
    # if more than one point is given
    for i in range(len(position_elements_to_refine)):
        
        # position of the k point selected in the list
        element_to_refine=old_k_list[position_elements_to_refine[i]]
        del old_k_list[position_elements_to_refine[i]]
        
        minimal_distance = np.min(np.abs(old_k_list[:,:3]-element_to_refine[i][:3]))
        refinment_spacing = minimal_distance/2
        
        local_refinment(
            refined_grid,
            element_to_refine,
            refinment_spacing,
            refinment_iteration,
            symmetries,
            threshold,
            brillouin_primitive_vectors
        )     

    # transforming from a list to array_like elements
    refined_grid = np.reshape(refined_grid, len(refined_grid)/4, 4)

    # eliminating indices from the old k list
    del old_k_list[:,4:]

    # adding the new points to the old list
    old_k_list=np.vstack([old_k_list,refined_grid])
    
    if ordering == False:
        return old_k_list
    else:
        # interpolating and ordering the 2D k points grid
        new_k_list, new_k_grid, n0, n1=interpolation_k_points_weights(old_k_list,brillouin_primitive_vectors,None)
        return (new_k_list,new_k_grid,n0,n1)

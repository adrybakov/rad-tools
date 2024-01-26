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

# function checking if the k points in the list k_list are inside the Brillouin zone
def check_inside_brillouin_zone(
    k_list,
    brillouin_primitive_vectors_3d
):
    r"""
        check if any point in a list is inside the first brillouin zone, in case shifting it back to the first brillouin zone
        Parameters
        ----------
        k_list: (N, 3) :list kx,ky,kz (k points are given in cartesian coordinates)
        brillouin_primitive_vectors_3d: (3,3) :list columns: kx,ky,kz, rows: b1,b2,b3 
        Returns
        -------
        new_k_list: (N, 3) :list kx,ky,kz (k points are given in cartesian coordinates)
        """
    #k_list=np.array(k_list)

    matrix_transformation_crystal_to_cartesian = np.zeros((3, 3))
    matrix_transformation_cartesian_to_crystal = np.zeros((3, 3))
    number_elements = len(k_list.shape[0])
    transformed_k_list = np.zeros(number_elements, 3)

    matrix_transformation_crystal_to_cartesian[:,0] = brillouin_primitive_vectors_3d[0, :]
    matrix_transformation_crystal_to_cartesian[:,1] = brillouin_primitive_vectors_3d[1, :]
    matrix_transformation_crystal_to_cartesian[:,2] = brillouin_primitive_vectors_3d[2, :]
    matrix_transformation_cartesian_to_crystal = np.linalg.inv(np.matrix(matrix_transformation_crystal_to_cartesian))

    # writing k points in crystal coordinates
    transformed_k_list = matrix_transformation_cartesian_to_crystal @ k_list
    # checking if any point of the list is outside the first brillouin zone
    for r in range(3):
        modules=np.absolute(transformed_k_list[:,r])
        for i in range(number_elements):
            if modules[i] >=1:
                k_list[i,:]=k_list[i,:]-int(transformed_k_list[i,:])*brillouin_primitive_vectors_3d[r,:]
    return k_list

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
    Given a list of k points "k_list" and an origin "k_origin" with a certain weight "k_origin_weight", different symmetry operations are
    applied to the list of k points, keeping the origin fixed, and the origin weight propelry distributed between the points
    The orbits of the symmetry operations acquire a certain percentage of the origin weight, which takes into account how many orbits are found 
    Of the orbits one k point is chosen as representative, and the weight of the orbit assigned only to it
    Parameters
    ----------
    k_origin:(3) |array| fixed point coordinates (kx,ky,kz), fixed point considered in the point-group symmetry operations
    k_origin_weight : float  weight to distribute between the k points of the list 
    k_list : (N, 3) |array| list of k points considered in the symmetry analysis (kx,ky,kz)
    symmetries : list of lists
        a symmetry is a list of 3 elements
        (the versor is the axis of rotation, while the modulus is the angle).
    threshold : float
        threshold to recognize a symmetry.
    Returns
    -------
    new_k_list: (N,4) :|array| list of k points with respcetive weight from the symmetry analysis
    """

    #k_list = np.array(k_list)
    #k_origin = np.array(k_origin)

    number_elements = k_list.shape[0]  

    new_k_list = np.c_[ k_list, (k_origin_weight / number_elements) * np.ones(number_elements, dtype=float)]
    new_k_list = np.array(new_k_list)

    # if there are no symmetry operations,
    # than the weight on each subgrid k point is exactly 1/N of the original weight
    if symmetries[0][0] == symmetries[0][1] == symmetries[0][2] == 0:
        return new_k_list.tolist()

    # defining a matrix to save the degeneracies, after each symmetry operation
    check_degeneracy = np.zeros((number_elements, number_elements), dtype=bool)
    # saving flag if no degeneracy detected
    no_degeneracy = True

    # weights given to the k list points
    k_list_weights=np.zeros(number_elements,type=float)

    # each k point is compared to each other k point after the symmetry operation has been applied
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
    # if no degeneracy is detected, the common result is given
    if no_degeneracy:
        return new_k_list
    else:
        # if degeneracy is detected, of the degenerate points only one is chosen
        
        # assuming all k points as not-degenerate, and creating a dictionary
        # where to each eigenspace is associated an eigenvector (k_point)
        list = dict([(str(i), [i]) for i in range(number_elements)])

        for i in range(0, number_elements - 1):
            for j in range(i + 1, number_elements):
                if len(list) != 1:
                    # reading the degeneracy matrix, the eigenspaces are properly populated,
                    # in particular degenerate k points are put in the same eigenspace
                    if check_degeneracy[i][j] == True:
                        for key, values in list.items():
                            for value in values:
                                if value == i:
                                    positioni = key
                                if value == j:
                                    positionj = key
                        new_list = {}
                        # new list is created with the proper eigenspaces, but at the same
                        # time is compared with the old one, to check if there are old
                        # degeneracies to take into account
                        if positioni != positionj:
                            count = 0
                            for key, values in list.items():
                                if key == positioni:
                                    new_list[str(count)] = values + list[positionj]
                                    count = count + 1
                                else:
                                    if key != positionj:
                                        new_list[str(count)] = values
                                        count = count + 1
                        # The new list is save to be used as an old one
                        list = {}
                        for key, values in new_list.items():
                            list[str(key)] = new_list[str(key)]

        # one of the equivalent k points is chosen 
        # the weight of the eigenspace is given to the chosen k point
        for key, values in new_list.items():
            new_list[str(key)] = sorted(set(new_list[str(key)]))

        if len(new_list) != 0:
            weight_eigenspace = k_origin_weight / len(new_list)
        else:
            weight_eigenspace = 0.0
    
        for key, values in new_list.items():
            # if the eigenspace has only one k point, the weight of the eigenspace is givent to the k point itself
            if len(values) == 1:
                k_list_weights[values] = weight_eigenspace
            # if there are more k points one is chosen, and the weight of the eigenspace is givent to the chosen k point=
            else:
                # chosing of the degenerate points only the first one
                k_list_weights[values[0]] = weight_eigenspace
        
        new_k_list = np.c_[ k_list, k_list_weights]
        # returning the matrix with the respective weights
        return new_k_list

# local refinment of the k list in the 2d plane (the generalization to a 3d parallelepiped is quite straighforward)
def local_refinment(
    new_k_list,
    old_k_list,
    refinment_spacing,
    refinment_iteration,
    symmetry,
    threshold,
    brillouin_primitive_vectors_3d,
    selected_2d_plane,
    normalized_brillouin_primitive_vectors_2d
):
    r"""
    starting from a set of k points (old_k_list) with certain weights, to each k point iteratively (refinment_iteration) a new subset of k points is associated,
    a symmetry analysis (symmetry) is performed on the subset, and consequently the initial weight of the k point is distributed in the subset of k points
    this procedure goes on for the number of refinment iterations 
    at each iteration for each k point 4 new points are generated along the reciprocal primitive vectors of the selcted 2d plane, at a distance from the original k point 
    equal to the refinment spacing (the refinment spacing is halved at each iteration)

    obviously no refinment is aplied to k points with null weight
    Parameters
    ----------
    new_k_list: (,4) |array_like| list of values produced by the refinment procedure
    old_k_list: (,4) |array_like| list of values give to the refinment procedure (kx,ky,kz,weight)
    refinment_spacing: |double| initial spacing given to generate the refinment, half of the preceding refinment spacing is considered at each refinment iteration
    refinment_iteration: |int| number of refinment iterations considered
    symmetry: |list of lists| a symmetry is a list of 3 elements (the versor is the axis of rotation, while the modulus is the angle)
    threshold: |double| therhold to recognize a symmetry
    normalized_brillouin_primitive_vectors_2d : (2,3) |array_like| normalized vectors definining the 2d plane in the reciprocal space, chosen for the refinment procedure
    brillouin_primitive_vectors_3d (3,3) |array_like| :list columns: kx,ky,kz, rows: b1,b2,b3 
    """

    length_old_k_list=old_k_list.shape[1]
    if refinment_iteration == 0:
        for r in range(length_old_k_list):
            if old_k_list[i,3]!=0:
                new_k_list.append(old_k_list[i,:])
                
    else:
        iter = 0
        while iter != length_old_k_list:
            # reading the k points inside the list and refine around them
            if length_old_k_list == 1:
                k_tmp = old_k_list[:3]
                k_tmp_weight = old_k_list[3]
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
                    symmetry,
                    threshold
                )
                
                k_tmp_subgrid=k_tmp_subgrid[np.all(k_tmp_subgrid[:,3]!=0, axis=1),:]

                k_tmp_subgrid[:,:3] = check_inside_brillouin_zone(
                    k_tmp_subgrid[:,:3],
                    brillouin_primitive_vectors_3d,
                    selected_2d_plane,
                )
                
                # applying to the remaining points of the subgrid the refinment procedure
                new_refinment_spacing = refinment_spacing / 2
                local_refinment(
                    new_k_list,
                    k_tmp_subgrid,
                    new_refinment_spacing,
                    refinment_iteration - 1,
                    symmetry,
                    threshold,
                    brillouin_primitive_vectors_3d,
                    selected_2d_plane,
                    normalized_brillouin_primitive_vectors_2d  
                )

def interpolation_k_points_weights_in_2D(
    k_points_list_with_weights,
    brillouin_primitive_vectors_2d,
    number_elements
    ):
    r"""
    through an interpolation, the given list of k points defined in a 2d plane (brillouin_primitive_vectors_2d) is used to produce a list of ordered k points, properly weighted
    the number of k points produced is equal to number_elements x number_elements
    Parameters 
    ----------
    k_points_list_with_weights : (,4) |array_like| : list of k points (kx,ky,kz,w) where w is the respective weight of the k point
    brillouin_primitive_vectors_2d: (2,3) |array_like| chosen 2D plane in the reciprocal space (the list of k points has to be generated from these reciprocal vectors, i.e. be in the 2D plane generated by these vectors)
    
    Return
    -------------
    new_k_points_list_with_weights_and_ordering : (:6) |matrix| : list of k points (kx,ky,kz,w,i,j) where w is the respective weight of the k point, and (i,j) the indices mapping the 2D grid 
    note: the weights, being added some points in the interpolation, are renormalized over the entire set of points (it is assumed that the initial list is already normalized)
    """
    # considering first a mapping from the list of k vectors in 3D with the respective weight, to a list of k vectors in 2D with the respective weight
    # the projection of the k vectors on the primitive reciprocal vectors is therefore considered, instead of the cartesian coordinates
    
    number_k_points=k_points_list_with_weights.shape[1]
    k_vectors=k_points_list_with_weights[:,:3]
    k_vectors_projections=k_vectors @ brillouin_primitive_vectors 
    
    #this pairs (the projection coordinates) are localized in a box [0,1)x[0,1), making an interpolation quite straightforward and reliable
    xmin,xmax=0,1
    ymin,ymax=0,1
    
    #number of points considered in the interpolation (a differentiation between the two directions can be done)
    if number_elements is None:
        ny,nx=number_k_points,number_k_points
        number_elements=number_k_points
    else:
        ny,nx=number_elements,number_elements

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
    new_k_points_list_with_weights_and_indices=np.zeros((number_elements,number_elements,4))
    coordinates=np.zeros(3)
    for i in range(number_elements):
        for j in range(number_elements):
            vector=np.asanyarray([xi[i,j],yi[i,j]])
            coordinates=(brillouin_primitive_vectors_2d.T)@vector
            new_k_points_list_with_weights_and_indices[i][j][:3]=coordinates
            new_k_points_list_with_weights_and_indices[i][j][3]=zi[i][j]

    return(new_k_points_list_with_weights_and_indices,number_elements,number_elements)

# generation of a 2D k points grid, through the refinment procedure
# due to the disordering of the k points by the refinment procedure, 
# an interpolation scheme has been added to obtain a propelry ordered k points grid
def k_points_grid_generator_2D(
    brillouin_primitive_vectors_3d,
    selected_plane_2d,
    initial_grid_spacing,
    shift_in_plane,
    shift_in_space,
    symmetry,
    refinment_spacing,
    refinment_iteration,
    threshold,
    number_elements
):
    r"""
    k points grid 2D
    Parameters
    ----------
    brillouin_primitive_vectors_3d (3,3) |array_like| :list columns: kx,ky,kz, rows: b1,b2,b3 
    selected_plane_2d : (,3) |array_like|_  ex. [0,1,1] the plane is the b2 x b3 in reciprocal space where the bi are the primitive vectors
    initial_grid_spacing: |double| spacing of the k grid before the refinment procedure
    shift_in_plane : (,2) |array_like|_  shift of the selected reciprocal 2D plane with respecet to the crystal coordinates
    shift_in_space : (,3) |array_like|_  shift of the selected reciprocal 2D plane with respecet to the cartesian coordinates
    symmetry: |list of lists| a symmetry is a list of 3 elements (the versor is the axis of rotation, while the modulus is the angle)
    refinment_spacing: |double| initial spacing given to generate the refinment, half of the preceding refinment spacing is considered at each refinment iteration
    refinment_iteration: |int| number of refinment iterations considered
    threshold: |double| therhold to recognize a symmetry
    Return
    -------------
    normalized_brillouin_primitive_vectors_2d: (2,3) |matrix| versor of the chosen 2D reciprocal plane
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
        if selected_plane_2d[i] != 0:
            brillouin_primitive_vectors_2d[count] = brillouin_primitive_vectors_3d[i]
            normalized_brillouin_primitive_vectors_2d[count] = brillouin_primitive_vectors_2d[count]/(brillouin_primitive_vectors_2d[count]@brillouin_primitive_vectors_2d[count])
            count += 1

    k0 = int(
        (brillouin_primitive_vectors_2d[0] @ brillouin_primitive_vectors_2d[0])/ initial_grid_spacing
    )
    k1 = int(
        (brillouin_primitive_vectors_2d[1] @ brillouin_primitive_vectors_2d[1]) / initial_grid_spacing
    )
    initial_weights = 1 / (k0 * k1)
    
    # building the starting 2d k points grid
    k_points_grid_2d = []
    k_point_tmp = np.zeros(4)
    for i in range(k0):
        for j in range(k1):
            k_point_tmp[:3] = (float(i + shift_in_plane[0]) / k0) * (
                shift_in_space + brillouin_primitive_vectors_2d[0]
            ) + (float(j + shift_in_plane[1]) / k1) * (
                shift_in_space + brillouin_primitive_vectors_2d[1]
            )
            k_point_tmp[3] = initial_weights
            k_points_grid_2d.extend([k_point_tmp])
    
    print(shift_in_plane)
    print(k_points_grid_2d)
    ##for i in range(len(k_points_grid_2d)):
      ##  print(k_points_grid_2d[i,:])
  
   ## # applying the refinment procedure to the 2d k points grid
   ## refined_grid = []
   ## local_refinment(
   ##     refined_grid,
   ##     k_points_grid_2d,
   ##     refinment_spacing,
   ##     refinment_iteration,
   ##     symmetry,
   ##     threshold,
   ##     brillouin_primitive_vectors_3d,
   ##     selected_plane_2d,
   ##     normalized_brillouin_primitive_vectors_2d  
   ## )     
   ## refined_grid = np.reshape(
   ##     refined_grid, len(refined_grid)/4, 4)
   ## 
   ## interpolating and ordering the 2d k points
   ## new_list_of_k_points,n0,n0=interpolation_k_points_weights_in_2D(refined_grid,brillouin_primitive_vectors_2d,number_elements)
    
    return (normalized_brillouin_primitive_vectors_2d,k0,k1,k_points_grid_2d)

# function considering one point in a list of k points, a refinment procedure is applied to the k point selected
# if a radius is indicated, all the near by ones k points are selected and their weight redistributed over a denser grid
#(the second option is quite straightforward to implement...)
def dynamical_refinment(
    brillouin_primitive_vectors_3d,
    selected_plane_2d,
    brillouin_primitive_vectors_2d,
    normalized_brillouin_primitive_vectors_2d,
    symmetry,
    refinment_iteration,
    threshold,
    k_point_selected,
    old_k_list,
    ordering,
    number_elements
):
    # position of the k point selected in the list
    position_k_point = [i for i,val in enumerate(old_k_list) if val == k_point_selected]
    del old_k_list[position_k_point]
    
    minimal_distance = np.min(np.abs(old_k_list[:,:3]-k_point_selected[:3]))
    refinment_spacing = minimal_distance
    # applying the refinment procedure to the k point selected
    refined_grid = []
    local_refinment(
        refined_grid,
        k_point_selected,
        refinment_spacing,
        refinment_iteration,
        symmetry,
        threshold,
        brillouin_primitive_vectors_3d,
        selected_plane_2d,
        normalized_brillouin_primitive_vectors_2d  
    )     
    refined_grid = np.reshape(
        refined_grid, len(refined_grid)/4, 4)

    # adding the new points to the old list
    old_k_list.append(refined_grid)

    if ordering == False:
        return old_k_list
    else:
    # interpolating and ordering the 2d k points
        new_list_of_k_points,n0,n0=interpolation_k_points_weights_in_2D(refined_grid,brillouin_primitive_vectors_2d,number_elements)
        return new_list_of_k_points,n0,n0

#TESTING INPUT
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

    brillouin_primitive_vectors_3d=np.zeros((3,3))
    brillouin_primitive_vectors_3d[0]=[1,0,0]
    brillouin_primitive_vectors_3d[1]=[0,1,0]
    brillouin_primitive_vectors_3d[2]=[0,0,1]
    selected_plane_2d=[1,1,0]
    #threshold for understanding if there is a degeneracy
    threshold=0.0000001
    ##shift_in_plane=[0,0]
    ###shift_in_space=[0,0,0]
    symmetry=[[0,0,0]]
    initial_grid_spacing=0.1
    refinment_iteration=0
    ### refinment_spacing=0.05

    
    normalized_brillouin_primitive_vectors_2d,n0,n0,new_list_of_k_points=k_points_grid_generator_2D(
        brillouin_primitive_vectors_3d,
        selected_plane_2d,
        initial_grid_spacing,
        None,
        None,
        symmetry,
        None,
        refinment_iteration,
        threshold,
        None
    )

   ## for i in range(len(new_list_of_k_points)):
   ##     print(new_list_of_k_points[i])

    file1="k_point_normal_grid.txt"
    with open(file1, 'w') as file1:
        for i in range(n0):
            for j in range(n0):
                for r in range(len(new_list_of_k_points)):
                    file1.write(str(new_list_of_k_points[r,0])+' '+str(new_list_of_k_points[r,1])+' '+str(new_list_of_k_points[r,2])+' '+str(new_list_of_k_points[r,3])+'\n')

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

#function checking if the k points in the list k_list are inside the Brillouin zone
def downfold_inside_brillouin_zone(
    k_point_list,
    brillouin_primitive_vectors_3d
):
    r"""
        check if any point in a list is inside the first brillouin zone, in case the point is outside it is shifted back to the first brillouin zone
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
    flag_inside_out=np.zeros((number_k_point_elements),dtype=int)

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
                flag_inside_out[i]=1
                for d in range(3):
                    k_point_list[i,d]-=int(transformed_k_point_list[i,d])*brillouin_primitive_vectors_3d[r,d]
    return k_point_list, flag_inside_out

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
                    max_value=max(i,j)
                    k_eigenspaces[k_eigenspaces==max_value]=min_value

        # counting number of different eigenspaces
        number_eigenspaces=len(np.unique(k_eigenspaces))
        weights=k_origin_weight/number_eigenspaces

        # ordering eigenspaces values
        ordered_k_eigenspaces=np.zeros(number_elements)
        counting=0
        assigned_indices=[]
        for i in range(number_elements):
            if k_eigenspaces[i] not in assigned_indices:
                if counting<number_eigenspaces:
                    ordered_k_eigenspaces[k_eigenspaces==k_eigenspaces[i]]=counting
                    assigned_indices.extend([r for r in range(number_elements) if k_eigenspaces[r]==k_eigenspaces[i]])
                    counting+=1
        # considering a representative for each eigenspace
        print(k_eigenspaces)
        for i in range(number_eigenspaces):
            list_positions=[r for r in range(number_elements) if ordered_k_eigenspaces[r]==i]
            new_k_list[list_positions[0],3]=weights
            new_k_list[list_positions[1:],3]=0

        # returning the matrix with the respective weights
        return new_k_list

# local refinment of the k list in the plane (parallelepiped) pointed out by the brillouin primitive vectors given (it can be 2D, 3D or ..)
def local_refinment_with_symmetry_analysis(
    new_k_list,
    old_k_list,
    refinment_spacing,
    refinment_iteration,
    symmetries,
    threshold,
    brillouin_primitive_vectors_3d,
    brillouin_primitive_vectors_2d,
    normalized_brillouin_primitive_vectors_2d,
    downfold
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
                new_k_list.extend(old_k_list[i,:])
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
                
                if downfold==True:
                    k_tmp_subgrid[:,:3],flag_inside_out= downfold_inside_brillouin_zone(
                        k_tmp_subgrid[:,:3],
                        brillouin_primitive_vectors_3d,
                    )

                # applying to the points of the subgrid to the refinment procedure
                new_refinment_spacing = refinment_spacing / 2
                
                local_refinment_with_symmetry_analysis(
                    new_k_list,
                    k_tmp_subgrid,
                    new_refinment_spacing,
                    refinment_iteration-1,
                    symmetries,
                    threshold,
                    brillouin_primitive_vectors_3d,
                    brillouin_primitive_vectors_2d,
                    normalized_brillouin_primitive_vectors_2d,
                    downfold
                )

def k_points_triangulation_2d(
    k_points_list,
    brillouin_primitive_vectors_2d,
    count
):
    r"""
    the sparse k points (kx,ky,kz,w) are linked togheter through a triangulation procedure (Delaunay procedure)
    the so-built triangles are listed, each triangle is characterized by three indices pointing to the vertices position in the k poins list  
    the ordering of the vertices is so to assure a cloack-wise path along the triangles

    the periodicity of the BZ is properly taken into account:
        the points at the left border or bottom border are considered twice in the triangulation
        once on their proper border, then on the opposite border (translation of a reciprocal vector)
        the index of the triangulation is the same both times, this assures the periodicity condition
    
    Parameters
    ---------
    k_points_list: (n,4) |array_like| (kx,ky,kz,w)
    brillouin_primitive_vectors_2d: (2,3) |array_like|
    count: (int) taking into account if the k points need to be added to a list, so numbering properly the vertices, in order to be able to use these triangles with the new list
    border_count: (int) taking into account if the first k points in k_points_list are from the border and are not neeeded to be added to the list  

    Returns
    ---------
    k_points_list: (n,4) |array_like| (kx,ky,kz,w)
    triangles: (s,3) |array_like| 
    left_border_points, bottom_border_points indices of the points at the border of the BZ
    """

    number_elements=k_points_list.shape[0]
    k_points_list_projections=np.zeros((number_elements,2),dtype=float)

    #choosing one point as an origin in the 2d brillouin plane
    origin=k_points_list[0,:3]
    #calculating the projections of the different k points on the primitive vectors
    vector_0=np.zeros((number_elements,2),dtype=float)
    vector_1=np.zeros((number_elements,2),dtype=float)
    for i in range(number_elements):
        vector_0[i,0]=1
        vector_1[i,1]=1
        for j in range(2):
            k_points_list_projections[i,j]=(k_points_list[i,:3]-origin)@brillouin_primitive_vectors_2d[j,:]
    
    #considering the periodic images
    k_points_list_projections_periodic_images_0=np.zeros((number_elements,2),dtype=float)
    k_points_list_projections_periodic_images_1=np.zeros((number_elements,2),dtype=float)
    k_points_list_projections_periodic_images_0=k_points_list_projections+vector_0
    k_points_list_projections_periodic_images_1=k_points_list_projections+vector_1

    #triangulation of the set of points
    old_triangles_tmp = Delaunay(k_points_list_projections).simplices
    old_triangles=[]
    for i in range(len(old_triangles_tmp)):
        old_triangles.append([[element,0] for element in old_triangles_tmp[i]])

    #to consider the border is equivalent to consider the trianglues between k_points_list_projections_periodic_images_0/1 and k_points_list_projections
    old_triangles_border_1=[]
    k_points_list_projections_periodic_images_1=np.vstack([k_points_list_projections,k_points_list_projections_periodic_images_1])
    old_triangles_border_1_tmp = Delaunay(k_points_list_projections_periodic_images_1).simplices
    for i in range(len(old_triangles_border_1_tmp)):
        if np.any(old_triangles_border_1_tmp[i]<number_elements) and  np.any(old_triangles_border_1_tmp[i]>number_elements):
            old_triangles_border_1.append([[element,0] if element<number_elements else [element-number_elements,2] for element in old_triangles_border_1_tmp[i]])
    old_triangles_border_0=[]
    k_points_list_projections_periodic_images_0=np.vstack([k_points_list_projections,k_points_list_projections_periodic_images_0])
    old_triangles_border_0_tmp = Delaunay(k_points_list_projections_periodic_images_0).simplices
    for i in range(len(old_triangles_border_0_tmp)):
        if np.any(old_triangles_border_0_tmp[i]<number_elements) and  np.any(old_triangles_border_0_tmp[i]>number_elements):
            old_triangles_border_0.append([[element,0] if element<number_elements else [element-number_elements,1] for element in old_triangles_border_0_tmp[i]])

    #concatenating the triangles totally inside the FBZ (0) and the ones going through the border (1/2)
    old_triangles=np.vstack([old_triangles,old_triangles_border_0])
    old_triangles=np.vstack([old_triangles,old_triangles_border_1])
    #these are the different triangles, where each element has two numbers representing the respective element of the list and one element distinguishing if the point is passing trough the border of the FBZ
    
    #the projections along the two axis can be used to properly order the vertices of the triangles
    new_triangles = []
    vertices_projections = np.zeros((3,2))
    number_of_triangles = len(old_triangles) 
    for i in range(number_of_triangles):
        vertices=old_triangles[i,:]
        one=[1,0]
        two=[0,1]
        for r in range(3):
            if vertices[r,1]==0:
                vertices_projections[r]=k_points_list_projections[vertices[r,0]]
            elif vertices[r,1]==1:
                vertices_projections[r]=k_points_list_projections[vertices[r,0]]+one
            else:
                vertices_projections[r]=k_points_list_projections[vertices[r,0]]+two
        indices=[tup[0] for tup in sorted(enumerate(vertices_projections[:,1]))]
        vertices_projections=vertices_projections[indices]
        if vertices_projections[0,1]==np.min(vertices_projections[:,0]):
            indices=[tup[0] for tup in sorted(enumerate(vertices_projections[:,0]))]
            vertices_projections=vertices_projections[indices]
        vertices=vertices[indices]
        #considering the case the function is called on an existing list (the existing list has length equal to border_count)
        for s in range(3):
            if  vertices[s,0]> count:
                vertices[s,0]+=count
        new_triangles.append(vertices)
       
    #for each triangle the ordering of the vertices is now clock-wise
    return k_points_list,new_triangles

#generation of a 2D k points grid 
def k_points_generator_2D(
    brillouin_primitive_vectors_3d,
    initial_grid_spacing,
    chosen_plane,
    brillouin_primitive_vectors_2d,
    symmetry_analysis_flag,
    refinment_spacing=0,
    symmetries=[[0,0,0]],
    threshold=0.0001,
    refinment_iteration=0,
    shift_in_plane=[0,0],
    shift_in_space=[0,0,0],
    count=0 
):
    r"""
    k points generator in a 2d plane
    
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
    symmetry_analysis: |boolean| local_refinment_with_symmetry_analysis or not
    Return
    -------------
    "local_refinment_with_symmetry_analysis" not considered:
        k_points_list: (1,k0*k1) |array_like| k points (kx,ky,kz,w) where w is the weight of the k point
        parallelograms : (n,4) |array_like| each row is a parallelogram, where the 4 vertices are integers, which correspond to the list positions of the given k points 
    "local_refinment_with_symmetry_analysis" considered:
        k_points_list: (1,k0*k1) |array_like| k points (kx,ky,kz,w) where w is the weight of the k point
        triangles: (n,3) |array_like| each row is a triangle, where the 3 vertices are integers, which correspond to the list positions of the given k points 
    """
    normalized_brillouin_primitive_vectors_2d = np.zeros((2, 3),dtype=float)

    if  chosen_plane==[0,0,0]:
        chosen_plane=[]
        for i in range(2):
            normalized_brillouin_primitive_vectors_2d[i] = brillouin_primitive_vectors_2d[i]/(brillouin_primitive_vectors_2d[i]@brillouin_primitive_vectors_2d[i])
            # redefining the brillouin_primitive_vectors_3d in order to take into account the input brillouin_primitive_vectors_2d
            for j in range(3):
                if (brillouin_primitive_vectors_3d[j,:]@brillouin_primitive_vectors_2d[i,:])!=0:
                    chosen_plane.append(j)
        count = 0
        for i in chosen_plane:
            brillouin_primitive_vectors_3d[i,:]=brillouin_primitive_vectors_2d[count,:]
            count +=1                    
    else:
        #saving the chosen 2d plane 
        count = 0
        for i in range(3):
            if chosen_plane[i] != 0:
                brillouin_primitive_vectors_2d[count] = brillouin_primitive_vectors_3d[i]
                normalized_brillouin_primitive_vectors_2d[count] = brillouin_primitive_vectors_2d[count]/(brillouin_primitive_vectors_2d[count]@brillouin_primitive_vectors_2d[count])
                count += 1

    #initial 2D grid dimensions
    k0 = int(
        (brillouin_primitive_vectors_2d[0] @ brillouin_primitive_vectors_2d[0])/ initial_grid_spacing
        )
    k1 = int(
        (brillouin_primitive_vectors_2d[1] @ brillouin_primitive_vectors_2d[1]) / initial_grid_spacing
        )
    initial_weights = 1 / (k0 * k1)
    #building the initial 2D k points grid
    k_points_grid = np.zeros((k0*k1,4),dtype=float)
    for i in range(k0):
        for j in range(k1):
            k_points_grid[i*k1+j,:3] = (float(i + shift_in_plane[0]) / k0) * (shift_in_space + brillouin_primitive_vectors_2d[0,:]) + (float(j + shift_in_plane[1]) / k1) * (shift_in_space + brillouin_primitive_vectors_2d[1,:])
            k_points_grid[i*k1+j,3] = initial_weights
    
    if symmetry_analysis_flag == True:
        if refinment_spacing==0:
            refinment_spacing=initial_grid_spacing/2

        #applying the refinment procedure to the 2d k points grid
        #obtaining a refined list
        refined_k_points_list = []
        local_refinment_with_symmetry_analysis(
            refined_k_points_list,
            k_points_grid,
            refinment_spacing,
            refinment_iteration,
            symmetries,
            threshold,
            brillouin_primitive_vectors_3d,
            brillouin_primitive_vectors_2d,
            normalized_brillouin_primitive_vectors_2d,
            True
        )
        #transforming a list into an array-like element
        k0=int(len(refined_k_points_list)/4)
        refined_k_points_list = np.reshape(refined_k_points_list,(int(len(refined_k_points_list)/4), 4))
        #applying a triangulation to order the k points
        #the periodicity of the BZ is properly considered
        refined_k_points_list,triangles=k_points_triangulation_2d(refined_k_points_list,brillouin_primitive_vectors_2d,0)
        return refined_k_points_list,triangles
    else:
        not_refined_k_points_list = np.zeros((k0*k1,4))
        parallelograms=np.zeros((k0*k1,4))
        left_border_points=[]
        bottom_border_points=[]
        for i in range(k0):
            for j in range(k1):
                not_refined_k_points_list[i*k1+j][:3] = (float(i + shift_in_plane[0]) / k0) * (shift_in_space + brillouin_primitive_vectors_2d[0]) \
                   + (float(j + shift_in_plane[1]) / k1) * (shift_in_space + brillouin_primitive_vectors_2d[1])
                not_refined_k_points_list[i*k1+j][3] = initial_weights
                #periodicity of the BZ
                if i<k0-1 and j<k1-1:
                    parallelograms[i*k1+j]=[i*k1+j,(i+1)*k1+j,(i+1)*k1+j+1,i*k1+j+1]
                elif i<k0-1 and j==k1-1:
                    left_border_points.append([i*k1])
                    parallelograms[i*k1+j]=[i*k1+j,(i+1)*k1+j,(i+1)*k1+0,i*k1+0]
                elif i==k0-1 and j<k1-1:
                    bottom_border_points.append([j])
                    parallelograms[i*k1+j]=[i*k1+j,0*k1+j,0*k1+j+1,i*k1+j+1]
                else:
                    parallelograms[i*k1+j]=[i*k1+j,0*k1+j,0*k1+0,i*k1+0]
                parallelograms[i*k1+j]=[[count]*4]+parallelograms[i*k1+j]

        return not_refined_k_points_list,parallelograms,left_border_points,bottom_border_points,k0,k1

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
    k_point_list: (k-s,3) |array_like| the k points that resulted to be inside one of the triangles
    """
    
    number_elements=k_point_list.shape[0]
    number_triangles=subset_triangles.shape[0]

    # for each triangle the vecotrial products between the different vertices, and between the vertices and the to-be-checked k points
    # are used to calculate the barycentric weights of the to-be-checked k points, which give information about the to-be-checked k point being inside or outside the triangle itself
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

def dynamical_refinment_little_paths_2d(
        little_paths,
        position_littles_paths_to_refine,
        number_vertices,
        all_k_points_list,
        refinment_iteration,
        symmetries,
        threshold,
        brillouin_primitive_vectors_3d,
        brillouin_primitive_vectors_2d,
        normalized_brillouin_primitive_vectors_2d
):
    r"""
    some k points of the "all k points list" are selected (the ones in the selected little paths "position_littles_paths_to_refine"), 
    and a local refinment is applied to them, in case of triangles (number_vertices=3) a refinment with symmetry analysis is applied,
    in case of parallelograms (number_vertices=4) a refinment without symmetry analysis is applied 
    points on the border or nearby the border of the BZ are properly considered

    a refinmente without symmetry analysis is done in the case of parallelograms 

    Parameters
    ---------
    little_paths: (s,number_vertices) |array_like| (v1,v2,v3...) int (v=vertex)
    all_k_points_list: (n,4) |array_like| (kx,ky,kz,w)
    positions_little_paths_to_refine: |list| (int) positions of the little paths to which the refinment procedure is applied
    number_vertices: |int| distinguishing between traingles and parallelograms

    (parameters for the refinment...)

    Returns
    ---------
    new_k_points_list: (n,4) |array_like| (kx,ky,kz,w)
    added_little_paths: (,3) |array_like| (v1,v2,v3)
    """

    number_little_paths=len(little_paths)
    # canceling any repetition
    position_littles_paths_to_refine=list(np.unique(np.asarray(position_littles_paths_to_refine)))
    number_littles_paths_to_refine=len(position_littles_paths_to_refine)
    # considering the selected little paths
    little_paths_to_refine = little_paths[position_littles_paths_to_refine]
    # in the list of all little paths, the ones selected are substituted with a little path with vertices equal to -1
    for i in position_littles_paths_to_refine:
        little_paths[i]=[-1 for j in range(number_vertices)]
    # the little paths around the to-refine little paths are associated to each to-refine little path
    little_paths_around=[]
    for i in range(number_littles_paths_to_refine):
        little_paths_around.append(list(set([ [int(j),i] for r in range(number_little_paths) for j in little_paths[r] if j in little_paths_to_refine[i]])))
    
    FINIRE DA QUI

    print(little_paths_around)
    # checking if two to-refine little paths are nearby (they have one side in common, i.e. one of the little_paths_around is equal to position_littles_paths_to_refine)
    check_degeneracy = np.zeros((number_littles_paths_to_refine, number_littles_paths_to_refine), dtype=bool)
    any_degeneracy = 0
    for i in range(number_littles_paths_to_refine-1):
        for j in range(i,number_littles_paths_to_refine):
            for r in little_paths_around[i]:
                if r in little_paths_around[j]:
                    check_degeneracy[i,j]=True
                    any_degeneracy+=1
                else:
                    check_degeneracy[i,j]=False
    if any_degeneracy != 0:
        for i in range(number_littles_paths_to_refine-1):
            for j in range(i,number_littles_paths_to_refine):
                if check_degeneracy[i][j]==True:
                    min_value=min(i,j)
                    max_value=max(i,j)
                    little_paths_around[min_value].extend(little_paths_around[max_value])
                    
        
        print(little_paths_around)
        print(check_degeneracy)
        # counting number of new different little paths
        number_littles_paths_to_refine=len(np.unique(little_paths_around ))

        ### ordering eigenspaces values
        ##ordered_k_eigenspaces=np.zeros(number_elements)
        ##counting=0
        ##assigned_indices=[]
        ##for i in range(number_elements):
        ##    if k_eigenspaces[i] not in assigned_indices:
        ##        if counting<number_eigenspaces:
        ##            ordered_k_eigenspaces[k_eigenspaces==k_eigenspaces[i]]=counting
        ##            assigned_indices.extend([r for r in range(number_elements) if k_eigenspaces[r]==k_eigenspaces[i]])
        ##            counting+=1
        
    # checking if two to refine-little paths have a little path around in common         



    #print(little_paths_around)                            
    # the k points of the little paths around are associated to each to-refine little path
    border_k_points={}
    for key,element in little_paths_around.items():
        border_k_points[key]=border_k_points[key]+[i for i in little_paths[element] if i!=key]
    # procede to a refinement of the k points constituing the to-refine little paths
    if number_vertices==3:
        all_new_triangles=[]
        # finding the minimum distance between the to-refine k points and the border k points
        for i in range(number_littles_paths_to_refine):
            distances=np.zeros(number_vertices*len(border_k_points[i]))
            for j in range(number_vertices):
                for s in border_k_points[i]:
                    distances[j*3+s]=np.linalg.norm(all_k_points_list[little_paths_to_refine[i,j]]-all_k_points_list[s])
            minimal_distance=np.min(distances)
            refined_k_points_list=[]
            local_refinment_with_symmetry_analysis(
                refined_k_points_list,
                all_k_points_list[little_paths_to_refine[i,:]],
                minimal_distance,
                refinment_iteration,
                symmetries,
                threshold,
                brillouin_primitive_vectors_3d,
                brillouin_primitive_vectors_2d,
                normalized_brillouin_primitive_vectors_2d,
                False
            )
            #transforming from a list to array_like elements
            refined_k_points_list = np.reshape(refined_k_points_list,(int(len(refined_k_points_list)/4), 4))
            #check if the points are inside the border 
            refined_k_points_list=check_inside_closed_shape(
                np.concatenate(little_paths_around[i],little_paths_to_refine[i]),
                refined_k_points_list,
                all_k_points_list)
            #triangulation
            border_plus_refined_k_points_list=np.concatenate(all_k_points_list[border_k_points[i]],refined_k_points_list)
            border_plus_refined_k_points_list,added_triangles=k_points_triangulation_2d(
                                                                    border_plus_refined_k_points_list,
                                                                    brillouin_primitive_vectors_2d,
                                                                    count,
                                                                    len(border_k_points[i]))
            all_new_triangles.append(added_triangles)
            return refined_k_points_list, added_triangles
           
        ## possono stare al bordo i punti da raffinare, o uno dei punti del bordo dei punti da raffinare...
    else:
        number_border_k_points=len(border_k_points)
        projections_border_k_points=np.zeros((border_k_points,2))
        origin=all_k_points_list[border_k_points[0]]
        #calculating the projections of the different k points on the primitive vectors
        for i in range(number_border_k_points):
            projections_border_k_points[i]=np.dot(all_k_points_list[border_k_points[i]]-origin,brillouin_primitive_vectors_2d)
        #calculatting the origin-point
        origin=np.argmin(np.asarray((list(map(lambda x,y: np.sqrt(x*x+y*y), projections_border_k_points[:,0],projections_border_k_points[:,1])))))
        #calculating the projections of the different k points on the primitive vectors (with the new origin)
        for i in range(number_border_k_points):
            projections_border_k_points[i]=np.dot(all_k_points_list[border_k_points[i]]-all_k_points_list[border_k_points[origin]],brillouin_primitive_vectors_2d)
        elementy=sorted(projections_border_k_points,key=lambda x: x[0])[-1]
        elementx=sorted(projections_border_k_points,key=lambda x: x[1])[-1]
        for i in range(number_border_k_points):
            if elementx == projections_border_k_points[i]:
                x_i=i
            if elementy == projections_border_k_points[i]:
                y_i=i
        #defining the primitive vectors of the subspace
        little_brillouin_primitive_vectors_2d=np.zeros((2,3))
        little_brillouin_primitive_vectors_2d[0,:]=all_k_points_list[border_k_points[x_i]]-all_k_points_list[border_k_points[origin]]
        little_brillouin_primitive_vectors_2d[1,:]=all_k_points_list[border_k_points[y_i]]-all_k_points_list[border_k_points[origin]]
        not_refined_k_points_list,parallelograms,_bt,_bl=k_points_generator_2D(
            brillouin_primitive_vectors_3d,
            None,
            None,
            little_brillouin_primitive_vectors_2d,
            None,
            None,
            None,
            None,
            None,
            len(all_k_points_list))
        ###NECESSARIO STUDIARE IL CASO AL BORDO ...
        return not_refined_k_points_list, parallelograms


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

    brillouin_primitive_vectors_3d=np.zeros((3,3),dtype=float)
    brillouin_primitive_vectors_3d[0]=[2,0,0]
    brillouin_primitive_vectors_3d[1]=[0,2,0]
    brillouin_primitive_vectors_3d[2]=[0,0,1]
    brillouin_primitive_vectors_2d=np.zeros((2,3),dtype=float)
    chosen_plane=[1,1,0]
    symmetries=[[0,0,np.pi/2]]
    initial_grid_spacing=0.1
    refinment_iteration=5
    refinment_spacing=0
    threshold=0.001
    count=0
    shift_in_plane=[1,0]
    shift_in_space=[0,0,0]
   
    not_refined_k_points_list,parallelograms,left_border_points,bottom_border_points,k0,k1=k_points_generator_2D(
        brillouin_primitive_vectors_3d,
        initial_grid_spacing,
        chosen_plane,
        brillouin_primitive_vectors_2d,
        False,
        refinment_spacing,
        symmetries,
        threshold,
        refinment_iteration,
        shift_in_space,
        shift_in_space,
        count
        )
    print(not_refined_k_points_list)
    print(parallelograms)
    indices=[0,1,2,3,4,4,4,5]
    subset_parallelograms=parallelograms[indices]
    print(subset_parallelograms)

    normalized_brillouin_primitive_vectors_2d=np.zeros((2,3))
    for i in range(2):
        for r in range(3):
            normalized_brillouin_primitive_vectors_2d[i,r]=brillouin_primitive_vectors_2d[i,r]/np.dot(brillouin_primitive_vectors_2d[i],brillouin_primitive_vectors_2d[i])
    
    not_refined_k_points_list, parallelograms=dynamical_refinment_little_paths_2d(
        parallelograms,
        indices,
        4,
        not_refined_k_points_list,
        refinment_iteration,
        symmetries,
        threshold,
        brillouin_primitive_vectors_3d,
        brillouin_primitive_vectors_2d,
        normalized_brillouin_primitive_vectors_2d
    )
    
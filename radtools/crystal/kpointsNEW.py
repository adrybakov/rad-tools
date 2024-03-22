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
import scipy
from scipy.interpolate import griddata
from scipy.spatial.transform import Rotation
from scipy.spatial import Delaunay
from radtools.geometry import absolute_to_relative
from scipy.spatial.distance import cdist
import math
from operator import itemgetter

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

###function checking if the k points in the list k list are inside the BZ
def downfold_inside_brillouin_zone(
    k_points_list,
    brillouin_primitive_vectors_3D):
   
    matrix_transformation_crystal_to_cartesian=np.zeros((3,3),dtype=float)
    matrix_transformation_cartesian_to_crystal=np.zeros((3,3),dtype=float)
    
    number_k_points = k_points_list.shape[0]
    transformed_k_points_list = np.zeros((number_k_points,3),dtype=float)
    
    matrix_transformation_crystal_to_cartesian=np.matrix(brillouin_primitive_vectors_3D).transpose()
    matrix_transformation_cartesian_to_crystal=np.linalg.inv(matrix_transformation_crystal_to_cartesian)

    ### writing k points in crystal coordinates
    for i in range(number_k_points):
        transformed_k_points_list[i,:]=matrix_transformation_cartesian_to_crystal@k_points_list[i,:]

    ### checking if any point of the list is outside the first brillouin zone
    ### to each k point is associated a int 0,1,2,3 if it goes over the boundaries
    indices=[]
    for i in range(number_k_points):
        count=0
        for r in range(3):
            if transformed_k_points_list[i,r]>=1.0 or transformed_k_points_list[i,r]<0:
                count+=r
                for d in range(3):
                    k_points_list[i,d]=k_points_list[i,d]-brillouin_primitive_vectors_3D[r,d]*np.sign(transformed_k_points_list[i,r])
        indices.append(count)
    return k_points_list,indices

### function applying symmetry analysis to a list of k points, the symmetry operations considered are the point group ones
### the fixed point is given; the analysis aim is to redistribute the fixed point weight between the list of k points
### the number of found symmetry orbits are counted in the redistribution of the fixed point weight
def symmetry_analysis(
    k_fixed_point,
    k_fixed_point_weight,
    k_points_list,
    symmetries,
    threshold_symmetry
):
    number_elements = k_points_list.shape[0]  

    new_k_points_list = np.c_[ k_points_list, (k_fixed_point_weight / number_elements) * np.ones(number_elements, dtype=float)]
    
    ### if there are no symmetry operations, the weight on each k point is exactly 1/(number_elements) of the original weight
    if ((symmetries is None) or (symmetries[0][0] == symmetries[0][1] == symmetries[0][2] == 0)):
        return new_k_points_list

    ### defining a matrix to save the degeneracies, after each symmetry operation
    check_degeneracy = np.zeros((number_elements, number_elements), dtype=bool)
    ### saving flag if no degeneracy is detected
    no_degeneracies = True
    
    eigenspaces=np.zeros(number_elements,dtype=int)
    for i in range(number_elements):
        eigenspaces[i]=i
    ### each k point is compared to each other k point after each of the symmetry operations has been applied
    for i in range(number_elements-1):
        for j in range(i+1,number_elements):
            ### for each pair considering all the symmetry operations
            for rotvec in symmetries:
                transformed_k_point = (
                    Rotation.from_rotvec(rotvec).as_matrix()
                    @ (k_points_list[j] - k_fixed_point)
                    + k_fixed_point
                )
                ### if any degeneracy is detected, the check on all the symmetry operations given is stopped
                if np.allclose(
                    k_points_list[i], transformed_k_point, atol=threshold_symmetry
                ):
                    check_degeneracy[i][j] = True
                    no_degeneracies = False
                    break

    ### if no degeneracy is detected, the same result as in the case of no symmetry operation is given
    if no_degeneracies:
        return new_k_points_list
    else:
        ### if degeneracy is detected, between the degenerate points only one is chosen as representative of the respective orbit (the first to appear in the k points list)
        ### at first assuming all k points as not-degenerate
        ### to each eigenspace (orbit) is associated a not null k_point
        for i in range(number_elements-1):
            for j in range(i+1,number_elements):
                if check_degeneracy[i][j]==True:
                    min_value=min(i,j)
                    max_value=max(i,j)
                    eigenspaces[eigenspaces==max_value]=min_value

        ### counting number of different eigenspaces
        number_eigenspaces=len(np.unique(eigenspaces))
        weights=k_fixed_point_weight/number_eigenspaces

        ### ordering eigenspaces values
        ordered_eigenspaces=np.zeros(number_elements)
        counting=0
        assigned_indices=[]
        for i in range(number_elements):
            if eigenspaces[i] not in assigned_indices:
                if counting<number_eigenspaces:
                    ordered_eigenspaces[eigenspaces==eigenspaces[i]]=counting
                    assigned_indices.extend([r for r in range(number_elements) if eigenspaces[r]==eigenspaces[i]])
                    counting+=1
        ### considering a representative for each eigenspace
        for i in range(number_eigenspaces):
            list_positions=[r for r in range(number_elements) if ordered_eigenspaces[r]==i]
            if len(list_positions)!=0:
                new_k_points_list[list_positions[0],3]=weights
                new_k_points_list[list_positions[1:],3]=0
            else:
                break

        ### returning the k points and their respective weights
        return new_k_points_list

def k_points_triangulation_2D(
    k_points_list,
    brillouin_primitive_vectors_2D,
    precedent_count
):
    
    number_elements=k_points_list.shape[0]
    k_points_list_projections=np.zeros((number_elements,2),dtype=float)
    ###choosing one point as an origin in the 2D brillouin plane
    origin=np.zeros(3,dtype=float)
    ###calculating the projections of the different k points on the primitive vectors
    vector_0=np.zeros((number_elements,2),dtype=float)
    vector_1=np.zeros((number_elements,2),dtype=float)
    vector_2=np.ones((number_elements,2),dtype=float)
    for i in range(number_elements):
        vector_0[i,0]=1
        vector_1[i,1]=1
        for j in range(2):
            k_points_list_projections[i,j]=(k_points_list[i,:3]-origin)@brillouin_primitive_vectors_2D[j,:]
    ###considering the periodic images as well
    k_points_list_projections_with_periodic_images=np.zeros((4*number_elements,2),dtype=float)
    k_points_list_projections_with_periodic_images[:number_elements,:]=k_points_list_projections
    k_points_list_projections_with_periodic_images[number_elements:2*number_elements,:]=k_points_list_projections+vector_0
    k_points_list_projections_with_periodic_images[2*number_elements:3*number_elements,:]=k_points_list_projections+vector_1
    k_points_list_projections_with_periodic_images[3*number_elements:,:]=k_points_list_projections+vector_2
    
    ###triangulation of the set of points
    triangles_tmp = Delaunay(k_points_list_projections_with_periodic_images,qhull_options="QJ")
    triangles_tmp=triangles_tmp.simplices
    ###each vertex has a pair of elements: the nuber of the k point associated, and if it the traingle is there passing through the BZ border
    number_triangles = int(len(triangles_tmp))
    triangles=[]
    for i in range(number_triangles):
        single_triangle=np.zeros((3,2),dtype=int)
        count_element=0
        count_outside=0
        for element in triangles_tmp[i]:
            if element<number_elements:
                single_triangle[count_element,1]=0
                single_triangle[count_element,0]=element
            elif element>=number_elements and element<2*number_elements:
                single_triangle[count_element,1]=1
                single_triangle[count_element,0]=element-number_elements
                count_outside+=1
            elif element>=2*number_elements and element<3*number_elements:
                single_triangle[count_element,1]=2
                single_triangle[count_element,0]=element-2*number_elements
                count_outside+=1
            else:
                single_triangle[count_element,1]=3
                single_triangle[count_element,0]=element-3*number_elements
                count_outside+=1
            count_element+=1
        if count_outside < 3:
            triangles.append(single_triangle)
    ###counting the proper number of triangles
    number_triangles = int(len(triangles))

    ###TESTING NEEDED HERE
    ###the projections along the two axis can be used to properly order the vertices of the triangls
    #ordered_triangles = []
    #vertices_projections = np.zeros((3,2),dtype=float)
    #one=[1,0]
    #two=[0,1]
    #three=[1,1]
    #for i in range(number_triangles):
    #    vertices=triangles[i]
    #    for r in range(3):
    #        if vertices[r,1]==0:
    #            vertices_projections[r]=k_points_list_projections[vertices[r,0]]
    #        elif vertices[r,1]==1:
    #            vertices_projections[r]=k_points_list_projections[vertices[r,0]]+one
    #        elif vertices[r,1]==2:
    #            vertices_projections[r]=k_points_list_projections[vertices[r,0]]+two
    #        else:
    #            vertices_projections[r]=k_points_list_projections[vertices[r,0]]+three
    #    indices=[tup[0] for tup in sorted(enumerate(vertices_projections[:,1]))]
    #    vertices_projections=vertices_projections[indices]
    #    if vertices_projections[0,1]==np.min(vertices_projections[:,0]):
    #        indices=[tup[0] for tup in sorted(enumerate(vertices_projections[:,0]))]
    #        vertices_projections=vertices_projections[indices]
    #    vertices=vertices[indices]
    #    ###considering the case the function is called on an existing list (the existing list has length equal to precedent_count)
    #    for s in range(3):
    #        if  vertices[s,0]> precedent_count:
    #            vertices[s,0]+=precedent_count
    #    ordered_triangles.append(vertices)
       
    ###for each triangle the ordering of the vertices is now clock-wise
    #return k_points_list,ordered_triangles
    return k_points_list,triangles

### local refinment of the k points list in the 2D plane
def local_refinment_with_symmetry_analysis(
    old_k_points_list,
    new_k_points_list,
    refinment_spacing,
    refinment_iterations,
    symmetries,
    threshold_symmetry,
    brillouin_primitive_vectors_3d,
    brillouin_primitive_vectors_2d,
    normalized_brillouin_primitive_vectors_2d,
    downfold_procedure=False
):
    length_old_k_points_list=old_k_points_list.shape[0]
    if  refinment_iterations == 0:
        if downfold_procedure == True:
            old_k_points_list[:,:3],indices= downfold_inside_brillouin_zone(
                old_k_points_list[:,:3],
                brillouin_primitive_vectors_3d
            )
        for i in range(length_old_k_points_list):
            if old_k_points_list[i,3]!=0.0:
                new_k_points_list.extend(old_k_points_list[i,:])
    else:
        iter = 0
        while iter != length_old_k_points_list:
            ### reading the k points inside the list and refine around them
            if length_old_k_points_list == 1:
                k_point_tmp = old_k_points_list[0,:3]
                k_point_tmp_weight = old_k_points_list[0,3]
                iter = length_old_k_points_list
            else:
                k_point_tmp = old_k_points_list[iter,:3]
                k_point_tmp_weight = old_k_points_list[iter,3]
                iter = iter + 1
            ### refinment procedure is applied to the single k point if its weight is not zero
            if k_point_tmp_weight != 0:
                ### a sub_grid of points is associated to each k point
                k_point_tmp_subgrid = np.zeros((4, 4),dtype=float)
                k_point_tmp_subgrid[:,:3] = k_point_tmp
                k_point_tmp_subgrid[0,:3] +=  normalized_brillouin_primitive_vectors_2d[0,:]*refinment_spacing
                k_point_tmp_subgrid[1,:3] -=  normalized_brillouin_primitive_vectors_2d[0,:]*refinment_spacing
                k_point_tmp_subgrid[2,:3] +=  normalized_brillouin_primitive_vectors_2d[1,:]*refinment_spacing
                k_point_tmp_subgrid[3,:3] -=  normalized_brillouin_primitive_vectors_2d[1,:]*refinment_spacing
                ### the weight of the k point is redistributed on the subgrid taking into account a symmetry analysis
                k_point_tmp_subgrid = symmetry_analysis(
                    k_point_tmp,
                    k_point_tmp_weight,
                    k_point_tmp_subgrid[:,:3],
                    symmetries,
                    threshold_symmetry
                )
                ### deleting the k points of the subgrid having weight equal to zero
                k_point_tmp_subgrid=np.delete(k_point_tmp_subgrid, np.where(k_point_tmp_subgrid[:,3]==0),axis=0)
                
                ### applying downfolding procedure
                if downfold_procedure == True:
                    k_point_tmp_subgrid[:,:3],indices=downfold_inside_brillouin_zone(
                        k_point_tmp_subgrid[:,:3],
                        brillouin_primitive_vectors_3d,
                    )

                ### applying reiteratively the refinment procedure
                new_refinment_spacing = refinment_spacing / 2

                local_refinment_with_symmetry_analysis(
                    k_point_tmp_subgrid,
                    new_k_points_list,
                    new_refinment_spacing,
                    refinment_iterations-1,
                    symmetries,
                    threshold_symmetry,
                    brillouin_primitive_vectors_3d,
                    brillouin_primitive_vectors_2d,
                    normalized_brillouin_primitive_vectors_2d,
                    downfold_procedure
                )

#generation of a 2D k points grid 
def k_points_generator_2D(
    new_brillouin_primitive_vectors_3d,
    old_brillouin_primitive_vectors_3d,
    chosen_plane,
    grid_spacing,
    covering_BZ=False,
    default_gridding=100,
    symmetry_analysis=False,
    refinment_spacing=0,
    refinment_iterations=0,
    symmetries=[[0,0,0]],
    threshold_symmetry=0.0001,
    shift_in_plane=[0,0],
    shift_in_space=[0,0,0],
    precedent_count=0
):
    ###individuating the 2D plane
    brillouin_primitive_vectors_2d=np.zeros((2, 3),dtype=float)
    count = 0
    for i in range(3):
        if chosen_plane[i]!=0:
            brillouin_primitive_vectors_2d[count] = new_brillouin_primitive_vectors_3d[i]
            count += 1

    ###initial 2D grid number of points (considering the grid_spacing input)
    if grid_spacing!=0.0:
        k0 = int(
            (brillouin_primitive_vectors_2d[0] @ brillouin_primitive_vectors_2d[0])/grid_spacing
            )
        k1 = int(
            (brillouin_primitive_vectors_2d[1] @ brillouin_primitive_vectors_2d[1])/grid_spacing
            )
        if k0==0 or k1==0:
            k0=default_gridding
            k1=default_gridding
    else:
        k0=default_gridding
        k1=default_gridding
    initial_weights = 1 / (k0 * k1)

    ###building the initial 2D k points grid
    k_points_grid=np.zeros((k0*k1,4),dtype=float)
    for i in range(k0):
        for j in range(k1):
            k_points_grid[i*k1+j,:3] =(float(i + shift_in_plane[0]) / k0) * (shift_in_space + brillouin_primitive_vectors_2d[0])\
                       + (float(j + shift_in_plane[1]) / k1) * (shift_in_space + brillouin_primitive_vectors_2d[1])
            k_points_grid[i*k1+j,3] = initial_weights
    
    ###symmetry analysis
    if symmetry_analysis == True:
        ###if not refinment spacing is given, the one from the building of the k points grid is used
        if refinment_spacing==0:
            refinment_spacing=grid_spacing/2.0

        ###normalized 2D plange generators, needed in the symmetry-refinment procedure
        normalized_brillouin_primitive_vectors_2d = np.zeros((2, 3),dtype=float)
        for i in range(2):
            normalized_brillouin_primitive_vectors_2d[i]=brillouin_primitive_vectors_2d[i]/(brillouin_primitive_vectors_2d[i]@brillouin_primitive_vectors_2d[i])
        
        ###applying the symmetry-refinment procedure to the 2D k points grid
        ###obtaining a refined list as an output
        ###loosing the ordering of the k points in the 2D plane
        refined_k_points_list = []
        local_refinment_with_symmetry_analysis(
            k_points_grid,
            refined_k_points_list,
            refinment_spacing,
            refinment_iterations,
            symmetries,
            threshold_symmetry,
            new_brillouin_primitive_vectors_3d,
            brillouin_primitive_vectors_2d,
            normalized_brillouin_primitive_vectors_2d,
            True
        )
        ###transforming the list into an array-like object (calling it list as well)
        refined_k_points_list = np.reshape(refined_k_points_list,(int(len(refined_k_points_list)/4), 4))
        ###if a covering of the BZ is not requested
        if covering_BZ == True:
            ###applying a triangulation procedure to order the k points list 
            ###the periodicity of the BZ is properly considered in the triangulation procedure
            refined_k_points_list,triangles=k_points_triangulation_2D(refined_k_points_list,brillouin_primitive_vectors_2d,precedent_count)
            return refined_k_points_list,triangles
        else:
            refined_k_points_list
    ###if no symmetry analysis is applied, building parallelograms to cover the BZ
    else:
        if covering_BZ == False:
            k_points_grid
        else:
            parallelograms=[]
            ###in case the downfolding procedure is requested with respect to a different BZ
            if np.array_equal(new_brillouin_primitive_vectors_3d,old_brillouin_primitive_vectors_3d):
                for i in range(0,k0):
                    for j in range(0,k1):
                        #periodicity of the BZ
                        parallelogram=[]
                        if i<k0-1 and j<k1-1:
                            parallelogram=[[i*k1+j,0],[(i+1)*k1+j,0],[(i+1)*k1+j+1,0],[i*k1+j+1,0]]
                        elif i<k0-1 and j==k1-1:
                            parallelogram=[[i*k1+j,0],[(i+1)*k1+j,0],[(i+1)*k1+0,2],[i*k1+0,2]]
                        elif i==k0-1 and j<k1-1:
                            parallelogram=[[i*k1+j,0],[0*k1+j,1],[0*k1+j+1,1],[i*k1+j+1,0]]
                        else:
                            parallelogram=[[i*k1+j,0],[0*k1+j,1],[0*k1+0,3],[i*k1+0,2]] 
                        parallelograms.append(parallelogram)
                parallelograms=np.reshape(parallelograms,(int(k0*k1),4,2))
            else:
                k_points_grid[:,:3],indices=downfold_inside_brillouin_zone(
                    k_points_grid[:,:3],
                    old_brillouin_primitive_vectors_3d)
                parallelogram=[]
                for i in range(k0):
                    for j in range(k1):
                        if i<k0-1 and j<k1-1:
                            parallelogram=[[i*k1+j,indices[i*k1+j]],[(i+1)*k1+j,indices[(i+1)*k1+j]],[(i+1)*k1+j+1,indices[(i+1)*k1+j+1]],[i*k1+j+1,indices[i*k1+j+1]]]
                            parallelograms.append(parallelogram)
                parallelograms=np.reshape(parallelograms,(len(parallelograms),4,2))        
            return k_points_grid,parallelograms

def check_inside_closed_shape_2d(
    eigenvectors_around_array,
    refined_k_points_list,
    brillouin_primitive_vectors_2d
):
    number_elements_border=eigenvectors_around_array.shape[0]
    number_elements=refined_k_points_list.shape[0]

    k_points_list_projections_border=np.zeros((number_elements_border,2),dtype=float)
    k_points_list_projections=np.zeros((number_elements_border,2),dtype=float)

    origin=np.zeros(3,dtype=float)
    ### calculating the projections of the different k points on the primitive vectors
    for i in range(number_elements_border):
        for j in range(2):
            k_points_list_projections_border[i,j]=(eigenvectors_around_array[i,:3]-origin)@brillouin_primitive_vectors_2d[j,:]
    for i in range(number_elements):
        for j in range(2):
            k_points_list_projections[i,j]=(refined_k_points_list[i,:3]-origin)@brillouin_primitive_vectors_2d[j,:]
    
    #point inclusion in polygon test W. Randolph Franklin (WRF) 
    def pnpoly(nvert,vert,test):
        i=0
        c=1
        while i < nvert:
            for j in range(i+1,nvert-1):
                if ( ((vert[i,1]>test[1]) != (vert[j,1]>test[1])) and (test[0]<(vert[j,0]-vert[i,0])*(test[1]-vert[i,1])/(vert[j,1]-vert[i,1])+vert[i,0]) ):    
                    c = -c
        return c
    
    indices=[]
    for i in range(number_elements):
        if pnpoly(number_elements_border,k_points_list_projections_border,k_points_list_projections[i])<0:
            indices.append(i)
    del refined_k_points_list[indices]

    return refined_k_points_list

def dynamical_refinment_little_paths_2D(
        little_paths,
        k_points_list,
        position_little_paths_to_refine,
        number_vertices,
        refinment_iterations,
        symmetries,
        threshold_symmetry,
        brillouin_primitive_vectors_3d,
        chosen_plane,
        brillouin_primitive_vectors_2d,
        normalized_brillouin_primitive_vectors_2d,
        threshold_minimal_refinment=0.00000000001,
        default_gridding=100
):
    new_k_points_list=k_points_list
    counting_offset=len(k_points_list)

    only_new_little_paths=[]
    number_little_paths=len(little_paths)
    
    eigenspaces_little_paths_to_refine=[[i] for i in position_little_paths_to_refine]

    flag_checking_indipendence=True

    while flag_checking_indipendence == True:
        number_eigenspaces=len(eigenspaces_little_paths_to_refine)
        number_total_eigenvectors=0
        number_total_little_paths_around=0
        little_paths_around_each_eigenspace=[]
        for eigenspace in eigenspaces_little_paths_to_refine:
            number_eigenvectors=len(eigenspace)
            number_total_eigenvectors+=number_eigenvectors
            ### considering the to-refine little paths in the selected eigenspace
            eigenvectors = little_paths[eigenspace]
            ### the little paths around each eigenvector are associated to it
            little_paths_around_each_eigenvector=[]
            for i in range(number_eigenvectors):
                little_paths_around_each_eigenvector.extend([r for r in range(number_little_paths) for j in little_paths[r][:,0] if j in eigenvectors[i][:,0] and r not in eigenspace])
            little_paths_around_each_eigenvector=list(np.unique(little_paths_around_each_eigenvector))
            number_total_little_paths_around+=len(little_paths_around_each_eigenvector)
            little_paths_around_each_eigenspace.append(little_paths_around_each_eigenvector)
        ### checking if two eigenspaces are in reality the same eigenspace
        ### between the "little_paths_around" of one of the two eigenspaces the index of one of the eigenvectors of the other eigenspace appears
        check_degeneracy = np.zeros((number_eigenspaces, number_eigenspaces), dtype=bool)
        any_degeneracy = 0
        for i in range(number_eigenspaces-1):
            for j in range(i+1,number_eigenspaces):
                for r in range(len(eigenspaces_little_paths_to_refine[i])):
                    if eigenspaces_little_paths_to_refine[i][r] in little_paths_around_each_eigenspace[j]:
                        check_degeneracy[i,j]=True
                        any_degeneracy+=1
                        break
                    else:
                        check_degeneracy[i,j]=False
        if any_degeneracy!=0:
            ## unifying the eigenspaces properly
            ## after having written them in an array form to facilitate the unification
            ## unifying as well the "little paths around"
            eigenspaces_little_paths_to_refine_array=np.zeros((number_total_eigenvectors,2),dtype=int)
            little_paths_around_each_eigenspace_array=np.zeros((number_total_little_paths_around,2),dtype=int)
            count1=0
            count2=0
            count3=0
            for eigenspace in eigenspaces_little_paths_to_refine:
                for eigenvector in eigenspace:
                    eigenspaces_little_paths_to_refine_array[count1,0]=eigenvector
                    eigenspaces_little_paths_to_refine_array[count1,1]=count2
                    count1+=1
                for around in little_paths_around_each_eigenspace[count2]:
                    little_paths_around_each_eigenspace_array[count3,0]=around
                    little_paths_around_each_eigenspace_array[count3,1]=count2
                    count3+=1
                count2+=1
            for i in range(number_eigenspaces-1):
                for j in range(i+1,number_eigenspaces):
                    if check_degeneracy[i][j]==True:
                        min_value=min(i,j)
                        max_value=max(i,j)
                        eigenvectors_positions=[r for r in range(number_total_eigenvectors) if eigenspaces_little_paths_to_refine_array[r,1]==max_value]
                        around_positions=[r for r in range(number_total_little_paths_around) if little_paths_around_each_eigenspace_array[r,1]==max_value]
                        eigenspaces_little_paths_to_refine_array[eigenvectors_positions,1]=min_value
                        little_paths_around_each_eigenspace_array[around_positions,1]=min_value
            ### redefining eigenspaces in order to take into account the found degeneracies
            values_eigenspaces=list(np.unique(eigenspaces_little_paths_to_refine_array[:,1]))
            eigenspaces_little_paths_to_refine=[]
            little_paths_around_each_eigenspace=[]
            for value in values_eigenspaces:
                eigenspaces_little_paths_to_refine.append([eigenspaces_little_paths_to_refine_array[r,0] for r in range(number_total_eigenvectors) if eigenspaces_little_paths_to_refine_array[r,1]==value])
                little_paths_around_each_eigenspace.append([little_paths_around_each_eigenspace_array[r,0] for r in range(number_total_little_paths_around) if little_paths_around_each_eigenspace_array[r,1]==value])
            ### eliminating in each eigenspace the presence of eigenvectors in other eigenvectors little paths around
            number_eigenspaces=len(eigenspaces_little_paths_to_refine)
            for i in range(number_eigenspaces):
                indices=[]
                count4=0
                for r in little_paths_around_each_eigenspace[i]:
                    for j in eigenspaces_little_paths_to_refine[i]:
                        if j == r:
                            indices.append(count4)
                    count4+=1
                indices.sort()
                if len(indices)!=0:
                    count4=0
                    for index in indices:
                        print(index)
                        del little_paths_around_each_eigenspace[i][index-count4]
                        count4+=1
            ### checking if two eigenspaces have one little path around in common 
            ### this in common little path is inserted into the to refine little paths and the procedure is repeated
            new_little_paths_to_refine=[]
            for i in range(number_eigenspaces-1):
                for j in range(i+1,number_eigenspaces):
                    for element in little_paths_around_each_eigenspace[i]:
                        if element in little_paths_around_each_eigenspace[j]:
                            new_little_paths_to_refine.append(element)
            if len(new_little_paths_to_refine)==0:
                flag_checking_indipendence=True
            else:
                new_little_paths_to_refine=list(np.unique(new_little_paths_to_refine))
                for s in new_little_paths_to_refine:
                    eigenspaces_little_paths_to_refine.append([s])
        else:
            flag_checking_indipendence=False
        print("cycling",eigenspaces_little_paths_to_refine,little_paths_around_each_eigenspace)
    print("little paths to refine:",eigenspaces_little_paths_to_refine)
    print("little paths around the ones to refine:",little_paths_around_each_eigenspace)
    ### substituting to the little paths positions, the k points positions in the list of k points
    number_eigenspaces=len(eigenspaces_little_paths_to_refine)
    eigenspaces_k_points_to_refine=[]
    eigenspaces_k_points_around=[]
    for eigenspace in eigenspaces_little_paths_to_refine:
        eigenspaces_k_points_to_refine.append(np.unique([vertex for eigenvector in eigenspace for vertex in little_paths[eigenvector]],axis=0))
    for eigenspace in little_paths_around_each_eigenspace:
        eigenspaces_k_points_around.append(np.unique([vertex for eigenvector in eigenspace for vertex in little_paths[eigenvector]],axis=0))    
    print("k points to refine:", np.shape(eigenspaces_k_points_to_refine), eigenspaces_k_points_to_refine)
    print("k points around the ones to refine:", np.shape(eigenspaces_k_points_around),eigenspaces_k_points_around)
    ### procede to a refinement of the k points constituing the to-refine little paths
    if number_vertices==3:
        k_points_list_tmp=[]
        # finding the minimum distance between the to-refine k points and the border k points
        for i in range(number_eigenspaces):
            number_eigenvectors=len(eigenspaces_k_points_to_refine[i])
            eigenvectors_array=np.zeros((number_eigenvectors,4),dtype=float)
            number_eigenvectors_around=len(eigenspaces_k_points_around[i])
            eigenvectors_around_array=np.zeros((number_eigenvectors_around,3),dtype=float)
            for j in range(number_eigenvectors):
                if eigenspaces_k_points_to_refine[i][j,1]==0:
                    eigenvectors_array[j]=k_points_list[eigenspaces_k_points_to_refine[i][j,0]]
                elif eigenspaces_k_points_to_refine[i][j,1]==1:
                    eigenvectors_array[j]=k_points_list[eigenspaces_k_points_to_refine[i][j,0]]+brillouin_primitive_vectors_2d[0]
                else:
                    eigenvectors_array[j]=k_points_list[eigenspaces_k_points_to_refine[i][j,0]]+brillouin_primitive_vectors_2d[1]
            
            minimal_distance = (cdist(eigenvectors_array-eigenvectors_around_array)).min()
            refined_k_points_list=[]
            
            local_refinment_with_symmetry_analysis(
                    eigenvectors_array,
                    refined_k_points_list,
                    minimal_distance,
                    refinment_iterations,
                    symmetries,
                    threshold_symmetry,
                    brillouin_primitive_vectors_3d,
                    brillouin_primitive_vectors_2d,
                    normalized_brillouin_primitive_vectors_2d,
                    False
                    )

            ### check if the points are inside the border 
            refined_k_points_list[:,:3] = check_inside_closed_shape_2d(
                                        eigenvectors_around_array,
                                        refined_k_points_list[:,:3],
                                        brillouin_primitive_vectors_2d)
            ### triangulation
            refined_k_points_list,triangles=k_points_triangulation_2D(
                                                    refined_k_points_list,
                                                    brillouin_primitive_vectors_2d,
                                                    counting_offset)
            counting_offset+=len(refined_k_points_list)
            k_points_list_tmp.append(refined_k_points_list)
            only_new_little_paths.append(triangles)
        new_k_points_list.append(k_points_list_tmp)
    else:
        ### here it is sufficient to find for each eigenspace the extremal k points
        ### these k points are then used to draw a simil-BZ, which is refined using the grid-generation 2D methods
        ### first up-folding the k points 
        for i in range(number_eigenspaces):
            number_eigenvectors_around=len(eigenspaces_k_points_around[i])
            eigenvectors_around_array=np.zeros((number_eigenvectors_around,4),dtype=float)
            for j in range(number_eigenvectors_around):
                if eigenspaces_k_points_around[i][j][1]==3:
                    eigenvectors_around_array[j,:3]=k_points_list[eigenspaces_k_points_around[i][j][0],:3]+brillouin_primitive_vectors_2d[1]+brillouin_primitive_vectors_2d[0]
                elif eigenspaces_k_points_around[i][j][1]==2:
                    eigenvectors_around_array[j,:3]=k_points_list[eigenspaces_k_points_around[i][j][0],:3]+brillouin_primitive_vectors_2d[1]
                elif eigenspaces_k_points_around[i][j][1]==1:
                    eigenvectors_around_array[j,:3]=k_points_list[eigenspaces_k_points_around[i][j][0],:3]+brillouin_primitive_vectors_2d[0]
                else:
                    eigenvectors_around_array[j]=k_points_list[eigenspaces_k_points_around[i][j][0]]
                eigenvectors_around_array[j,3]=k_points_list[eigenspaces_k_points_around[i][j][0],3]
            ### calculating the projections of the different k points on the primitive vectors
            k_points_list_projections_border=np.zeros((number_eigenvectors_around,2),dtype=float)
            k_points_list_projections_border_sum=np.zeros((number_eigenvectors_around),dtype=float)
            for i in range(number_eigenvectors_around):
                for j in range(2):
                    k_points_list_projections_border[i,j]=eigenvectors_around_array[i,:3]@brillouin_primitive_vectors_2d[j,:]
                k_points_list_projections_border_sum[i]=np.sum(k_points_list_projections_border[i])
            ### sorting in order to find the simil-BZ
            origin_ij=np.argmin(k_points_list_projections_border_sum)
            origin=eigenvectors_around_array[origin_ij,:3]
            for i in range(number_eigenvectors_around):
                for j in range(2):
                    k_points_list_projections_border[i,j]=(eigenvectors_around_array[i,:3]-origin)@brillouin_primitive_vectors_2d[j,:]
            ### finding the bottom-right point
            min_0=k_points_list_projections_border[0,0]
            max_1=k_points_list_projections_border[0,1]
            position_pair_1=0
            for i in range(number_eigenvectors_around):
                if k_points_list_projections_border[i,0]<min_0 and  k_points_list_projections_border[i,1]>max_1:
                    min_0=k_points_list_projections_border[i,0]
                    max_1=k_points_list_projections_border[i,1]
                    position_pair_1=i
            ### finding the top-left point
            max_0=k_points_list_projections_border[0,0]
            min_1=k_points_list_projections_border[0,1]
            position_pair_0=0
            for i in range(number_eigenvectors_around):
                if k_points_list_projections_border[i,0]>max_0 and  k_points_list_projections_border[i,1]<min_1:
                    max_0=k_points_list_projections_border[i,0]
                    min_1=k_points_list_projections_border[i,1]
                    position_pair_0=i
            ### defining the new bordening
            new_brillouin_primitive_vectors_2d=np.zeros((2,3),dtype=float)
            new_brillouin_primitive_vectors_2d[0]=eigenvectors_around_array[position_pair_1,:3]-origin
            new_brillouin_primitive_vectors_2d[1]=eigenvectors_around_array[position_pair_0,:3]-origin
            new_brillouin_primitive_vectors_3d=np.zeros((3,3),dtype=float)
            count = 0
            for i in range(3):
                if chosen_plane[i] != 0:
                    new_brillouin_primitive_vectors_3d[i] = new_brillouin_primitive_vectors_2d[count]
                    count += 1
                else:
                    new_brillouin_primitive_vectors_3d[i] = brillouin_primitive_vectors_3d[i]
            ### finding the minimal grid spacing
            max_value=np.linalg.norm(new_brillouin_primitive_vectors_2d[0]+new_brillouin_primitive_vectors_2d[1],2)
            moduli=np.zeros(number_eigenvectors_around,dtype=float)
            for i in range(number_eigenvectors_around):
                moduli[i]=np.linalg.norm(eigenvectors_around_array[i,:3]-origin,2)
                if moduli[i] < threshold_minimal_refinment:
                    moduli[i]=max_value
            initial_grid_spacing=np.min(moduli)/2
            new_k_points_list_tmp,parallelograms=k_points_generator_2D(
                                                    new_brillouin_primitive_vectors_3d,
                                                    brillouin_primitive_vectors_3d,
                                                    chosen_plane,
                                                    initial_grid_spacing,
                                                    True,
                                                    default_gridding,
                                                    False,
                                                    0,
                                                    0,
                                                    [[0,0,0]],
                                                    0.0,
                                                    [0,0],
                                                    origin,
                                                    counting_offset
                                                )                        
            for parallelogram in parallelograms:
                for vertex in parallelogram:
                    vertex[1]+=counting_offset
            counting_offset+=len(new_k_points_list_tmp)

            new_k_points_list=np.append(new_k_points_list,new_k_points_list_tmp,axis=0)
            only_new_little_paths.extend(parallelograms)

    return new_k_points_list,only_new_little_paths

def printing_covering_BZ_2D(
    brillouin_primitive_vectors_3d,
    chosen_plane,
    k_points_list,
    little_paths,
    number_vertices
    ):
    
    brillouin_primitive_vectors_2d=np.zeros((2,3),dtype=float)
    count = 0
    for i in range(3):
        if chosen_plane[i]!=0:
            brillouin_primitive_vectors_2d[count] = brillouin_primitive_vectors_3d[i]
            count += 1
    
    number_little_paths=len(little_paths)
    little_path=[0]*number_vertices
    ### substituing to the vertices the k points coordinates
    new_little_paths=[]
    for i in range(number_little_paths):
        indx_1=[little_paths[i][r][0] for r in range(number_vertices)]
        indx_2=[little_paths[i][r][1] for r in range(number_vertices)]
        count=0
        for indx in indx_1:
            little_path[count]=[k_points_list[indx,s] for s in range(3)]
            if indx_2[count]!=0:
                for s in range(3):
                    if indx_2[count]==1:
                        little_path[count][s]+=brillouin_primitive_vectors_2d[0,s]
                    elif indx_2[count]==2:
                        little_path[count][s]+=brillouin_primitive_vectors_2d[1,s]
                    else:
                        little_path[count][s]+=brillouin_primitive_vectors_2d[1,s]+brillouin_primitive_vectors_2d[0,s]
            count+=1
        new_little_paths.extend(little_path)

    ### reordering the list
    number_little_paths=int(len(new_little_paths)/number_vertices)
    new_little_paths=np.reshape(new_little_paths,(number_little_paths,number_vertices,3))
    ### considering only coordinates in the chosen plane
    new_little_paths_projected=np.zeros((number_little_paths,number_vertices,2),dtype=float)
    for i in range(number_little_paths):
        for j in range(number_vertices):
            for r in range(2):
                new_little_paths_projected[i,j,r]+=(new_little_paths[i,j]@brillouin_primitive_vectors_2d[r])
    fig,ax=plt.subplots()
    patches=[]
    for i in range(number_little_paths):
        polygon=Polygon(new_little_paths_projected[i],closed=True,fill=None,edgecolor='c')
        patches.append(polygon)
    p = PatchCollection(patches,cmap=matplotlib.cm.jet,alpha=0.4)
    colors = 100*np.random.rand(len(patches))
    p.set_array(np.array(colors))
    ax.add_collection(p)
    fig=plt.figure()
    ax=fig.add_subplot(projection='3d')
    ax.scatter(k_points_list[:,0],k_points_list[:,1],k_points_list[:,3])
    plt.show()

### TESTING INPUT
if __name__ == "__main__":
    import timeit
    import os
    import matplotlib
    from matplotlib.artist import Artist
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
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
    brillouin_primitive_vectors_3d[0]=[1,0,0]
    brillouin_primitive_vectors_3d[1]=[0,1,0]
    brillouin_primitive_vectors_3d[2]=[0,0,1]
    chosen_plane=[1,1,0]
    symmetries=[[0,0,np.pi/2]]
    grid_spacing=0.1
    refinment_iterations=3
    refinment_spacing=0.1
    threshold_symmetry=0.001
    epsilon=0.000000001
    default_gridding=1000
    count=0
    shift_in_plane=[0,0]
    shift_in_space=[0,0,0]
    covering_BZ=True

    ### LITTLE PARALLELOGRAMS
    ##number_vertices=4
    ##not_refined_k_points_list,parallelograms=k_points_generator_2D(
    ##    brillouin_primitive_vectors_3d,
    ##    brillouin_primitive_vectors_3d,
    ##    chosen_plane,
    ##    grid_spacing,
    ##    covering_BZ,
    ##    default_gridding,
    ##    False,
    ##    refinment_spacing,
    ##    refinment_iterations,
    ##    symmetries,
    ##    threshold_symmetry,
    ##    shift_in_plane,
    ##    shift_in_space,
    ##    0
    ##)
    ##printing_covering_BZ_2D(
    ##    brillouin_primitive_vectors_3d,
    ##    chosen_plane,
    ##    not_refined_k_points_list,
    ##    parallelograms,
    ##    number_vertices
    ##)

    ###LITTLE TRIANGLES
    number_vertices=3
    refined_k_points_list,triangles=k_points_generator_2D(
        brillouin_primitive_vectors_3d,
        brillouin_primitive_vectors_3d,
        chosen_plane,
        grid_spacing,
        covering_BZ,
        default_gridding,
        True,
        refinment_spacing,
        refinment_iterations,
        symmetries,
        threshold_symmetry,
        shift_in_plane,
        shift_in_space,
        0
        )
    printing_covering_BZ_2D(
        brillouin_primitive_vectors_3d,
        chosen_plane,
        refined_k_points_list,
        triangles,
        number_vertices
    )
    
    
    

  
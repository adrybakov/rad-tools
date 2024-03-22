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

### function applying symmetry analysis to a list of k points, the symmetry operations considered are the point group ones
### the fixed point is given; the analysis aim is to redistribute the fixed point weight between the list of k points
### the found symmetry orbits are counted in the redistribution of the fixed point weight a number of times equal to the number of k points inside the orbit itself
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
        for i in range(number_eigenspaces):
            list_positions=[r for r in range(number_elements) if ordered_k_eigenspaces[r]==i]
            if len(list_positions)!=0:
                new_k_list[list_positions[0],3]=weights
                new_k_list[list_positions[1:],3]=0
            else:
                break

        # returning the matrix with the respective weights
        return new_k_list






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
        for i in range(length_old_k_points_list):
            if old_k_points_list[i,4]!=0.0:
                new_k_points_list.extend(old_k_points_list[i,:])
        if downfold_procedure == True:
            new_k_points_list = np.reshape(new_k_points_list,(int(len(new_k_points_list)/4), 4))
            new_k_points_list[:,:3]= downfold_inside_brillouin_zone(
                new_k_points_list[:,:3],
                brillouin_primitive_vectors_3d
            )
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
                    k_point_tmp_subgrid[:,:3]=downfold_inside_brillouin_zone(
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
    shift_in_space=[0,0,0]
):
    ###individuating the 2D plane
    brillouin_primitive_vectors_2d=np.zeros((2, 3),dtype=float)
    count = 0
    for i in chosen_plane:
        if i!=0:
            brillouin_primitive_vectors_2d[count] = new_brillouin_primitive_vectors_3d[i]
            count += 1

    ###initial 2D grid number of points (considering the grid_spacing input)
    if grid_spacing!=0.0:
        k0 = int(
            (brillouin_primitive_vectors_2d[0] @ brillouin_primitive_vectors_2d[0])/grid_spacing
            )
        k1 = int(
            (brillouin_primitive_vectors_2d[1] @ brillouin_primitive_vectors_2d[1]) /grid_spacing
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
            k_points_grid[i*k1+j,:3] = (float(i + shift_in_plane[0]) / k0) * (shift_in_space + brillouin_primitive_vectors_2d[0,:]) + (float(j + shift_in_plane[1]) / k1) * (shift_in_space + brillouin_primitive_vectors_2d[1,:])
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
        ###the list should be already ordered
        ###transforming the list into an array-like object (calling it list as well)
        ###refined_k_points_list = np.reshape(refined_k_points_list,(int(len(refined_k_points_list)/4), 4))
        ###if a covering of the BZ is not requested
        if covering_BZ == True:
            ###applying a triangulation procedure to order the k points list 
            ###the periodicity of the BZ is properly considered in the triangulation procedure
            refined_k_points_list,triangles=k_points_triangulation_2D(refined_k_points_list,brillouin_primitive_vectors_2d,0)
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
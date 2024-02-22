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

from copy import deepcopy
from math import sqrt
import numpy as np
import os 
import ipyparallel as ipp
from radtools import MagnonDispersion
from scipy.spatial.transform import Rotation
from radtools.crystal.kpoints import dynamical_refinment_triangulation_2d
from radtools.crystal.kpoints import k_points_generator_2D
from radtools.geometry import span_orthonormal_set
from radtools.magnons.diagonalization import ColpaFailed, solve_via_colpa
from radtools.spinham.hamiltonian import SpinHamiltonian
from radtools.magnons import dispersion


def update_magnonic_branches_degeneracies(
        degeneracy_matrix, number_magnons, magnonic_branches, number_eigenspaces
    ):
    r"""
        Considering the degeneracy matrix, the knowledge in the magnonic branches subdivision is updated
        
        Parameters
        ---------
        number_magnons: n |int| number of different magnonic branches (number of magnetic atoms in the unit cell) 
        degeneracy_matrix: (n,n) |array_like| (Boolean Values) degeneracies detected between the different magnonic branches
        old_magnonic_branches: (n) |array_like| to each element is associated a different magnonic branch, the value of each element is equal to the respective eigenspace
        (in case of a degenerate eigenspace two or more elements are expected to be equal
        i.e. considering the magnonic branches 1 and 5, if they are degenerate we expect to have the elements 0 and 4 equal to the same number, which is the respective eigenspace)
        number_eigenspaces: |int| number of eigenspaces before updating the knwoledge about degeneracies

        Return
        ------
        magnonic_branches: (n) |array_like|
        number_eigenspace: |int| number of eigenspace after updating the knowledge about degeneracies        
    """ 
    # if there is already only one eigenspace, then no study of degeneracies is really needed
    if number_eigenspaces == 1:
        return magnonic_branches, number_eigenspaces
    
    # studying the degeneracy matrix and updating the respective knowledge
    number_eigenspaces_tmp=number_eigenspaces
    for i in range(0, number_magnons - 1):
        for j in range(i + 1, number_magnons):
            if number_eigenspaces_tmp == 1:
                break
            else:
                if degeneracy_matrix[i][j]==True:
                    if magnonic_branches[i]==magnonic_branches[j]:
                        break
                    else:
                        min_eigenspace=min(magnonic_branches[i],magnonic_branches[j])
                        max_eigenspace=max(magnonic_branches[i],magnonic_branches[j])
                        magnonic_branches[magnonic_branches==max_eigenspace]=min_eigenspace
                        number_eigenspaces_tmp==number_eigenspaces_tmp-1
    
    # the eigenspaces have to be properly numerated in case there is any change
    ordered_magnonic_branches=np.zeros(number_magnons)
    counting=0
    for i in range(number_magnons):
        if counting < number_eigenspaces:
            ordered_magnonic_branches[magnonic_branches==magnonic_branches[i]]=counting
            counting+=1

    return magnonic_branches, number_eigenspaces

class Berry_curvature:
    r"""
       Berry Curvature in a 2D plane
       the Berry curvature is calculated locally as the line integral of the Berry connection over little parallelograms or triangles
       this formulation preserve the Berry curvature gauge invariance

       Parameters
       ----------
       the ones needed in the production of the 2D k points list (through the refinment procedure) plus the ones needed in the Berry curvature calculation 
    """
    def __init__(
    self,
    spinham,
    nodmi,
    noaniso,
    ) -> None:
        self.spinham=spinham
        self.nodmi=nodmi
        self.noaniso=noaniso
        self.dispersion=MagnonDispersion(self.spinham,nodmi=self.nodmi,noaniso=self.noaniso)
        self.N=self.dispersion.N

    def magnonic_surfaces(
        self, k_points_list, noeigenvectors, nocheckdegeneracy, threshold_omega, grid
    ):
        if grid == True:
            np.reshape(k_points_list,k_points_list.shape[0]*k_points_list.shape[1],order='C')
        number_k_points=len(k_points_list)
        omegas=np.zeros((number_k_points,self.N),dtype=complex)
        us=np.zeros((number_k_points,self.N,self.N),dtype=complex)
        number_eigenspaces = self.N
        magnonic_branches = np.asarray([i for i in range(self.N)])
        
        if (nocheckdegeneracy == True) and (noeigenvectors == True):
            for i in range(number_k_points):
                omegas[i]=self.dispersion(k_points_list[i][:3],False)
                omegas[i]=omegas[i]*k_points_list[i][3]
        elif (nocheckdegeneracy == True) and (noeigenvectors == False):
            us=np.zeros((number_k_points,self.N,self.N),dtype=complex)
            for i in range(number_k_points):
                omegas[i],us[i]=self.dispersion(k_points_list[i][:3],True)
                omegas[i]=omegas[i]*k_points_list[i][3]
                us[i]=us[i]*k_points_list[i][3]
        elif (nocheckdegeneracy == False) and (noeigenvectors == True):
            degeneracy_matrix = np.zeros((self.N,self.N),dtype=bool)  
            for i in range(number_k_points):
                omegas[i]=self.dispersion(k_points_list[i][:3],False)
                omegas[i]=omegas[i]*k_points_list[i][3]
                degeneracy=0
                for r in range(0,self.N-1):
                    for s in range(r+1,self.N):
                        degeneracy_matrix[i] = False
                        degeneracy_matrix[i] = np.isclose(omegas[i][r],omegas[i][s],ato=threshold_omega)
                        if degeneracy_matrix[i] == True:
                            degeneracy+=1
                if degeneracy != 0:
                    magnonic_branches,number_eigenspaces=update_magnonic_branches_degeneracies(
                        degeneracy_matrix, self.N, magnonic_branches, number_eigenspaces)     
        else:
            degeneracy_matrix = np.zeros((self.N,self.N),dtype=bool)
            for i in range(number_k_points):
                omegas[i],us[i]=self.dispersion(k_points_list[i][:3],True)
                omegas[i]=omegas[i]*k_points_list[i][3]
                us[i]=us[i]*k_points_list[i][3]
                degeneracy=0
                for r in range(0,self.N-1):
                    for s in range(r+1,self.N):
                        degeneracy_matrix[i] = False
                        degeneracy_matrix[i] = np.isclose(omegas[i][r],omegas[i][s],ato=threshold_omega)
                        if degeneracy_matrix[i] == True:
                            degeneracy+=1
                if degeneracy != 0:
                    magnonic_branches,number_eigenspaces=update_magnonic_branches_degeneracies(
                        degeneracy_matrix, self.N, magnonic_branches, number_eigenspaces)     
        
        return omegas,us,magnonic_branches,number_eigenspaces

    def circuitation_over_a_path(
        u_values_along_the_path,k_point_weights_along_the_path,number_modes,number_eigenspaces,magnonic_branches,number_vertices
    ):
        phases=np.zeros(number_eigenspaces,dtype=complex)
        if number_eigenspaces == number_modes:
            count=0
            while count<number_vertices-1:
                phases=phases+np.angle(np.asarray(u_values_along_the_path[count]@u_values_along_the_path[count+1]))*k_point_weights_along_the_path[count]
                count+=1
            phases=phases+np.angle(np.asarray(u_values_along_the_path[number_vertices-1]@u_values_along_the_path[0]))*k_point_weights_along_the_path[number_vertices-1]  
        ##In the degenerate case the Berry curvature is calculated using the Non-Abelian formulation
        else:
            for n in range(number_eigenspaces):
                bands=[magnonic_branches==n]
                count=0
                while count<number_vertices-1:
                    phases[n]=phases[n]+np.angle(np.asarray(u_values_along_the_path[count][np.ix_(bands,bands)] @ u_values_along_the_path[count+1][np.ix_(bands,bands)]))*k_point_weights_along_the_path[count]
                    count+=1
                phases[n]=phases[n]+np.angle(np.asarray(u_values_along_the_path[number_vertices-1][np.ix_(bands,bands)] @ u_values_along_the_path[0][np.ix_(bands,bands)]))*k_point_weights_along_the_path[number_vertices-1]
        return phases

    def berry_curvature(
        self,
        triangulation,
        brillouin_primitive_vectors,
        chosen_plane,
        initial_grid_spacing,
        shift_in_plane,
        shift_in_space,
        symmetries,
        threshold_k_point,
        refinment,
        refinment_spacing,
        refinment_iteration,
        threshold_omega,
        dynamical_refinment,
        dynamical_refinment_iteration,
        threshold_dynamical_refinment
    ):  
        refined_k_points_list,little_paths=k_points_generator_2D(
                                brillouin_primitive_vectors,
                                chosen_plane,
                                initial_grid_spacing,
                                shift_in_plane,
                                shift_in_space,
                                symmetries,
                                threshold_omega,
                                refinment,
                                refinment_spacing,
                                refinment_iteration,
                                triangulation
                            )
        _,us,magnonic_branches,number_eigenspaces=self.magnonic_surfaces(refined_k_points_list,False,False,threshold_omega,False)

        if triangulation==True:
            number_vertices=3
        else:
            number_vertices=4

        u=np.zeros((number_vertices,self.N,self.N),dtype=complex)
        weights=np.zeros(number_vertices,dtype=float)
        old_chern=np.zeros(number_eigenspaces,dtype=float)
        chern=np.zeros(number_eigenspaces,dtype=float)

        while dynamical_refinment==True:
            number_little_paths=little_paths.shape[0]
            _,us,magnonic_branches,number_eigenspaces=self.magnonic_surfaces(refined_k_points_list,False,False,threshold_omega,False)
            i=0
            phases=np.zeros((number_little_paths,number_eigenspaces),dtype=complex)
            while i < number_little_paths:
                for j in range(number_vertices):
                    u[j]=us[little_paths[i][j]]
                    weights[j]=refined_k_points_list[little_paths[i][j],3]
                phases[i]=self.circuitation_over_a_path(u,weights,self.N,number_eigenspaces,magnonic_branches,triangulation)
                chern+=phases
            # condition for dynamical refinment
            if np.any(np.asarray(list(map(lambda x: float(x).is_integer(),chern)))) == False:
                old_chern+=chern
                max_indices=[]
                max_phases=np.max(phases,axis=0)
                max_index=np.argmax(max_phases)
                for j in range(number_little_paths):
                    if j != max_index:
                        if (max_phases[j]<max_phases[max_index]+threshold_dynamical_refinment) and (max_phases[j]>max_phases[max_index]-threshold_dynamical_refinment):
                            max_indices.append(j)
                selected_little_paths=little_paths[max_indices,:]
                refined_little_paths,little_paths=dynamical_refinment_2d(
                                                        little_paths,
                                                        selected_little_paths,
                                                        refined_k_points_list,
                                                        dynamical_refinment,
                                                        dynamical_refinment_iteration,
                                                        symmetries,
                                                        threshold_omega,
                                                        brillouin_primitive_vectors,
                                                        chosen_plane,
                                                        triangulation
                                                    )
            else:
                dynamical_refinment=False

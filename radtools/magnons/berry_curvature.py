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
                omegas[i]=omegas[i,j]*k_points_list[i][3]
                us[i,j]=us[i,j]*k_points_list[i][3]
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
        u_values_along_the_path,k_point_weights_along_the_path,number_modes,number_eigenspaces,magnonic_branches,triangulation
    ):
        if triangulation == True:
            if number_eigenspaces == number_modes:
                phase = np.asarray(list(map(
                    lambda x,xw,y,yw,z,zw: 
                    np.angle(np.asarray(x))*xw
                    + np.angle(np.asarray(y))*yw
                    + np.angle(np.asarray(z))*zw,
                    (u_values_along_the_path[0] @ u_values_along_the_path[1]),k_point_weights_along_the_path[0],
                    (u_values_along_the_path[1] @ u_values_along_the_path[2]),k_point_weights_along_the_path[1],
                    (u_values_along_the_path[2] @ u_values_along_the_path[0]),k_point_weights_along_the_path[2])))

             ##in the degenerate case the Berry curvature is calculated using the non-abelian formulation
            else:
                for n in range(number_eigenspaces):
                    bands=[magnonic_branches==n]
                    phase[n] = np.asarray(list(map(
                    lambda x,xw,y,yw,z,zw: 
                    np.angle(np.asarray(x))*xw
                    + np.angle(np.asarray(y))*yw
                    + np.angle(np.asarray(z))*zw,
                    (u_values_along_the_path[0][np.ix_(bands,bands)] @ u_values_along_the_path[1][np.ix_(bands,bands)]),k_point_weights_along_the_path[0],
                    (u_values_along_the_path[1][np.ix_(bands,bands)] @ u_values_along_the_path[2][np.ix_(bands,bands)]),k_point_weights_along_the_path[1],
                    (u_values_along_the_path[2][np.ix_(bands,bands)] @ u_values_along_the_path[0][np.ix_(bands,bands)]),k_point_weights_along_the_path[2])))
        else:
            if number_eigenspaces == number_modes:
                phase = np.asarray(list(map(
                    lambda x,xw,y,yw,z,zw,w,ww: 
                    np.angle(np.asarray(x))*xw
                    + np.angle(np.asarray(y))*yw
                    - np.angle(np.asarray(z))*zw
                    - np.angle(np.asarray(w))*ww,
                    (u_values_along_the_path[0] @ u_values_along_the_path[1]),k_point_weights_along_the_path[0],
                    (u_values_along_the_path[1] @ u_values_along_the_path[2]),k_point_weights_along_the_path[1],
                    (u_values_along_the_path[3] @ u_values_along_the_path[2]),k_point_weights_along_the_path[3],
                    (u_values_along_the_path[0] @ u_values_along_the_path[3]),k_point_weights_along_the_path[0])))
             ##in the degenerate case the Berry curvature is calculated using the non-abelian formulation
            else:
                for n in range(number_eigenspaces):
                    bands=[magnonic_branches==n]
                    phase[n] = np.asarray(list(map(
                    lambda x,xw,y,yw,z,zw, w,ww: 
                    np.angle(np.asarray(x))*xw
                    + np.angle(np.asarray(y))*yw
                    - np.angle(np.asarray(z))*zw
                    - np.angle(np.asarray(w))*ww,
                    (u_values_along_the_path[0][np.ix_(bands,bands)] @ u_values_along_the_path[1][np.ix_(bands,bands)]),k_point_weights_along_the_path[0],
                    (u_values_along_the_path[1][np.ix_(bands,bands)] @ u_values_along_the_path[2][np.ix_(bands,bands)]),k_point_weights_along_the_path[1],
                    (u_values_along_the_path[3][np.ix_(bands,bands)] @ u_values_along_the_path[2][np.ix_(bands,bands)]),k_point_weights_along_the_path[3],
                    (u_values_along_the_path[0][np.ix_(bands,bands)] @ u_values_along_the_path[3][np.ix_(bands,bands)]),k_point_weights_along_the_path[0])))
        return phase

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
        if triangulation == True:
            refined_k_points_list,triangles=k_points_generator_2D(
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
                                    False,
                                    True
                                )
            
            _,us,magnonic_branches,number_eigenspaces=self.magnonic_surfaces(k_points_list,False,False,threshold_omega,False)
            
            u=np.zeros((3,self.N,self.N),dtype=complex)
            weight_a=np.zeros(3,dtype=float)
            phase=np.zeros(number_eigenspaces,dtype=complex)
            chern_number=np.zeros(number_eigenspaces)
           
            position_elements_to_refine_list=[]
            ##taking into account possible degeneracies is equivalent to adopt the non-Abelian formulation of the Berry Curvature
            number_triangles=triangles.shape[0]

            while i < number_triangles:
                position_elements_to_refine=[]
                for j in range(3):
                    u[j]=us[triangles[i][j]]
                    weight_a=k_points_list[triangles[i][j],3]
                    position_elements_to_refine.append(triangles[i][j])
                phase=self.circuitation_over_a_path(u,weight_a,self.N,number_eigenspaces,magnonic_branches,triangulation)
                
                #saving all divergent points: if the refinment procedure is applied, the grid will be refined around these points
                if np.max(phase)>threshold_dynamical_refinment:
                    position_elements_to_refine_list.extend(position_elements_to_refine)
                else:
                    chern_number+=phase 
                
                ###dynamical refinment
                ###the degeneracy is not checked again
                ###the circuitation is calculated only on the new triangles
                if dynamical_refinment==True:
                    new_bias=len(refined_k_points_list)
                    refined_k_points_list,triangles=dynamical_refinment_triangulation_2d(
                                                        triangles,
                                                        position_elements_to_refine_list,
                                                        refined_k_points_list,
                                                        dynamical_refinment_iteration,
                                                        symmetries,
                                                        threshold_omega,
                                                        brillouin_primitive_vectors,
                                                        chosen_plane
                                                    )
                    new_triangles=triangles[new_bias:,:]
                    new_k_points=refined_k_points_list[new_bias:,:]
                    _,us,magnonic_branches,number_eigenspaces=self.magnonic_surfaces(new_k_points,False,True,threshold_omega)
                    number_new_triangles=new_triangles.shape[0]
                    
                    while i < number_new_triangles:
                        for j in range(3):
                            u[j]=us[new_triangles[i][j]]
                            weight_a=k_points_list[new_triangles[i][j],3]
                        phase=self.circuitation_over_a_path(u,weight_a,self.N,number_eigenspaces,magnonic_branches,triangulation)
                        chern_number+=phase 
                    

        AGGIUNGERE CIRCUITAZIONE LUNGO REETTANGOLI


                NOT NEEDED TO RECALCULATE EVERYTHING AGAIN.....
                PROPER INTEGRATION OF THE DYNAMICAL ROUTINE....



            u=np.zeros((4,self.N,self.N),dtype=complex)
            weight_a=np.zeros(4,dtype=float)
            # the dynamical refinment as formulated requires to interpolate all the k grid again
            # and to recalculate all the magnonic surfaces
            cycling_through=True
            while cycling_through==True:
                position_elements_to_refine_list=[]
                ##taking into account possible degeneracies is equivalent to adopt the non-Abelian formulation of the Berry Curvature
                phase=np.zeros(number_eigenspaces,dtype=complex)
                chern_number=np.zeros(number_eigenspaces)
                i=0
                j=0
                while i < n0:
                    while j < n1:
                        print(
                            "berry grid",
                            int(i / n0 * 100),
                            "%",
                            end="\r",
                            flush=True,
                        )
                    ## considering the little parallepiped path
                    if i < n0 -1 and j < n1 -1:
                        u[0]=us[i][j]
                        u[3]=us[i][j+1]
                        u[2]=us[i+1][j+1]
                        u[1]=us[i+1][j]
                        weight_a[0]=k_points_grid[i,j][3]+k_points_grid[i+1,j][3]
                        weight_a[1]=k_points_grid[i+1,j][3]+k_points_grid[i+1,j+1][3]
                        weight_a[2]=k_points_grid[i,j+1][3]+k_points_grid[i+1,j+1][3]
                        weight_a[3]=k_points_grid[i,j][3]+k_points_grid[i,j+1][3]
                        position_elements_to_refine=[i*n1+j,i*n1+j+1,(i+1)*n1+j+1,(i+1)*n1+j]
                    elif i < n0 -1 and j == n1-1:
                        u[0]=us[i][j]
                        u[3]=us[i][0]
                        u[2]=us[i+1][0]
                        u[1]=us[i+1][j]
                        weight_a[0]=k_points_grid[i,j][3]+k_points_grid[i+1,j][3]
                        weight_a[1]=k_points_grid[i+1,j][3]+k_points_grid[i+1,0][3]
                        weight_a[2]=k_points_grid[i,0][3]+k_points_grid[i+1,0][3]
                        weight_a[3]=k_points_grid[i,j][3]+k_points_grid[i,0][3]
                        position_elements_to_refine=[i*n1+j,i*n1,(i+1)*n1,(i+1)*n1+j]
                    elif i == n0 -1 and j < n1-1:
                        u[0]=us[i][j]
                        u[3]=us[i][j+1]
                        u[2]=us[0][j+1]
                        u[1]=us[0][j]
                        weight_a[0]=k_points_grid[i,j][3]+k_points_grid[0,j][3]
                        weight_a[1]=k_points_grid[0,j][3]+k_points_grid[0,j+1][3]
                        weight_a[2]=k_points_grid[i,j+1][3]+k_points_grid[0,j+1][3]
                        weight_a[3]=k_points_grid[i,j][3]+k_points_grid[i,j+1][3]
                        position_elements_to_refine=[i*n1+j,i*n1+j+1,j+1,j]
                    else:
                        u[0]=us[i][j]
                        u[3]=us[i][0]
                        u[2]=us[0][0]
                        u[1]=us[0][j]
                        weight_a[0]=k_points_grid[i,j][3]+k_points_grid[0,j][3]
                        weight_a[1]=k_points_grid[0,j][3]+k_points_grid[0,0][3]
                        weight_a[2]=k_points_grid[i,0][3]+k_points_grid[0,0][3]
                        weight_a[3]=k_points_grid[i,j][3]+k_points_grid[i,0][3]
                        position_elements_to_refine=[i*n1+j,i*n1,0,j]
                    ##in the not-degenerate case the Berry curvature is straightforward
                    if number_eigenspaces == self.N:
                        phase = np.asarray(list(map(
                            lambda x,xw,y,yw,z,zw,w,ww: 
                            np.angle(np.asarray(x))*xw
                            + np.angle(np.asarray(y))*yw
                            - np.angle(np.asarray(z))*zw
                            - np.angle(np.asarray(w))*ww,
                            (u[0] @ u[1]),weight_a[0],
                            (u[1] @ u[2]),weight_a[1],
                            (u[3] @ u[2]),weight_a[3],
                            (u[0] @ u[3]),weight_a[0])))
                        chern_number+=phase
                        ##in the degenerate case the Berry curvature is calculated using the non-abelian formulation
                    else:
                        for n in range(number_eigenspaces):
                            bands=[magnonic_branches==n]
                            phase[n] = np.asarray(list(map(
                            lambda x,xw,y,yw,z,zw,w,ww: 
                            np.angle(np.asarray(x))*xw
                            + np.angle(np.asarray(y))*yw
                            - np.angle(np.asarray(z))*zw
                            - np.angle(np.asarray(w))*ww,
                            (u[0][np.ix_(bands,bands)] @ u[1][np.ix_(bands,bands)]),weight_a[0],
                            (u[1][np.ix_(bands,bands)] @ u[2][np.ix_(bands,bands)]),weight_a[1],
                            (u[3][np.ix_(bands,bands)] @ u[2][np.ix_(bands,bands)]),weight_a[3],
                            (u[0][np.ix_(bands,bands)] @ u[3][np.ix_(bands,bands)]),weight_a[0])))
                        chern_number+=phase
    
                        #saving all divergent points: if the refinment procedure is applied, the grid will be refined around these points
                        if np.max(phase)>threshold_dynamical_refinment:
                            position_elements_to_refine_list.extend(position_elements_to_refine)
    
                ###dynamical refinment
                ###the degeneracy is not checked again
                if dynamical_refinment==True:
                    (k_points_list,new_k_grid,n0,n1)=dynamical_refinment(k_points_list,position_elements_to_refine_list,brillouin_primitive_vectors,dynamical_refinment_iteration,symmetries,threshold_k_point,True,None)
                    _,us,magnonic_branches,number_eigenspaces=self.magnonic_surfaces(n0,n1,new_k_grid,False,True,threshold_omega)
                    cycling_through=False
                else:
                    cycling_through=False
    
            return chern_number,magnonic_branches,number_eigenspaces
    
        
    AGGI    UNGERE CALCOLO ATTRAVERSO TRIANGOLARIZZAZIONE
    
    

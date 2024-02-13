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
from radtools.crystal.kpoints import k_points_grid_generator_2D
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
       the Berry Curvature is calculated locally as the line integral of the Berry Connection over a little parallelogram
       this formulation preserve the Berry Curvature Gauge invariance

       Parameters
       ----------
       the ones needed in the production of the 2D k points grid (through the refinment procedure) plus the ones for the Berry Curvature calculation itself

       
       
       threshold_omega |float| the threshold to check degeneracies in magnonic bands
       added_refinment_iteration |int| number of refinment iteration in case of a phase increase over the phase threshold
       number_iterations_dynamical_refinment |int| number of times the dynamical refinment is applied
       Returns
       ----
       chern numbers double |array_like|_
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
        self, n0, n1, k_points_grid, noeigenvectors, nocheckdegeneracy, threshold_omega
    ):
        omegas=np.zeros((n0,n1,self.N),dtype=complex)
        us=np.zeros((n0,n1,self.N,self.N),dtype=complex)
        number_eigenspaces = self.N
        magnonic_branches = np.asarray([i for i in range(self.N)])
        
        if (nocheckdegeneracy == True) and (noeigenvectors == True):
            for i in range(n0):
                for j in range(n1):
                    omegas[i,j]=self.dispersion(k_points_grid[i,j][:3],False)
                    omegas[i,j]=omegas[i,j]*k_points_grid[i,j][3]
        elif (nocheckdegeneracy == True) and (noeigenvectors == False):
            us=np.zeros((n0,n1,self.N,self.N),dtype=complex)
            for i in range(n0):
                for j in range(n1):
                    omegas[i,j],us[i,j]=self.dispersion(k_points_grid[i,j][:3],True)
                    omegas[i,j]=omegas[i,j]*k_points_grid[i,j][3]
                    us[i,j]=us[i,j]*k_points_grid[i,j][3]
        elif (nocheckdegeneracy == False) and (noeigenvectors == True):
            degeneracy_matrix = np.zeros((self.N,self.N),dtype=bool)  
            for i in range(n0):
                for j in range(n1):
                    omegas[i,j]=self.dispersion(k_points_grid[i,j][:3],False)
                    omegas[i,j]=omegas[i,j]*k_points_grid[i,j][3]
                    degeneracy=0
                    for r in range(0,self.N-1):
                        for s in range(r+1,self.N):
                            degeneracy_matrix[i][j] = False
                            degeneracy_matrix[i][j] = np.isclose(omegas[i,j][r],omegas[i,j][s],ato=threshold_omega)
                            if degeneracy_matrix[i][j] == True:
                                degeneracy+=1
                    if degeneracy != 0:
                        magnonic_branches,number_eigenspaces=update_magnonic_branches_degeneracies(
                            degeneracy_matrix, self.N, magnonic_branches, number_eigenspaces)     
        else:
            degeneracy_matrix = np.zeros((self.N,self.N),dtype=bool)
            for i in range(n0):
                for j in range(n1):
                    omegas[i,j],us[i,j]=self.dispersion(k_points_grid[i,j][:3],True)
                    omegas[i,j]=omegas[i,j]*k_points_grid[i,j][3]
                    us[i,j]=us[i,j]*k_points_grid[i,j][3]
                    degeneracy=0
                    for r in range(0,self.N-1):
                        for s in range(r+1,self.N):
                            degeneracy_matrix[i][j] = False
                            degeneracy_matrix[i][j] = np.isclose(omegas[i,j][r],omegas[i,j][s],ato=threshold_omega)
                            if degeneracy_matrix[i][j] == True:
                                degeneracy+=1
                    if degeneracy != 0:
                        magnonic_branches,number_eigenspaces=update_magnonic_branches_degeneracies(
                            degeneracy_matrix, self.N, magnonic_branches, number_eigenspaces)     
        
        return omegas,us,magnonic_branches,number_eigenspaces

    def berry_curvature(self,brillouin_primitive_vectors,
        chosen_plane,
        initial_grid_spacing,
        shift_in_plane,
        shift_in_space,
        symmetries,
        threshold_k_point,
        refinment,
        refinment_spacing,
        refinment_iteration,
        number_elements_interpolation,
        threshold_omega,
        ###dynamical refinment still requires some work...
        dynamical_refinment,
        dynamical_refinment_iteration,
        threshold_dynamical_refinment
    ):
        ## to use a refined grid can produce non integer Chern numbers
        ## however the singularity of the Berry Curvature is sufficient to guess the corresponding Chern number
        k_points_list,k_points_grid,n0,n1=k_points_grid_generator_2D(
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
            number_elements_interpolation,
            threshold_omega
        )

        _,us,magnonic_branches,number_eigenspaces=self.magnonic_surfaces(n0,n1,k_points_grid,False,False,threshold_omega)
        
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

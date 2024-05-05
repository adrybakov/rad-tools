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
from radtools import MagnonDispersion
from scipy.spatial.transform import Rotation
###from radtools.crystal.kpoints import dynamical_refinment_little_paths_2D
from radtools.crystal.kpoints import k_points_generator_2D
from radtools.crystal.kpoints import printing_covering_BZ_2D
from radtools.crystal.kpoints import dynamical_refinment_little_paths_2D
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

def printing_berry_curvature(
        number_vertices,
        phases,
        number_eigespaces,
        refined_k_points_list,
        little_paths,
        brillouin_primitive_vectors_2d,
        file_writing
    ):
        number_little_paths=little_paths.shape[0]
       
        with open(file_writing,'w') as file:
            for i in range(number_little_paths):
                print("litte paths ",i/number_little_paths)
                k_points_around=np.zeros(3,dtype=float)
                for j in range(number_vertices):
                    k_points_around+=refined_k_points_list[little_paths[i][j][0],:3]
                    if little_paths[i][j][1]!=0.0:
                        if little_paths[i][j][1]==3.0:
                            k_points_around+=brillouin_primitive_vectors_2d[0]+brillouin_primitive_vectors_2d[1]
                        else:
                            k_points_around+=brillouin_primitive_vectors_2d[int(little_paths[i][j][1]-1)]
                k_points_around=k_points_around/number_vertices
                file.write(str(k_points_around[0])+" "+str(k_points_around[1])+" "+str(k_points_around[2])+" ")
                for j in range(number_eigespaces):
                    file.write(str(phases[i][j])+" ")
                file.write("\n")

    
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
        self.dispersion=MagnonDispersion(self.spinham,None,None,nodmi=self.nodmi,noaniso=self.noaniso,custom_mask=None)
        self.N=self.dispersion.N

    def magnonic_surfaces(
        self, k_points_list, noeigenvectors, nocheckdegeneracy, threshold_omega
    ):
        number_k_points=len(k_points_list)
        omegas=np.zeros((number_k_points,self.N),dtype=complex)
        us=np.zeros((number_k_points,self.N,self.N),dtype=complex)
        number_eigenspaces = self.N
        magnonic_branches = np.asarray([i for i in range(self.N)])
    
        if (nocheckdegeneracy == True) and (noeigenvectors == True):
            for i in range(number_k_points):
                omegas[i]=self.dispersion.omega(k_points_list[i][:3],True,False)
        elif (nocheckdegeneracy == True) and (noeigenvectors == False):
            us=np.zeros((number_k_points,self.N,self.N),dtype=complex)
            for i in range(number_k_points):
                print("magn",i/number_k_points)
                omegas[i],us[i]=self.dispersion.omega(k_points_list[i][:3],False,False)
        elif (nocheckdegeneracy == False) and (noeigenvectors == True):
            degeneracy_matrix = np.zeros((self.N,self.N),dtype=bool)  
            for i in range(number_k_points):
                omegas[i]=self.dispersion.omega(k_points_list[i][:3],True,False)
                degeneracy=0                
                for r in range(0,self.N-1):
                    for s in range(r+1,self.N):
                        degeneracy_matrix[i] = False
                        degeneracy_matrix[i] = np.isclose(omegas[i][r],omegas[i][s],atol=threshold_omega)
                        if degeneracy_matrix[i] == True:
                            degeneracy+=1
                if degeneracy != 0:
                    magnonic_branches,number_eigenspaces=update_magnonic_branches_degeneracies(
                        degeneracy_matrix, self.N, magnonic_branches, number_eigenspaces)     
        else:
            degeneracy_matrix = np.zeros((self.N,self.N),dtype=bool)
            for r in range(0,self.N):
                    for s in range(self.N):
                        degeneracy_matrix[r][s] = False
            for i in range(number_k_points):
                omegas[i],us[i]=self.dispersion.omega(k_points_list[i][:3],False,False)
                degeneracy=0
                for r in range(0,self.N-1):
                    for s in range(r+1,self.N):
                        degeneracy_matrix[r][s] = False
                        degeneracy_matrix[r][s] = np.isclose(omegas[i][r],omegas[i][s],atol=threshold_omega)
                        if degeneracy_matrix[r][s] == True:
                            degeneracy+=1
                if degeneracy != 0:
                    magnonic_branches,number_eigenspaces=update_magnonic_branches_degeneracies(
                        degeneracy_matrix, self.N, magnonic_branches, number_eigenspaces)     
        
        return omegas,us,magnonic_branches,number_eigenspaces

    def circuitation_over_a_path(self,u_values_along_the_path,k_point_weights_along_the_path,number_modes,number_eigenspaces,magnonic_branches,number_vertices
    ):
        print(u_values_along_the_path,number_eigenspaces,magnonic_branches,k_point_weights_along_the_path)
        print(u_values_along_the_path[0][0])
        phases=np.zeros(number_eigenspaces,dtype=complex)
        if number_eigenspaces == number_modes:
            for n in range(number_eigenspaces):
                count=0
                while count<number_vertices-1:
                    phases[n]+=np.angle(np.asarray(u_values_along_the_path[count,n]@u_values_along_the_path[count+1,n]))*k_point_weights_along_the_path[count]
                    count+=1
                phases[n]+=np.angle(np.asarray(u_values_along_the_path[number_vertices-1,n]@u_values_along_the_path[0,n]))*k_point_weights_along_the_path[number_vertices-1]  
                ###print(phases)
        ### in the degenerate case the Berry curvature is calculated using the Non-Abelian formulation
        else:
           ## print("dentro degenerate procedure")
            for n in range(number_eigenspaces):
                bands=[magnonic_branches==n]
                count=0
                while count<number_vertices-1:
                    phases[n]+=np.angle(np.asarray(u_values_along_the_path[count][np.ix_(bands,bands)] @ u_values_along_the_path[count+1][np.ix_(bands,bands)]))*k_point_weights_along_the_path[count]
                    count+=1
                phases[n]+=np.angle(np.asarray(u_values_along_the_path[number_vertices-1][np.ix_(bands,bands)] @ u_values_along_the_path[0][np.ix_(bands,bands)]))*k_point_weights_along_the_path[number_vertices-1]
        return phases

    
    def berry_curvature(
        self,
        triangulation,
        brillouin_primitive_vectors_3d,
        chosen_plane,
        grid_spacing,
        default_gridding,
        shift_in_plane,
        shift_in_space,
        symmetries,
        threshold_symmetry,
        refinment_spacing,
        refinment_iterations,
        threshold_omega,
        dynamical_refinment_flag,
        threshold_dynamical_refinment
    ):  
        number_modes=self.N

        refined_k_points_list,little_paths=k_points_generator_2D(
                                                    brillouin_primitive_vectors_3d,
                                                    brillouin_primitive_vectors_3d,
                                                    chosen_plane,
                                                    grid_spacing,
                                                    True,
                                                    default_gridding,
                                                    False,
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
                little_paths,
                4
            )
        
        brillouin_primitive_vectors_2d=np.zeros((2,3),dtype=float)
        normalized_brillouin_primitive_vectors_2d=np.zeros((2,3),dtype=float)
        count = 0
        for i in range(3):
            if chosen_plane[i]!=0:
                brillouin_primitive_vectors_2d[count] = brillouin_primitive_vectors_3d[i]
                normalized_brillouin_primitive_vectors_2d[count] = brillouin_primitive_vectors_2d[count]/(brillouin_primitive_vectors_2d[count]@brillouin_primitive_vectors_2d[count])
                count += 1
        
        print("magnonic_surfaces ")
        _,us,magnonic_branches,number_eigenspaces=self.magnonic_surfaces(refined_k_points_list,False,True,threshold_omega)
        print(magnonic_branches)
        
        number_k_points=refined_k_points_list.shape[0]
        with open("magnonic_surfaces.data",'w') as file:
            for i in range(number_k_points):
                file.write(str(refined_k_points_list[i][0])+" "+str(refined_k_points_list[i][1])+" "+str(refined_k_points_list[i][2])+" ")
                for j in range(number_eigenspaces):
                    file.write(str(_[i][j].real)+" ")
                file.write("\n")

        if triangulation==True:
            number_vertices=3
        else:
            number_vertices=4

        u=np.zeros((number_vertices,number_modes,number_modes),dtype=complex)
        weights=np.zeros(number_vertices,dtype=float)
        old_chern=np.zeros(number_eigenspaces,dtype=float)
        chern=np.zeros(number_eigenspaces,dtype=float)

        dynamical_refinment=True
        while dynamical_refinment==True:
            number_little_paths=little_paths.shape[0]
            _,us,magnonic_branches,number_eigenspaces=self.magnonic_surfaces(refined_k_points_list,False,False,threshold_omega)
            i=0
            phases=np.zeros((number_little_paths,number_eigenspaces),dtype=float)
            while i < number_little_paths:
               # print(i/number_little_paths)
                for j in range(number_vertices):
                    print(i,j,little_paths[i][j])
                    u[j]=us[little_paths[i][j][0]]
                    weights[j]=1.0/number_k_points
                phases[i]=self.circuitation_over_a_path(u,weights,number_modes,number_eigenspaces,magnonic_branches,number_vertices)
                chern+=phases[i]
                i+=1
            
            dynamical_refinment=False

            #### condition for dynamical refinment (checking that the Chern numbers are integer numbers)
            if np.any(np.asarray(list(map(lambda x: float(x).is_integer(),chern)))) == False:
                average_phase=np.mean(phases)
                old_chern+=chern
                max_indices=[i for i in range(number_little_paths) if np.mean(phases[i])>average_phase-threshold_dynamical_refinment or np.mean(phases[i])<average_phase+threshold_dynamical_refinment]
                selected_little_paths=little_paths[max_indices,:]
              
                ###refined_k_points_list,little_paths=dynamical_refinment_little_paths_2D(
                ###                                                little_paths,
                ###                                                refined_k_points_list,
                ###                                                selected_little_paths,
                ###                                                number_vertices,
                ###                                                refinment_iterations,
                ###                                                symmetries,
                ###                                                threshold_symmetry,
                ###                                                brillouin_primitive_vectors_3d,
                ###                                                chosen_plane,
                ###                                                brillouin_primitive_vectors_2d,
                ###                                                normalized_brillouin_primitive_vectors_2d,
                ###                                                0.0000000000001,
                ###                                                default_gridding
                ###                                            )
                if dynamical_refinment_flag==True:
                    dynamical_refinment=False
            else:
                dynamical_refinment=False
        
        return number_eigenspaces,phases,refined_k_points_list,little_paths,brillouin_primitive_vectors_2d,chern

### TESTING INPUT
if __name__ == "__main__":
    import timeit
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    from termcolor import cprint
    from radtools.io.internal import load_template
    from radtools.io.tb2j import load_tb2j_model
    from radtools.magnons.dispersion import MagnonDispersion
    from radtools.decorate.stats import logo
    from radtools.spinham.constants import TXT_FLAGS
    from radtools.decorate.array import print_2d_array
    from radtools.decorate.axes import plot_hlines

    spin={"Cr1":[0,0,1],"Cr2":[0,0,1]}
    input_filename="/home/marco-marino/rad-tools/exchange.out"
    spinham = load_tb2j_model(input_filename)
    for key,values in spin.items():
        atom_name=key
        atom=spinham.get_atom(atom_name)
        atom_spin=list(values)
        atom.spin_vector=atom_spin

    primitive_vectors_3d=np.zeros((3,3),dtype=float)
    primitive_vectors_3d[0]=[7.176,0,0]
    primitive_vectors_3d[1]=[-3.588,6.215,0]
    primitive_vectors_3d[2]=[0,0,21.000]
    volume=np.dot(np.cross(primitive_vectors_3d[0],primitive_vectors_3d[1]),primitive_vectors_3d[2])/(2*np.pi)

    brillouin_primitive_vectors_3d=np.zeros((3,3),dtype=float)
    brillouin_primitive_vectors_3d[0]=np.cross(primitive_vectors_3d[1],primitive_vectors_3d[2])/volume
    brillouin_primitive_vectors_3d[1]=np.cross(primitive_vectors_3d[2],primitive_vectors_3d[0])/volume
    brillouin_primitive_vectors_3d[2]=np.cross(primitive_vectors_3d[0],primitive_vectors_3d[1])/volume
    print(brillouin_primitive_vectors_3d)

    chosen_plane=[1,1,0]
    symmetries=[[0,0,0]]
    grid_spacing=0.1
    refinment_iterations=3
    refinment_spacing=0.001
    threshold_symmetry=0.001
    threshold_minimal_refinment=0.000000001
    threshold_degeneracy=0.000001
    default_gridding=100
    count=0
    shift_in_plane=[0,0]
    shift_in_space=[0,0,0]
    covering_BZ=True
    nodmi=False
    noaniso=False

    kp = spinham.kpoints
    fig, ax = plt.subplots()
    dispersion = MagnonDispersion(spinham,None,None,nodmi,noaniso,None)
    omegas = dispersion(kp,True)
    ax.set_xticks(kp.coordinates(), kp.labels, fontsize=15)
    ax.set_ylabel("E, meV", fontsize=15)
    for omega_tmp in omegas:
        ax.plot(kp.flatten_points(), omega_tmp)
    
    ####plt.show()
    triangulation=False
    dynamical_refinment_flag=False
    threshold_dynamical_refinment=0.000001
    berry_curvature=Berry_curvature(spinham,nodmi,noaniso)
    number_vertices=4
    number_eigenspaces,phases,refined_k_points_list,little_paths,brillouin_primitive_vectors_2d,chern=berry_curvature.berry_curvature(triangulation,brillouin_primitive_vectors_3d,chosen_plane,
                                    grid_spacing,default_gridding,shift_in_plane,shift_in_space,
                                    symmetries,threshold_symmetry,refinment_spacing,refinment_iterations,
                                    threshold_degeneracy,dynamical_refinment_flag,threshold_dynamical_refinment)
    
   ## file_writing="berry_curvature.dat"
   ## printing_berry_curvature(number_vertices,phases,number_eigenspaces,refined_k_points_list,little_paths,brillouin_primitive_vectors_2d,file_writing)
   ## print(number_eigenspaces,phases,chern)
   ###set hideen3d
   ###set dgrid3d 50,50 qnorm 2
   ###splot '' w lines

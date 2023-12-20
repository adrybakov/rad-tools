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
from radtools import MagnonDispersion
from scipy.spatial.transform import Rotation
from radtools.crystal.kpoints import *
from radtools.crystal.kpoints import k_points_grid_2d_refinment_and_symmetry
from radtools.crystal.kpoints import dynamical_refinment
from radtools.crystal.kpoints import check_inside_brillouin_zone
from radtools.crystal.kpoints import mapping_to_square_grid
from radtools.geometry import span_orthonormal_set
from radtools.magnons.diagonalization import ColpaFailed, solve_via_colpa
from radtools.spinham.hamiltonian import SpinHamiltonian
from radtools.magnons import dispersion

def function_update_list_magnons(
        check_degeneracy, number_magnons, **magnonic_branches
    ):
        ##if there is only one element, a study of degeneracies does not make any sense
        if len(magnonic_branches) == 1:
            return magnonic_branches
        magnonic_branches_tmp = {}
        for eigenspace, eigenvalues in magnonic_branches.items():
            magnonic_branches_tmp[str(eigenspace)] = magnonic_branches[str(eigenspace)]
        
        for i in range(0, number_magnons - 1):
            for j in range(i + 1, number_magnons):
                ##if during the checks all the eigenvalues result degenerate, then a study of degeneracies is not anymore needed
                if len(magnonic_branches_tmp) == 1:
                    break
                else:
                    check_changes_branches = 0  
                    magnonic_branches_new = {}
                    if check_degeneracy[i][j] == True:
                        check_degeneracy_ij = 0 
                        ##associating the two degenerate eigenvalues to their eigenspaces 
                        for eigenspace, eigenvalues in magnonic_branches_tmp.items():
                            for eigenvalue in eigenvalues:
                                if eigenvalue == i:
                                    positioni = eigenspace
                                    check_degeneracy_ij = check_degeneracy_ij + 1
                                if eigenvalue == j:
                                    positionj = eigenspace
                                    check_degeneracy_ij = check_degeneracy_ij + 1
                        if check_degeneracy_ij == 2:
                            ### checking if the eigenvalues are already in the same eigenspace
                            if positioni != positionj:
                                count = 0
                                ### saving the new eigenspaces
                                for (
                                    eigenspace,
                                    eigenvalues,
                                ) in magnonic_branches_tmp.items():
                                    if eigenspace == positioni:
                                        magnonic_branches_new[str(count)] = (
                                            eigenvalues
                                            + magnonic_branches_tmp[
                                                str(positionj)
                                            ]
                                        )
                                        count = count + 1
                                    else:
                                        if eigenspace != positionj:
                                            magnonic_branches_new[
                                                str(count)
                                            ] = eigenvalues
                                            count = count + 1
                            else:
                                check_changes_branches = 1
                        else:
                            check_changes_branches = 1
                    else:
                        check_changes_branches = 1
                    ### saving the new eigenspaces in order to compare with the new degeneracies at the next cycle
                    if check_changes_branches != 1:
                        magnonic_branches_tmp = {}
                        for eigenspace, eigenvalues in magnonic_branches_new.items():
                            magnonic_branches_tmp[
                                str(eigenspace)
                            ] = magnonic_branches_new[str(eigenspace)]
                    else:
                        for eigenspace, eigenvalues in magnonic_branches_tmp.items():
                            magnonic_branches_new[
                                str(eigenspace)
                            ] = magnonic_branches_tmp[str(eigenspace)]
        ### ordering the eigenvalues inside each eigenspace
        for eigenspace, eigenvalues in magnonic_branches_new.items():
            magnonic_branches_new[str(eigenspace)] = sorted(
                set(magnonic_branches_new[str(eigenspace)])
            )
        return magnonic_branches_new

       # Following a few easy functions for easy use

# function returning the minimum difference between adjacent elements in an array of n elements
def function_minimum_difference_pair(a, n):
    a = sorted(a)
    diff = 10**20
    for i in range(n - 1):
        if a[i + 1] - a[i] < diff:
            diff = a[i + 1] - a[i]
    return diff
# function: scalar product between rows of square matrices
# a and b are matrices, c is a matrix which rows are the product of the a and b rows
def function_product_of_rows_eigenvector_matrix(a, b):
    c = np.zeros(len(a[:, 0]), dtype=complex)
    # its tricky, but try to use numpy here
    for i in range(0, a.shape[0]):
        for j in range(0, a.shape[1]):
            c[i] = c[i] + np.conj(a[i, j]) * b[i, j]
    return c

def function_generate_indipendent_list(n):
            list_magnonic_branches_tmp = {}
            for i in range(0, n):
                list_magnonic_branches_tmp[str(i)] = [i]
            return list_magnonic_branches_tmp

# function: scalar product and phase extraction in a eigenspace with multiplicity greater than zero
# a and b are matrices, values is a list of rows of a and b
# the rows are the generator of the degenerate eigenspace
def function_phase_eigenspace_indicated(a, b, values):
    a_tmp = np.zeros((len(values), len(values)), dtype=complex)
    b_tmp = np.zeros((len(values), len(values)), dtype=complex)
    # No need to initialize c
    c = np.zeros((len(values), len(values)), dtype=complex)
    # building submatrices
    for i, eli in enumerate(values):
        for j, elj in enumerate(values):
            a_tmp[i, j] = a[eli, elj]
            b_tmp[i, j] = b[eli, elj]
    # multiplying submatrices
    c = np.conj(a) * b
    return np.angle(np.linalg.det(c))
# Use docstrings to document your code
# function: updating information about magnonic branches taking into acount the new/old found degeneracies
# the old list of magnonic branches is given, togheter with the degenerate eigenvalues and the number of magnetic atoms (or magnonic branches)
# grouping magnonic branches with respect to degeneracies
# dictionary: each key is associated to an eigenspace, each value is a list of degenerate branches
# list_magnonic_branches.={ 1: [1,2,3], 2: [4,5]..} magnonic branches 1 2 and 3 are degenerate...

class Berry_curvature:
    r"""
       berry_curvature in a 2D plane
       the berry curvature is calculated locally as the line integral of the berry connection
       this formulation preserve the berry curvature gauge invariance

       Parameters
       ----------
       the ones needed in the production of the 2D k-grid plane plus the ones for the Berry curvature calculation itself
       
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
    brillouin_primitive_vectors,
    plane_2d,
    grid_spacing,
    shift_in_plane,
    shift_in_space,
    symmetry,
    refinment_spacing,
    refinment_iteration,
    threshold_k_grid,
    threshold_omega,
    added_refinment_iteration,
    number_iterations_dynamical_refinment,
    ) -> None:
        self.spinham=spinham
        self.nodmi=nodmi
        self.noaniso=noaniso
        self.brillouin_primitive_vectors=brillouin_primitive_vectors
        self.plane_2d=plane_2d
        self.grid_spacing=grid_spacing
        self.shift_in_plane=shift_in_plane
        self.shift_in_space=shift_in_space
        self.symmetry=symmetry
        self.refinment_spacing=refinment_spacing
        self.refinment_iteration=refinment_iteration
        self.threshold_k_grid=threshold_k_grid
        self.threshold_omega=threshold_omega
        self.added_refinment_iteration=added_refinment_iteration
        self.number_iterations_dynamical_refinment=number_iterations_dynamical_refinment
        self.dispersion=MagnonDispersion(self.spinham,nodmi=self.nodmi,noaniso=self.noaniso)
        self.N=self.dispersion.N
    
    ##function calculating mangonic energies and eigenvectors and checking any degeneracy; returning also the proper eigenspaces
    def calculation_omega_and_u_checking_degeneracy(
        self, k_point, threshold_omega, noeigenvectors, nocheckdegeneracy, **list_magnonic_branches
    ):
        # omit passing keyword arguments as positional ones
        # It reduces readability of the code
        # if eigenvectors are not required, the extraction is avoided
        if noeigenvectors==False:
            omega_k, u_k = self.dispersion.omega(k_point,False)
        else:
            omega_k = self.dispersion.omega(k_point,True)
        if nocheckdegeneracy == False:
            degeneracy = 0
            check_degeneracy = np.zeros((len(omega_k), len(omega_k)), dtype=bool)
            for i in range(0, len(omega_k) - 1):
                for j in range(i + 1, len(omega_k)):
                    check_degeneracy[i][j] = False
                    check_degeneracy[i][j] = np.isclose(
                        omega_k[i], omega_k[j], atol=threshold_omega
                    )
                    if check_degeneracy[i][j] == True:
                        # print(omega_k[i],omega_k[j])
                        degeneracy = degeneracy + 1
                        # print(degeneracy)
            if degeneracy != 0:
                list_magnonic_branches = function_update_list_magnons(
                    check_degeneracy, len(omega_k), **list_magnonic_branches
                )
            else:
                if not list_magnonic_branches:
                    for i in range(0, len(omega_k)):
                        list_magnonic_branches[str(i)] = [i]
        if noeigenvectors==False:
            return list_magnonic_branches, omega_k, u_k
        else:
            return list_magnonic_branches, omega_k


    def berry_curvauture_calculation(self):
        # generating k points grid with refinment as already pointed out
        (
            normalized_chosen_reciprocal_plane,
            k0,
            k1,
            k_points_grid_2d
        ) = k_points_grid_2d_refinment_and_symmetry(
            self.brillouin_primitive_vectors,
            self.plane_2d,
            self.grid_spacing,
            self.shift_in_plane,
            self.shift_in_space,
            self.symmetry,
            self.refinment_spacing,
            self.refinment_iteration,
            self.threshold_k_grid,
        )
        # reducing the k points grid to a square grid
        square_grid=mapping_to_square_grid(
            k_points_grid_2d,
            normalized_chosen_reciprocal_plane,
            self.refinment_spacing
        )
        dimx=len(square_grid)

        #checking if there are degeneracies and saving the magnonic surfaes in the file magnonic_surfaces.txt
        list_magnonic_branches = function_generate_indipendent_list(self.N)
        u_k=np.zeros((dimx,dimx,self.N,self.N), dtype=complex)
        weight_k=np.zeros((dimx,dimx))
        ##file_2 = open("magnonic_surfaces.txt", "w")
        for i in range(0, dimx):
            for j in range(0, dimx):
                ##print(
                ##    "degeneracy studies",
                ##    int(i / dimx * 100),
                ##    "%  and ",
                ##    int(j / dimx * 100),
                ##    "%",
                ##    end="\r",
                ##    flush=True,
                ##) 
                k_point_tmp = np.asarray(square_grid[i][j][:3])
                weight_k[i][j] = square_grid[i][j][3]
                (
                    list_magnonic_branches,
                    omega_k_tmp,
                    u_k[i][j]
                ) = self.calculation_omega_and_u_checking_degeneracy(
                    k_point_tmp,
                    self.threshold_omega,
                    True,
                    True,
                    **list_magnonic_branches,
                )
                u_k[i][j]=u_k[i][j]*weight_k[i][j]
                ##file_2.write(f"{k_point_tmp} {omega_k_tmp}   \n")
        print(list_magnonic_branches)
        
        ##calculating the berry curvature in each k point directly
        ##taking into account possible degeneracies (in the non-abelian formulation)
        phase=np.zeros((len(list_magnonic_branches)), dtype=complex)
        phase_tmp=np.zeros((len(list_magnonic_branches)), dtype=complex)
          
        ### here the berry calculation
        ### to each k point of the largest grid (i,j) is associated a refinment list (r)
        ### for each point of the refinment list (r) are considered 4 k points around it and on these 4 k points it is calculated the line integral of the berry connection
        ### if there is any degeneracy in the magnonic branches the line integral is properly formulated (degeneracy detected)
        ### the berry curvature (the line integral of the berry connection) is properly summed (k point weight) over the refinment list (r)
        ### if the berry curvature increases between two successive k points of the ij grid, a dynamical refinment is considered (to increase the precision around points with high increments in the berry curvature)
        ### the number of dynamical refinments is fixed by number_iterations_dynamical_refinment
        ### all these data are then saved in the respective files

        ### the procedure above introduces some phases that are not cancelled, so a better strategy is the following 
        ### mapping all the ij points to a square grid 
        ### the dynamical refinment is a little bit more complex here (but can be implemented as well)
       ## file_1 = open("berry_curvature.txt", "w") 
        ### here the eigenvectors are calculated on the way (the static procedure can be applied as well)
        for i in range(0, dimx):
            for j in range(0, dimx):
              ## print(
              ##     "berry grid",
              ##     int(i / k0 * 100),
              ##     "%  and ",
              ##     int(j / k1 * 100),
              ##     "%",
              ##     end="\r",
              ##     flush=True,
              ## )
              ## ### considering the little square path
                u_k_tmp=np.zeros((4,self.N,self.N))
                if i < dimx -1 and j < dimx -1:
                    u_k_tmp[0]=u_k[i][j]
                    u_k_tmp[3]=u_k[i][j+1]
                    u_k_tmp[2]=u_k[i+1][j+1]
                    u_k_tmp[1]=u_k[i+1][j]
                elif i < dimx -1 and j == dimx-1:
                    u_k_tmp[0]=u_k[i][j]
                    u_k_tmp[3]=u_k[i][0]
                    u_k_tmp[2]=u_k[i+1][0]
                    u_k_tmp[1]=u_k[i+1][j]
                elif i == dimx -1 and j < dimx-1:
                    u_k_tmp[0]=u_k[i][j]
                    u_k_tmp[3]=u_k[i][j+1]
                    u_k_tmp[2]=u_k[0][j+1]
                    u_k_tmp[1]=u_k[0][j]
                else:
                    u_k_tmp[0]=u_k[i][j]
                    u_k_tmp[3]=u_k[i][0]
                    u_k_tmp[2]=u_k[0][0]
                    u_k_tmp[1]=u_k[0][j]
                ##in the not-degenerate case the Berry curvature is straightforwardath
                if len(list_magnonic_branches) == self.N:
                    phase_tmp = np.asarray(
                        list(
                            map(
                                lambda x, y, z, w: np.angle(
                                    np.asarray(x)
                                )
                                + np.angle(np.asarray(y))
                                - np.angle(np.asarray(z))
                                - np.angle(np.asarray(w)),
                                function_product_of_rows_eigenvector_matrix(
                                    u_k_tmp[0], u_k_tmp[1]
                                ),
                                function_product_of_rows_eigenvector_matrix(
                                    u_k_tmp[1], u_k_tmp[2]
                                ),
                                function_product_of_rows_eigenvector_matrix(
                                    u_k_tmp[3], u_k_tmp[2]
                                ),
                                function_product_of_rows_eigenvector_matrix(
                                    u_k_tmp[0], u_k_tmp[3]
                                ),
                            )
                        )
                    )
                ##in the degenerate case the Berry curvature is calculated using the non-abelian formulation
                else:
                    for (
                        eigenspace,
                        eigenvalues,
                    ) in list_magnonic_branches.items():
                        phase_tmp[eigenspace] = (
                            np.asarray(
                                function_phase_eigenspace_indicated(
                                    u_k_tmp[0], u_k_tmp[1], eigenvalues
                                )
                            )
                            + np.asarray(
                                function_phase_eigenspace_indicated(
                                    u_k_tmp[1], u_k_tmp[2], eigenvalues
                                )
                            )
                            - np.asarray(
                                function_phase_eigenspace_indicated(
                                    u_k_tmp[3], u_k_tmp[2], eigenvalues
                                )
                            )
                            - np.asarray(
                                function_phase_eigenspace_indicated(
                                    u_k_tmp[0], u_k_tmp[3], eigenvalues
                                )
                            )
                        )
                phase=phase+phase_tmp
    
        return list_magnonic_branches, phase

    def __call__(self, *args, **kwargs):
        return self.omegas(*args, **kwargs)

if __name__ == "__main__":
    import timeit
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    from termcolor import cprint
    import radtools
    from radtools.io.internal import load_template
    from radtools.io.tb2j import load_tb2j_model
    from radtools.magnons.dispersion import MagnonDispersion
    from radtools.decorate.stats import logo
    from radtools.spinham.constants import TXT_FLAGS
    from radtools.decorate.array import print_2d_array
    from radtools.decorate.axes import plot_hlines

    spin={"Cr1":[0,0,1],"Cr2":[0,0,1]}
    input_filename='/home/marcomarino/TestBerry_curvature/exchange.2.out'
    spinham = load_tb2j_model(input_filename)
    for key,values in spin.items():
        atom_name=key
        atom=spinham.get_atom(atom_name)
        atom_spin=list(values)
        atom.spin_vector=atom_spin
    dispersion = MagnonDispersion(spinham,nodmi=False,noaniso=False)
    
    plane_2d=[1,1,0]
    threshold_k_grid=0.0000001
    shift_in_plane=[0,0]
    shift_in_space=[0,0,0]
    symmetry=[[0,0,0]]
    grid_spacing=0.1
    refinment_iteration=0
    refinment_spacing=0.5
    threshold_omega=0.6

    brillouin_primitive_vectors=np.zeros((3,3))
    brillouin_primitive_vectors[0]=[1,0,0]
    brillouin_primitive_vectors[1]=[-1/2,3,0]
    brillouin_primitive_vectors[2]=[0,0,1]
    print(brillouin_primitive_vectors)
    added_refinment_iteration=0
    noaniso=False
    nodmi=False
    number_iterations_dynamical_refinment=1
    berry=Berry_curvature(spinham,nodmi,noaniso,brillouin_primitive_vectors,plane_2d,grid_spacing,shift_in_plane,
                            shift_in_space,symmetry,refinment_spacing,refinment_iteration,threshold_k_grid,
                            threshold_omega,added_refinment_iteration,number_iterations_dynamical_refinment)
    print(berry.berry_curvauture_calculation())
 #sed 's/[][",'"'"']//g' berry_points.txt | awk 'NF>1{print}' > berry_points.cleaned.txt
    #sed 's/[][",'"'"']//g' magnonic_surfaces.txt | awk 'NF>1{print}' > magnonic_surfaces.cleaned.txt
     #sed 's/[][",'"'"']//g' k_points.txt | awk 'NF>1{print}' > k_points.cleaned.txt
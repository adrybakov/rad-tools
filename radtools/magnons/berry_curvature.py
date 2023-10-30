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
from scipy.spatial.transform import Rotation

from radtools.crystal.kpoints import *
from radtools.geometry import span_orthonormal_set
from radtools.magnons.diagonalization import ColpaFailed, solve_via_colpa
from radtools.spinham.hamiltonian import SpinHamiltonian


def berry_curvature(
        self,
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
    ):
        r"""
        2berry_curvature # short header

        # Detailed description

        Parameters
        ----------
        the ones needed in the production of the 2d k grid plut the ones for the berry phase calculation itself

        # There are no double in python
        # Try to follow the format of the docstrings from the package
        threshold_omega |double| the threshold to check degeneracies in magnonic bands
        added_refinment_iteration |int| number of refinment iteration in case of a phase increase over the phase threshold
        number_iterations_dynamical_refinment |int| number of times the dynamical refinment is applied

        Returns
        ----
        chern numbers double |array_like|_
        """
        # Use docstrings to document your code
        ####Following a few easy functions for easy use
        def function_minimum_difference_pair(a, n):
            a = sorted(a)
            diff = 10**20
            for i in range(n - 1):
                if a[i + 1] - a[i] < diff:
                    diff = a[i + 1] - a[i]
            return diff

        # Use docstrings to document your code
        # function: scalar product between rows of square matrices
        # a and b are matrices, c is a matrix which rows are the product of the a and b rows
        def function_product_of_rows_eigenvector_matrix(a, b):
            c = np.zeros(len(a[:, 0]), dtype=complex)
            # Its tricky, but try to use numpy here
            # a.shape[0] is the number of rows
            for i in range(0, len(a[:, 0])):
                # a.shape[1] is the number of columns
                for j in range(0, len(a[0, :])):
                    c[i] = c[i] + np.conj(a[i, j]) * b[i, j]
            return c

        # Use docstrings to document your code
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
        def function_update_list_magnons(
            check_degeneracy, number_magnons, **list_magnonic_branches
        ):
            # print(list_magnonic_branches)
            if len(list_magnonic_branches) == 1:
                return list_magnonic_branches

            # list is a dict? Bad naming
            list_magnonic_branches_tmp = {}
            # key, values - non-pythonic (pythonic example: for fruits, vegetable in coushion:)
            for key, values in list_magnonic_branches.items():
                list_magnonic_branches_tmp[str(key)] = list_magnonic_branches[str(key)]

            for i in range(0, number_magnons - 1):
                for j in range(i + 1, number_magnons):
                    if len(list_magnonic_branches_tmp) == 1:
                        break
                    else:
                        flag1 = 0  # Bad naming
                        list_magnonic_branches_new = {}
                        if check_degeneracy[i][j] == True:
                            flag2 = 0  # Bad naming
                            for key, values in list_magnonic_branches_tmp.items():
                                for value in values:
                                    if value == i:
                                        positioni = key
                                        flag2 = flag2 + 1
                                    if value == j:
                                        positionj = key
                                        flag2 = flag2 + 1
                            if flag2 == 2:
                                if positioni != positionj:
                                    count = 0
                                    for (
                                        key,
                                        values,
                                    ) in list_magnonic_branches_tmp.items():
                                        if key == positioni:
                                            list_magnonic_branches_new[str(count)] = (
                                                values
                                                + list_magnonic_branches_tmp[
                                                    str(positionj)
                                                ]
                                            )
                                            count = count + 1
                                        else:
                                            if key != positionj:
                                                list_magnonic_branches_new[
                                                    str(count)
                                                ] = values
                                                count = count + 1
                                else:
                                    flag1 = 1
                            else:
                                flag1 = 1
                        else:
                            flag1 = 1
                        if flag1 != 1:
                            list_magnonic_branches_tmp = {}
                            for key, values in list_magnonic_branches_new.items():
                                list_magnonic_branches_tmp[
                                    str(key)
                                ] = list_magnonic_branches_new[str(key)]
                        else:
                            for key, values in list_magnonic_branches_tmp.items():
                                list_magnonic_branches_new[
                                    str(key)
                                ] = list_magnonic_branches_tmp[str(key)]

            for key, values in list_magnonic_branches_new.items():
                list_magnonic_branches_new[str(key)] = sorted(
                    set(list_magnonic_branches_new[str(key)])
                )
            return list_magnonic_branches_new

        # Use docstrings to document your code
        ##function calculating mangonic energies and eigenvectors and checking any degeneracy; returning also the proper egenspaces
        def calculation_omega_and_u_checking_degeneracy(
            k_point, threshold_omega, **list_magnonic_branches
        ):
            degeneracy = 0
            # omit passing keyword arguments as positional ones
            # It reduces readability of the code
            omega_k, u_k = self.omega(k_point, True)
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
            return list_magnonic_branches, omega_k, u_k

        # generating k points grid with refinment as already pointed out
        (
            normalized_chosen_reciprocal_plane,
            k0,
            k1,
            k_points_grid_2d
        ) = self.k_points_grid_2d_refinment_and_symmetry(
            plane_2d,
            grid_spacing,
            shift_in_plane,
            shift_in_space,
            symmetry,
            refinment_spacing,
            refinment_iteration,
            threshold_k_grid,
        )

     ## # calculating the berry curvature in each k point directly
     ## # taking into account possible degeneracies (in the non-abelian formulation: pointwise) TO DISCUSS THE THEORETICAL VALIDITY OF THIS
        k_points_path = np.zeros((4, 3))
        u_k_tmp = np.zeros((4, self.N, self.N), dtype=complex)
        omega_k_tmp = np.zeros((4, self.N))

        chern_number_tmp = np.zeros((k0, k1, self.N), dtype=complex)
        chern_number = 0
        phase = np.zeros((self.N), dtype=complex)

        def function_generate_indipendent_list(n):
            list_magnonic_branches_tmp = {}
            for i in range(0, n):
                list_magnonic_branches_tmp[str(i)] = [i]
            return list_magnonic_branches_tmp

        list_magnonic_branches = function_generate_indipendent_list(self.N)
        local_phase = np.zeros((self.N), dtype=complex)
        old_total_phase = np.zeros((self.N), dtype=complex)
        check_tot = 0
        phase_threshold = 0.5
#   
         ## ##here the berry calculation
     ## # to each k point of the largest grid (i,j) is associated a refinment list (r)
     ## # for each point of the refinment list (r) are considered 4 k points around it and on these 4 k points it is calculated the line integral of the berry connection
     ## # if there is any degeneracy in the subset of 4 k points the line integral is properly formulated (degeneracy detected)
     ## # the berry curvature (the line integral of the berry connection) is properly summed (k point weight) over the refinment list (r)
     ## # if the berry curvature increases between two successive k points, a dynamical refinment is considered (to increase the precision around points with high increments in the berry curvature)
     ## # the number of dynamical refinment is fixed by number_iterations_dynamical_refinment
     ## # all these data are then saved in the respective files
#
     ## # Nested with statements reduce the readability of the code
     ## # use:
     ## # file_1 = open("berry_averaged_points.txt", "w")
     ## # file_2 = open("magnonic_surfaces.txt", "w")
     ## # file_3 = open("berry_detailed_points.txt", "w")
     ## # ...
     ## # file_1.close()
     ## # file_2.close()
     ## # file_3.close()

     ## with open("berry_averaged_points.txt", "w") as file_1:
     ##     with open("magnonic_surfaces.txt", "w") as file_2:
     ##         with open("berry_detailed_points.txt", "w") as file_3:
        file_1 = open("berry_averaged_points.txt", "w")
        file_2 = open("magnonic_surfaces.txt", "w")
        file_3 = open("berry_detailed_points.txt", "w")
        for i in range(0, k0):
            for j in range(0, k1):
                print(
                    "berry large grid",
                    int(i / k0 * 100),
                    "%  and ",
                    int(j / k1 * 100),
                    "%",
                )
                print(f"{i,j}")
                flag_local_evaluation = 0
                while (
                    flag_local_evaluation
                    < number_iterations_dynamical_refinment
                ):
                    k_points_list_tmp = k_points_grid_2d[i][j]
                    length_refinment_grid = len(k_points_list_tmp)
                    local_phase = 0
                    total_phase = 0
                    check_tmp = 0
                    phases_list = []
                    kpoints_list = []
                    omegas_list = []
                    for r in range(0, length_refinment_grid):
                        print(
                            "status: {}%".format(
                                int(r / length_refinment_grid * 100)
                            ),
                            end="\r",
                            flush=True,
                        )
                        list_magnonic_branches = (
                            function_generate_indipendent_list(self.N)
                        )
                        k_point_tmp = np.asarray(k_points_list_tmp[r, :3])
                        if length_refinment_grid > 1:
                            berry_spacing = (
                                function_minimum_difference_pair(
                                    k_point_tmp, len(k_point_tmp)
                                )
                                / 2
                            )
                        else:
                            berry_spacing = 1.0e-7
                        k_points_path[0, :] = (
                            k_point_tmp
                            + normalized_chosen_reciprocal_plane[0, :] * berry_spacing
                        )
                        k_points_path[1, :] = (
                            k_point_tmp
                            - normalized_chosen_reciprocal_plane[0, :] * berry_spacing
                        )
                        k_points_path[2, :] = (
                            k_point_tmp
                            + normalized_chosen_reciprocal_plane[1, :] * berry_spacing
                        )
                        k_points_path[3, :] = (
                            k_point_tmp
                            - normalized_chosen_reciprocal_plane[1, :] * berry_spacing
                        )

                        for s in range(0, 4):
                            (
                                list_magnonic_branches,
                                omega_k_tmp[s],
                                u_k_tmp[s],
                            ) = calculation_omega_and_u_checking_degeneracy(
                                k_points_path[s],
                                threshold_omega,
                                **list_magnonic_branches,
                            )
                                  ##in the not-degenerate case the Berry curvature is straightforwardath

                        if len(list_magnonic_branches) == self.N:
                            phase = np.asarray(
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
                                  ##in the degenerate case the Berry curvature is calculated using the Non-Abelian formulation
                        else:
                            print("degeneracy detected")
                            phase_degenerate = {}
                            for (
                                key,
                                values,
                            ) in list_magnonic_branches.items():
                                phase_degenerate[key] = (
                                    np.asarray(
                                        function_phase_eigenspace_indicated(
                                            u_k_tmp[0], u_k_tmp[1], values
                                        )
                                    )
                                    + np.asarray(
                                        function_phase_eigenspace_indicated(
                                            u_k_tmp[1], u_k_tmp[2], values
                                        )
                                    )
                                    - np.asarray(
                                        function_phase_eigenspace_indicated(
                                            u_k_tmp[3], u_k_tmp[2], values
                                        )
                                    )
                                    - np.asarray(
                                        function_phase_eigenspace_indicated(
                                            u_k_tmp[0], u_k_tmp[3], values
                                        )
                                    )
                                )
                                phase_degenerate[key] = (
                                    phase_degenerate[key]
                                    * k_points_list_tmp[r, 3]
                                )
                                for value in values:
                                    phase[value] = phase_degenerate[key]
                        check_tmp = check_tmp + k_points_list_tmp[r, 3]
                        total_phase = total_phase + phase
                        local_phase = (
                            local_phase + phase * k_points_list_tmp[r, 3]
                        )
                        phases_list.append(phase * k_points_list_tmp[r, 3])
                        omegas_list.append(
                                  list(
                                      map(
                                          lambda x, y, z, t: (x + y + z + t) / 4,
                                          omega_k_tmp[0],
                                          omega_k_tmp[1],
                                          omega_k_tmp[2],
                                          omega_k_tmp[3],
                                      )
                                  )
                              )
                        kpoints_list.append(k_points_list_tmp[r, :3])
                        ##we need at least a point to evaluate if the phase is increasing significantly, if the 00 point has some problems, you can consider a translation in the 2d k grid
                    if i != 0 or j != 0:
                          ##norm of a vector of type 2 (:max element)
                        max_difference_phase = (
                            np.abs(
                                np.asarray(total_phase)
                                - np.asarray(old_total_phase)
                            )
                        ).max()
                        max_phase = (
                            np.abs(np.asarray(old_total_phase))
                        ).max()
                        variation = max_difference_phase / max_phase
                              # print(variation)
                        if variation < phase_threshold:
                            chern_number_tmp[i][j] = local_phase
                            check_tot = check_tot + check_tmp
                            flag_local_evaluation = (
                                number_iterations_dynamical_refinment + 1
                            )
                        else:
                            ##a dynmical refinment is considered if the phase has a percentage increase of more than phase_threshold
                            print(
                                f"dynamical refinment {flag_local_evaluation}/{number_iterations_dynamical_refinment-1}"
                            )
                            k_points_grid_2d[i][j]=dynamical_refinment(
                                k_points_grid_2d[i][j],
                                normalized_chosen_reciprocal_plane,
                                refinment_spacing,
                                refinment_iteration
                                + added_refinment_iteration
                                * flag_local_evaluation,
                                symmetry,
                                threshold_k_grid
                            )
                            flag_local_evaluation = (
                                flag_local_evaluation + 1
                            )
                    else:
                        chern_number_tmp[i][j] = local_phase
                        check_tot = check_tot + check_tmp
                        flag_local_evaluation = (
                            number_iterations_dynamical_refinment + 1
                        )
                print("\n")
                old_total_phase = total_phase
                if (
                    flag_local_evaluation
                    == number_iterations_dynamical_refinment
                ):
                    chern_number_tmp[i][j] = local_phase
                    check_tot = check_tot + check_tmp
                chern_number = chern_number + chern_number_tmp[i][j]
                    ###printing
                k_point_coordinates = np.asarray(
                    k_points_grid_2d[i][j][:,:3]
                )
                chern_number_local = np.asarray(chern_number_tmp[i][j])
                file_1.write(
                    f"{k_point_coordinates} {chern_number_local}   \n"
                )
                kpoints_list = np.reshape(
                  kpoints_list, (len(kpoints_list), 3)
                )
                phases_list = np.reshape(
                        phases_list, (len(phases_list), self.N)
                )
                omegas_list = np.reshape(
                    omegas_list, (len(omegas_list), self.N)
                )
                for s in range(len(kpoints_list)):
                    file_2.write(f"{kpoints_list[s]} {omegas_list[s]}   \n")
                for t in range(len(kpoints_list)):
                    file_3.write(f"{kpoints_list[t]} {phases_list[t]}   \n")
                  ###to check if the weight of the different k points is summing to one
                print("check ", (check_tot / 1), " over 1\n")
     ##         file_3.close()
     ##     file_2.close()
     ## file_1.close()

    ##    ###printing of the new k point list after the dynamical refinment (to check)
    ##    with open("k_points_dynamically_refined.txt", "w") as file:
    ##        for i in range(0, k0):
    ##            for j in range(0, k1):
    ##                file.write(f"{i} {j} {k_points_grid_2d[i][j].refined_grid}\n")
    ##    file.close()
##
        return chern_number

    def __call__(self, *args, **kwargs):
        return self.omegas(*args, **kwargs)

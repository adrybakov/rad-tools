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

from radtools.crystal.kpoints import Kpoints
from radtools.geometry import span_orthonormal_set
from radtools.magnons.diagonalization import ColpaFailed, solve_via_colpa
from radtools.spinham.hamiltonian import SpinHamiltonian

__all__ = ["MagnonDispersion"]


class MagnonDispersion:
    r"""
    Magnon dispersion wrapper.

    Parameters
    ----------
    model : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian.
    Q : (3,) |array_like|_
        Ordering wave vector of the spin-spiral.
        In relative coordinates with respect to the model`s reciprocal cell.
        It rotates spins from unit cell to unit cell,
        but not from atom to atom in (0, 0, 0) unit cell.
    n : (3,) |array_like|_, optional
        Global rotational axis. If None provided, then it is set to the direction of ``Q``.
    nodmi : bool, default=False
        If True, then DMI is not included in the dispersion.
    noaniso : bool, default=False
        If True, then anisotropy is not included in the dispersion.
    custom_mask : func
        Custom mask for the exchange parameter. Function which take (3,3) :numpy:`ndarray`
        as an input and returns (3,3) :numpy:`ndarray` as an output.

    Attributes
    ----------
    Q : (3,) :numpy:`ndarray`
        Ordering wave vector of the spin-spiral. in absolute coordinates in reciprocal space.
    n : (3,) :numpy:`ndarray`
        Global rotational axis.
    N : int
        Number of magnetic atoms.
    J_matrices : (M,) :numpy:`ndarray`
        Exchange parameters.
    i : (M,) :numpy:`ndarray`
        Indices of the first atom in the exchange pair.
    j : (M,) :numpy:`ndarray`
        Indices of the second atom in the exchange pair.
    d : (M, 3) :numpy:`ndarray`
        Vectors from the first atom to the second atom in the exchange pair.
    S : (N, 3) :numpy:`ndarray`
        Spin vectors.
    u : (N, 3) :numpy:`ndarray`
        Defined from local spin directions.
    v : (N, 3) :numpy:`ndarray`
        Defined from local spin directions.
    """

    def __init__(
        self,
        model: SpinHamiltonian,
        Q=None,
        n=None,
        nodmi=False,
        noaniso=False,
        custom_mask=None,
    ):
        self._C = None
        # Store the exchange model, but privately
        self._model = deepcopy(model)
        self._model.notation = "SpinW"

        # Convert Q to absolute coordinates
        if Q is None:
            Q = [0, 0, 0]
        self.Q = np.array(Q, dtype=float) @ self._model.reciprocal_cell

        # Convert n to absolute coordinates, use Q if n is not provided
        if n is None:
            if np.allclose([0, 0, 0], Q):
                self.n = np.array([0, 0, 1])
            else:
                self.n = self.Q / np.linalg.norm(self.Q)
        else:
            self.n = np.array(n, dtype=float) / np.linalg.norm(n)

        # Get the number of magnetic atoms
        self.N = len(self._model.magnetic_atoms)

        # Get the exchange parameters, indices and vectors form the SpinHamiltonian
        (
            self.J_matrices,
            self.indices_i,
            self.indices_j,
            self.dis_vectors,
        ) = self._model.input_for_magnons(
            nodmi=nodmi, noaniso=noaniso, custom_mask=custom_mask
        )

        # Initialize spin vector, u vector and v vector arrays
        self.S = np.zeros((self.N, 3), dtype=float)
        self.u = np.zeros((self.N, 3), dtype=complex)
        self.v = np.zeros((self.N, 3), dtype=complex)

        # Get spin vectors and compute u and v vectors from local spin directions
        for a_i, atom in enumerate(self._model.magnetic_atoms):
            try:
                self.S[a_i] = atom.spin_vector
                e1, e2, e3 = span_orthonormal_set(self.S[a_i])
                self.v[a_i] = e3
                self.u[a_i] = e1 + 1j * e2
            except ValueError:
                raise ValueError(
                    f"Spin vector is not defined for {atom.fullname} atom."
                )

        # Rotate exchange matrices
        for i in range(len(self.J_matrices)):
            rotvec = self.n * (self.Q @ self.dis_vectors[i])
            R_nm = Rotation.from_rotvec(rotvec).as_matrix()
            self.J_matrices[i] = self.J_matrices[i] @ R_nm

    def J(self, k):
        r"""
        Computes J(k) matrix.

        .. math::

            \boldsymbol{J}_{i,j}(\boldsymbol{k}) = \sum_{\boldsymbol{d}}\boldsymbol{J}_{i,j}(\boldsymbol{d})e^{-i\boldsymbol{k}\boldsymbol{d}}

        where indices :math:`i` and :math:`j` correspond to the atoms in the exchange pair.

        Parameters
        ----------
        k : (3,) |array_like|_
            Reciprocal vector.
            In absolute coordinates.
        """
        # Initialize matrix
        result = np.zeros((self.N, self.N, 3, 3), dtype=complex)
        k = np.array(k)
        # Compute J(k)
        for index in range(len(self.J_matrices)):
            i = self.indices_i[index]
            j = self.indices_j[index]

            result[i][j] += self.J_matrices[index] * np.exp(
                -1j * (k @ self.dis_vectors[index])
            )
        return result

    def A(self, k):
        r"""
        Computes A(k) matrix.

        .. math::

            A(\boldsymbol{k})^{ij} = \dfrac{\sqrt{S_i\cdot S_j}}{2}\boldsymbol{u}^T_i\boldsymbol{J}_{i,j}(-\boldsymbol{k})\overline{\boldsymbol{u}}_j

        where indices :math:`i` and :math:`j` correspond to the atoms in the exchange pair.

        Parameters
        ----------
        k : (3,) |array_like|_
            Reciprocal vector.
            In absolute coordinates.
        """
        # Initialize matrix
        result = np.zeros((self.N, self.N), dtype=complex)
        k = np.array(k)

        # Compute A(k)
        J = self.J(-k)
        for i in range(len(J)):
            for j in range(len(J[i])):
                result[i][j] += (
                    sqrt(np.linalg.norm(self.S[i]) * np.linalg.norm(self.S[j]))
                    / 2
                    * (self.u[i] @ J[i][j] @ np.conjugate(self.u[j]))
                )
        return result

    def B(self, k):
        r"""
        Computes B(k) matrix.

        .. math::

            B(\boldsymbol{k})^{ij} = \dfrac{\sqrt{S_i\cdot S_j}}{2}\boldsymbol{u}^T_i\boldsymbol{J}_{i,j}(-\boldsymbol{k})\boldsymbol{u}_j

        where indices :math:`i` and :math:`j` correspond to the atoms in the exchange pair.

        Parameters
        ----------
        k : (3,) |array_like|_
            Reciprocal vector.
            In absolute coordinates.
        """
        # Initialize matrix
        result = np.zeros((self.N, self.N), dtype=complex)
        k = np.array(k)

        # Compute B(k)
        J = self.J(-k)
        for i in range(len(J)):
            for j in range(len(J[i])):
                result[i][j] += (
                    sqrt(np.linalg.norm(self.S[i]) * np.linalg.norm(self.S[j]))
                    / 2
                    * (self.u[i] @ J[i][j] @ self.u[j])
                )
        return result

    def C(self):
        r"""
        Computes C matrix.

        .. math::

            C^{i,j} = \delta_{i,j}\sum_{l}S_l \boldsymbol{v}^T_i\boldsymbol{J}_{i, l}(\boldsymbol{0})\boldsymbol{v}_l

        where indices :math:`i` and :math:`j` correspond to the atoms in the exchange pair.
        """

        if self._C is None:
            self._C = np.zeros((self.N, self.N), dtype=complex)

            # Compute C matrix, note: sum over l is hidden here
            J = self.J(np.zeros(3))
            for i in range(len(J)):
                for j in range(len(J[i])):
                    self._C[i][i] += (
                        np.linalg.norm(self.S[j]) * self.v[i] @ J[i][j] @ self.v[j]
                    )
        return self._C

    def h(self, k):
        k = np.array(k)
        # Compute h matrix
        left = np.concatenate(
            (2 * self.A(k) - 2 * self.C(), 2 * np.conjugate(self.B(k)).T), axis=0
        )
        right = np.concatenate(
            (2 * self.B(k), 2 * np.conjugate(self.A(-k)) - 2 * self.C()), axis=0
        )
        h = np.concatenate((left, right), axis=1)
        return h

    def omega(self, k, flag, zeros_to_none=False):
        r"""
        Computes magnon energies.

        Parameters
        ----------
        k : (3,) |array_like|_
            Reciprocal vector.
            In absolute coordinates.
        zeros_to_none : bool, default=False
            If True, then return ``None`` instead of 0 if Colpa fails.

        Returns
        -------
        omegas : (N,) :numpy:`ndarray`
            Magnon energies for the vector ``k``.
        """
        # Diagonalize h matrix via Colpa
        k = np.array(k)
        h = self.h(k)
        U_tmp = np.zeros((self.N,self.N),dtype=complex)
        try:
            omegas, U = solve_via_colpa(h)
            U_tmp = U
            omegas = omegas.real[: self.N]
        except ColpaFailed:
            # Try to solve for positive semidefinite matrix
            try:
                omegas, U = solve_via_colpa(h + np.diag(1e-8 * np.ones(2 * self.N)))
                U_tmp = U
                omegas = omegas.real[: self.N]
            except ColpaFailed:
                # Try to solve for negative defined matrix
                try:
                    omegas, U = solve_via_colpa(-h)
                    U_tmp = U
                    omegas = omegas.real[: self.N] * -1
                except ColpaFailed:
                    # Try to solve for negative semidefinite matrix
                    try:
                        omegas, U = solve_via_colpa(
                            -h - np.diag(1e-8 * np.ones(h.shape[0]))
                        )
                        U_tmp = U
                        omegas = omegas.real[: self.N] * -1
                    except ColpaFailed:
                        # If all fails, return None or 0
                        if zeros_to_none:
                            omegas = np.array([None] * self.N)
                        else:
                            omegas = np.zeros(self.N)
        omegas[np.abs(omegas) <= 1e-8] = 0
        #considering right number of eigenvectors
        U=U_tmp[:self.N,:self.N]
        #ordering the eigenvalues and the associated eigenvectors
        ordering=omegas.argsort()[::-1]
        omegas=omegas[ordering]
        U=U[:,ordering]
        #U=np.zeros((self.N,self.N))
        #normalizing the eigenvectors
        for i in range(0,len(omegas)):
            norm=np.dot(np.conj(U[i,:]),U[i,:])
            if norm !=0:
                U[i,:]=U[i,:]/norm
        if flag==0:
            return omegas
        else:
            return omegas,U

    def omegas(self, kpoints, zeros_to_none=False):
        r"""
        Dispersion spectra.

        Parameters
        ----------
        kpoints : (M, 3) |array_like|_
            K points in absolute coordinates.
        zeros_to_none : bool, default=False
            If True, then return ``None`` instead of 0 if Colpa fails.
        """

        flag=0
        data = []

        if isinstance(kpoints, Kpoints):
            kpoints = kpoints.points()

        for point in kpoints:
            data.append(self.omega(point, flag, zeros_to_none=zeros_to_none))

        return np.array(data).T
    
    
    def k_points_grid_2d_refinment_and_symmetry(self,plane_2d,grid_spacing,shift_in_plane,shift_in_space,symmetry,refinment_spacing,refinment_iteration,threshold_k_grid):
        r"""
        2d k points grid 
        Parameters
        ----------
        plane_2d : (1,3) |array_like|_  ex. [0,1,1] the plane is the b2xb3 in reciprocal space where the bi are the primitive vectors
        grid_spacing: |double| spacing of the largest k grid
        shift_in_plane : (1,2) |array_like|_  shift of the reciprocal plane_2d with respcet to the reciprocal(crystal) coordinates
        shift_in_space : (1,3) |array_like|_  shift of the reciprocal plane_2d with respcet to the cartesian coordinates
        symmetry : list [symmetry1,symmetry2...] each symmetry is a list of three elements: the versor is associated to the axis of rotation, while the modulus is the angle of rotation
        refinment_spacing : |double| first refinment iteration dimension, the next refinment itarations are multiple of this 
        refinment_iteration : |int| number of refinment iterations
        threshold_k_grid: |double| threshold for the symmetry checking
        
        Return
        -------------
        normalized_chosen_reciprocal_plane: (2,3) |matrix| versor of the chosen plane 
        k0,k1 (1,1) |int,int| the number of k points along the 2 direction for the largest k grid
        total_number_k_points |int|, total number of k points, taking into account also the refinment k points
        k_points_grid_2d k0xk1 grid |k_point class| the largest k grid, where each point is associated to a refinment grid
        """
         
        ##defining a k point class; the user has the possibility to associate to each k point a refinment grid
        ## to each k point is associated the coordinates, the weight of the point, the reciprocal plane in which the k point lives, and the threshold to check the symmetries
        ## moreover, to each k point is associated a refinment grid of k points
        class K_point:
            def __init__(self,coordinates=None,weight=None,reciprocal_vectors_2d=None,symmetry=None,threshold_k_grid=None):
                if coordinates is None:
                    self.coordinates=np.zeros(3)
                else:
                    self.coordinates=coordinates
                if weight is None:
                    self.weight=0
                else:    
                    self.weight=weight
                if symmetry is None:
                    self.symmetry=[]
                else:
                    self.symmetry=symmetry
                self.reciprocal_vectors_2d=reciprocal_vectors_2d
                if threshold_k_grid is None:
                    threshold_k_grid=0.001
                else:
                    self.threshold_k_grid=threshold_k_grid
                self.number_k_points=1
                self.refined_grid=[]
                self.refined_weight=[]
                
            # only rotations are considered 
            # rotation is 4-dimensional: axis of rotation + angle of rotation
            # this is a virtual-method like
            def symmetry_transformation(self,k_origin,k_point,axis):
                rotation=np.sqrt(np.dot(axis,axis))
                #print(k_point-k_origin,axis)
                axis=axis/rotation
                k_vector=np.zeros(3)
                for i in range(0,3):
                    k_vector[i]=k_point[i]-k_origin[i]
                module=np.sqrt(np.dot(k_vector,k_vector))
                unit_k_vector=k_vector/module
                W=np.zeros((3,3))
                I=np.zeros((3,3))
                R=np.zeros((3,3))
                I[0][0]=1
                I[1][1]=1
                I[2][2]=1
                W[0][1]=-axis[2]
                W[0][2]=axis[1]
                W[1][2]=-axis[0]
                W[1][0]=axis[2]
                W[2][0]=-axis[1]
                W[2][1]=axis[0]
                def mul_mat_by_scalar(matrix, scalar):
                    return [[scalar*j for j in i] for i in matrix]
                def mul_mat_by_mat(matrixa,matrixb):
                    product=np.zeros((3,3))
                    for k in range(0,3):
                        for l in range(0,3):
                            for r in range(0,3):
                                product[k][l]=product[k][l]+matrixa[k][r]*matrixb[r][l]
                    return product
                R=I+mul_mat_by_scalar(W,np.sin(rotation))+mul_mat_by_scalar(mul_mat_by_mat(W,W),2*pow(np.sin(rotation/2),2))
                k_point_transformed=np.zeros(3)
                for i in range(0,3):
                    for j in range(0,3):
                        k_point_transformed[i]=k_point_transformed[i]+R[i][j]*unit_k_vector[j]
                k_point_transformed=k_point_transformed*module+k_origin
                return k_point_transformed

            ##the symmetry analisys is applied on a subset of k points (k_points_subgrid), the origin of the subsystem is considered as well (k_origin)
            ##in our case the origin is the point to which the refinment procedure is applied
            ##from the symmetry annalysis is clear the distribution of the origin weight (k_origin_weight) between the different k points of the subset
            def symmetry_analysis(self,k_origin,k_origin_weight,k_points_subgrid,symmetry,threshold_k_grid):
                #print("Inside: ",k_points_subgrid)
                k_points_subgrid_weight_tmp=np.zeros(len(k_points_subgrid[:,0]))
                for i in range(0,len(k_points_subgrid[:,0])):
                    k_points_subgrid_weight_tmp[i]=k_origin_weight/4
                if symmetry[0][0]==symmetry[0][1]==symmetry[0][2]==0:
                    return  k_points_subgrid_weight_tmp
                check_degeneracy=np.zeros((len(k_points_subgrid[:,0]),len(k_points_subgrid[:,0])),dtype=bool)
                degeneracy=0
                for i in range(0,len(k_points_subgrid[:,0])-1):
                    for j in range(i+1,len(k_points_subgrid[:,0])):
                        r=0
                        while r < len(symmetry):
                            #print(symmetry[r])
                            #print(i,j,r,k_points_subgrid[i,:],k_points_subgrid[j,:],symmetry[r])
                            k_point_transformed=self.symmetry_transformation(k_origin,k_points_subgrid[j,:],symmetry[r]) 
                            #print(k_points_subgrid[j,:]-k_origin,k_point_transformed-k_origin)
                            #k_point_transformed=np.zeros(3)
                            flag=False
                            degeneracy_tmp=0
                            for s in range(len(k_point_transformed)):
                                flag=np.isclose(k_points_subgrid[i,s],k_point_transformed[s],atol=threshold_k_grid)
                                if flag==True:
                                    degeneracy_tmp=degeneracy_tmp+1
                            if degeneracy_tmp == len(k_point_transformed):
                                check_degeneracy[i][j]=True
                                degeneracy=degeneracy+1
                                r=len(symmetry)
                            else:
                                r=r+1
                            #print(i,j,r,check_degeneracy,degeneracy)
                if degeneracy==0:
                    for i in range(0,len(k_points_subgrid[:,0])):
                        k_points_subgrid_weight_tmp[i]=k_origin_weight/4
                    return k_points_subgrid_weight_tmp
                else:
                    #print(check_degeneracy)
                    list={}
                    for i in range(0,len(k_points_subgrid[:,0])):
                        list[str(i)]=[i]
                    #print(list)
                    for l in range(0,len(k_points_subgrid[:,0])-1):
                        for j in range(l+1,len(k_points_subgrid[:,0])):
                            #print(list)
                            if len(list) != 1:
                                if(check_degeneracy[l][j]==True):
                                    ### i < j (by definition)
                                    for key,values in list.items():
                                        for value in values:
                                            if value == l:
                                                positionl=key
                                            if value == j:
                                                positionj=key
                                    #print(l,j,positionl,positionj)
                                    list_new={}
                                    if positionl!=positionj:
                                        count=0
                                        for key,values in list.items():
                                            #print(i,j,positioni,positionj,list[key])
                                            if key==positionl:
                                                list_new[str(count)]=values+list[positionj]
                                                count=count+1
                                            else:
                                                if key!=positionj:
                                                    list_new[str(count)]=values
                                                    count=count+1
                                    list={}
                                    for key,values in list_new.items():
                                        list[str(key)]=list_new[str(key)]
                                    #print(i,j,list,list_new)
                    for key,values in list_new.items():
                        list_new[str(key)]=sorted(set(list_new[str(key)]))
                    #print(list,list_new)
                    if len(list_new)!=0:
                        k_weight_tmp=k_origin_weight/len(list_new)
                    else:
                        k_weight_tmp=0.0
                    ##print(list_new,k_weight_tmp)
                    for key,values in list_new.items():
                        if len(values)==1:
                            k_points_subgrid_weight_tmp[values]=k_weight_tmp
                        else:
                            count=0
                            for value in values:
                                if count==0:
                                    k_points_subgrid_weight_tmp[value]=k_weight_tmp
                                    count=count+1
                                else:
                                    k_points_subgrid_weight_tmp[value]=0
                    #print(k_points_subgrid_weight_tmp)
                    #print("END",k_origin_weight,k_points_subgrid_weight_tmp)
                    return k_points_subgrid_weight_tmp
            
            ##this is the local refinment applied to the particular k point (it is a recursive function), at each level of refinment a symmetry analysis is done
            def local_refinment(self,reciprocal_vectors_2d,k_origin,k_origin_weight,refinment_spacing,refinment_iteration,symmetry,threshold_k_grid):
                #print(refinment_iteration,k_origin,k_origin_weight)
                if refinment_iteration==0:
                    k_origin=list(k_origin.flatten())
                    if k_origin_weight != 0:
                        k_origin.append(k_origin_weight)
                        self.refined_grid.extend(k_origin)
                        self.number_k_points=self.number_k_points+len(k_origin)
                else:
                    if k_origin_weight != 0:
                        k_points_subgrid=np.zeros((4,3))
                        k_points_subgrid[0,:]=k_origin+reciprocal_vectors_2d[0,:]*refinment_spacing
                        k_points_subgrid[1,:]=k_origin-reciprocal_vectors_2d[0,:]*refinment_spacing
                        k_points_subgrid[2,:]=k_origin+reciprocal_vectors_2d[1,:]*refinment_spacing
                        k_points_subgrid[3,:]=k_origin-reciprocal_vectors_2d[1,:]*refinment_spacing
                        #print(k_origin,k_points_subgrid,reciprocal_vectors_2d[0,:]*refinment_spacing)
                        k_points_subgrid_weight=self.symmetry_analysis(k_origin,k_origin_weight,k_points_subgrid,symmetry,threshold_k_grid)
                        #test
                        #k_points_subgrid_weight=np.zeros(4)
                        for i in range(0,4):
                            new_refinment_spacing=refinment_spacing/2
                            #k_points_subgrid_weight[i]=k_origin_weight/4
                            self.local_refinment(reciprocal_vectors_2d,k_points_subgrid[i,:],k_points_subgrid_weight[i],new_refinment_spacing,refinment_iteration-1,symmetry,threshold_k_grid)

            ##this is another procedure to add detail to a preceding refinment analysis (it is like continuing the preceding local refinment)
            def dynamical_refinment(self,chosen_reciprocal_plane,refinment_spacing,symmetry,threshold_k_grid,added_refinment_iteration,refinment_iteration):
                new_refinment_spacing=(refinment_spacing/(pow(2,refinment_iteration+added_refinment_iteration)))
                print(f"new refinment spacing {new_refinment_spacing}")
                length_refined_grid=len(self.refined_grid)
                k_point_grid_dynamically_refined=[K_point(coordinates=self.refined_grid[i][:3],reciprocal_vectors_2d=chosen_reciprocal_plane,weight=self.refined_grid[i][3],threshold_k_grid=threshold_k_grid,symmetry=symmetry) for i in range(length_refined_grid)]
                #print("OLD:", len(self.refined_grid),self.refined_grid)
                self.refined_grid=[]
                for j in range(length_refined_grid):
                    #print(chosen_reciprocal_plane,k_point_grid_dynamically_refined[j].coordinates,k_point_grid_dynamically_refined[j].weight,new_refinment_spacing,added_refinment_iteration,symmetry,threshold_k_grid)
                    k_point_grid_dynamically_refined[j].local_refinment(chosen_reciprocal_plane,k_point_grid_dynamically_refined[j].coordinates,k_point_grid_dynamically_refined[j].weight,new_refinment_spacing,added_refinment_iteration,symmetry,threshold_k_grid)
                    self.refined_grid.extend(k_point_grid_dynamically_refined[j].refined_grid)
                refined_grid_number=int(len(self.refined_grid)/4)
                self.refined_grid=np.reshape(self.refined_grid,(refined_grid_number,4))
                #print("NEW: ", len(self.refined_grid),self.refined_grid)

        ###here the largest k grid is built, and at the same time the refinment procedure is applied
        count=0
        chosen_reciprocal_plane=np.zeros((2,3))
        for i in range(0,3):
            if plane_2d[i]!=0 and count<2:
                if i==0:
                    chosen_reciprocal_plane[count]=self._model.b1
                elif i ==1:
                    chosen_reciprocal_plane[count]=self._model.b2
                else:
                    chosen_reciprocal_plane[count]=self._model.b3
                count=count+1
        k0=int(np.dot(chosen_reciprocal_plane[0],chosen_reciprocal_plane[0])/grid_spacing)
        k1=int(np.dot(chosen_reciprocal_plane[1],chosen_reciprocal_plane[1])/grid_spacing)
        weight=1/(k0*k1)
        
        normalized_chosen_reciprocal_plane=np.zeros((2,3))
        for i in range(0,len(chosen_reciprocal_plane[:,0])):
            normalized_chosen_reciprocal_plane[i]=chosen_reciprocal_plane[i]/np.dot(chosen_reciprocal_plane[i],chosen_reciprocal_plane[i])

        k_points_grid_2d=[[K_point(coordinates=None,reciprocal_vectors_2d=chosen_reciprocal_plane,weight=weight,threshold_k_grid=threshold_k_grid,symmetry=symmetry) for s in range(k1)] for r in range(k0)]
        total_number_k_points=0
        with open("k_points.txt",'w') as file:
            for i in range(0,k0):
                for j in range(0,k1):
                    print("k grid: {}% and {}%".format(int(i/k0*100),int(j/k1*100)), end="\r", flush=True)
                    k_points_grid_2d[i][j].coordinates=((i+shift_in_plane[0])/k0)*(shift_in_space+chosen_reciprocal_plane[0])+((j+shift_in_plane[1])/k1)*(shift_in_space+chosen_reciprocal_plane[1])
                    k_points_grid_2d[i][j].local_refinment(normalized_chosen_reciprocal_plane,k_points_grid_2d[i][j].coordinates,weight,refinment_spacing,refinment_iteration,symmetry,threshold_k_grid)
                    refined_grid_tmp=k_points_grid_2d[i][j].refined_grid
                    #print(refined_grid_tmp)
                    refined_grid_number_tmp=int(len(refined_grid_tmp)/4)
                    k_points_grid_2d[i][j].refined_grid=np.reshape(refined_grid_tmp,(refined_grid_number_tmp,4))
                    total_number_k_points=total_number_k_points+k_points_grid_2d[i][j].number_k_points
                    for s in range(len(k_points_grid_2d[i][j].refined_grid[:,0])):
                        file.write(f"{i} {j} {k_points_grid_2d[i][j].coordinates} {k_points_grid_2d[i][j].weight} {k_points_grid_2d[i][j].refined_grid[s]}  \n")
        file.close()
        ##the k points are saved in the file k_points.txt
        
        return normalized_chosen_reciprocal_plane,k0,k1,total_number_k_points,k_points_grid_2d

    def berry_curvature(self,plane_2d,grid_spacing,shift_in_plane,shift_in_space,symmetry,refinment_spacing,refinment_iteration,threshold_k_grid,threshold_omega,added_refinment_iteration,number_iterations_dynamical_refinment):
        r"""
        2berry_curvature
        Parameters
        ----------
        the ones needed in the production of the 2d k grid plut the ones for the berry phase calculation itself
        
        threshold_omega |double| the threshold to check degeneracies in magnonic bands
        added_refinment_iteration |int| number of refinment iteration in case of a phase increase over the phase threshold
        number_iterations_dynamical_refinment |int| number of times the dynamical refinment is applied

        Returns
        ----
        chern numbers double |array_like|_ 
        """
        
        ####Following a few easy functions for easy use
        def function_minimum_difference_pair(a,n):
            a=sorted(a)
            diff=10**20
            for i in range(n-1):
                if a[i+1]-a[i]<diff:
                    diff=a[i+1]-a[i]
            return diff

        # function: scalar product between rows of square matrices
        # a and b are matrices, c is a matrix which rows are the product of the a and b rows
        def function_product_of_rows_eigenvector_matrix(a,b):
            c=np.zeros(len(a[:,0]),dtype=complex)
            for i in range(0,len(a[:,0])):
                for j in range(0,len(a[0,:])):
                    c[i]=c[i]+np.conj(a[i,j])*b[i,j]
            return c
        
        # function: scalar product and phase extraction in a eigenspace with multiplicity greater than zero
        # a and b are matrices, values is a list of rows of a and b
        # the rows are the generator of the degenerate eigenspace
        def function_phase_eigenspace_indicated(a,b,values):
            a_tmp=np.zeros((len(values),len(values)),dtype=complex)
            b_tmp=np.zeros((len(values),len(values)),dtype=complex)
            c=np.zeros((len(values),len(values)),dtype=complex)
            # building submatrices
            for i,eli in enumerate(values):
                for j,elj in enumerate(values):
                    a_tmp[i,j]=a[eli,elj]
                    b_tmp[i,j]=b[eli,elj]
            # multiplying submatrices
            c=np.conj(a)*b
            return np.angle(np.linalg.det(c))
        
        #function: updating information about magnonic branches taking into acount the new/old found degeneracies
        # the old list of magnonic branches is given, togheter with the degenerate eigenvalues and the number of magnetic atoms (or magnonic branches)
        # grouping magnonic branches with respect to degeneracies
        # dictionary: each key is associated to an eigenspace, each value is a list of degenerate branches 
        # list_magnonic_branches.={ 1: [1,2,3], 2: [4,5]..} magnonic branches 1 2 and 3 are degenerate...
        def function_update_list_magnons(check_degeneracy,number_magnons,**list_magnonic_branches):
            #print(list_magnonic_branches)
            if len(list_magnonic_branches)==1:
                return list_magnonic_branches
            
            list_magnonic_branches_tmp={}
            for key,values in list_magnonic_branches.items():
                list_magnonic_branches_tmp[str(key)]=list_magnonic_branches[str(key)]
            
            for i in range(0,number_magnons-1):
                for j in range(i+1,number_magnons):
                    if len(list_magnonic_branches_tmp)==1:
                        break
                    else:
                        flag1=0
                        list_magnonic_branches_new={}
                        if(check_degeneracy[i][j]==True):
                            flag2=0
                            for key,values in list_magnonic_branches_tmp.items():
                                for value in values:
                                    if value == i:
                                        positioni=key
                                        flag2=flag2+1
                                    if value == j:
                                        positionj=key
                                        flag2=flag2+1
                            if flag2==2:
                                if positioni!=positionj:
                                    count=0
                                    for key,values in list_magnonic_branches_tmp.items():
                                        if key==positioni:
                                            list_magnonic_branches_new[str(count)]=values+list_magnonic_branches_tmp[str(positionj)]
                                            count=count+1
                                        else:
                                            if key!=positionj:
                                                list_magnonic_branches_new[str(count)]=values
                                                count=count+1
                                else:
                                    flag1=1
                            else:
                               flag1=1
                        else:
                            flag1=1
                        if flag1!=1:
                            list_magnonic_branches_tmp={}
                            for key,values in list_magnonic_branches_new.items():
                                list_magnonic_branches_tmp[str(key)]=list_magnonic_branches_new[str(key)]
                        else:
                            for key,values in list_magnonic_branches_tmp.items():
                                list_magnonic_branches_new[str(key)]=list_magnonic_branches_tmp[str(key)]
                    
            for key,values in list_magnonic_branches_new.items():
                list_magnonic_branches_new[str(key)]=sorted(set(list_magnonic_branches_new[str(key)]))
            return list_magnonic_branches_new

        ##function calculating mangonic energies and eigenvectors and checking any degeneracy; returning also the proper egenspaces
        def calculation_omega_and_u_checking_degeneracy(k_point,threshold_omega,**list_magnonic_branches):      
            degeneracy=0
            omega_k,u_k=self.omega(k_point,True)   
            check_degeneracy=np.zeros((len(omega_k),len(omega_k)),dtype=bool)
            for i in range(0,len(omega_k)-1):
                for j in range(i+1,len(omega_k)):
                    check_degeneracy[i][j]=False
                    check_degeneracy[i][j]=np.isclose(omega_k[i],omega_k[j],atol=threshold_omega)      
                    if check_degeneracy[i][j]==True:
                        #print(omega_k[i],omega_k[j])
                        degeneracy=degeneracy+1
                        #print(degeneracy)
            if degeneracy!=0:
                list_magnonic_branches=function_update_list_magnons(check_degeneracy,len(omega_k),**list_magnonic_branches)
            else:
                if not list_magnonic_branches:
                    for i in range(0,len(omega_k)):
                        list_magnonic_branches[str(i)]=[i]
            return list_magnonic_branches,omega_k, u_k

        #generating k points grid with refinment as already pointed out
        chosen_reciprocal_plane,k0,k1,total_number_k_points,k_points_grid_2d=self.k_points_grid_2d_refinment_and_symmetry(plane_2d,grid_spacing,shift_in_plane,shift_in_space,symmetry,refinment_spacing,refinment_iteration,threshold_k_grid)
        
        # calculating the berry curvature in each k point directly
        # taking into account possible degeneracies (in the non-abelian formulation: pointwise) TO DISCUSS THE THEORETICAL VALIDITY OF THIS
        k_points_path=np.zeros((4,3))
        u_k_tmp=np.zeros((4,self.N,self.N),dtype=complex)
        omega_k_tmp=np.zeros((4,self.N))       

        chern_number_tmp=np.zeros((k0,k1,self.N),dtype=complex)
        chern_number=0
        phase=np.zeros((self.N),dtype=complex)
        
        def function_generate_indipendent_list(n):
            list_magnonic_branches_tmp={}
            for i in range(0,n):
                list_magnonic_branches_tmp[str(i)]=[i]
            return list_magnonic_branches_tmp
        
        list_magnonic_branches=function_generate_indipendent_list(self.N)
        local_phase=np.zeros((self.N),dtype=complex)
        old_total_phase=np.zeros((self.N),dtype=complex)
        check_tot=0
        phase_threshold=0.5
    
        ##here the berry calculation
        # to each k point of the largest grid (i,j) is associated a refinment list (r)
        # for each point of the refinment list (r) are considered 4 k points around it and on these 4 k points it is calculated the line integral of the berry connection
        # if there is any degeneracy in the subset of 4 k points the line integral is properly formulated (degeneracy detected)
        # the berry curvature (the line integral of the berry connection) is properly summed (k point weight) over the refinment list (r)
        # if the berry curvature increases between two successive k points, a dynamical refinment is considered (to increase the precision around points with high increments in the berry curvature)
        # the number of dynamical refinment is fixed by number_iterations_dynamical_refinment
        # all these data are then saved in the respective files
        with open("berry_averaged_points.txt",'w') as file_1:
            with open("magnonic_surfaces.txt", 'w') as file_2:
                with open("berry_detailed_points.txt",'w') as file_3:
                    for i in range(0,k0):
                        for j in range(0,k1):
                            print("berry large grid",int(i/k0*100),"%  and ",int(j/k1*100),"%")
                            print(f"{i,j}")
                            flag_local_evaluation=0
                            while flag_local_evaluation<number_iterations_dynamical_refinment:
                                k_points_list_tmp=k_points_grid_2d[i][j].refined_grid
                                length_refinment_grid=len(k_points_list_tmp)
                                local_phase=0
                                total_phase=0
                                check_tmp=0
                                phases_list=[]
                                kpoints_list=[]
                                omegas_list=[]
                                for r in range(0,length_refinment_grid):
                                    print("status: {}%".format(int(r/length_refinment_grid*100)), end="\r", flush=True)
                                    list_magnonic_branches=function_generate_indipendent_list(self.N)
                                    k_point_tmp=np.asarray(k_points_list_tmp[r,:3])
                                    
                                    if length_refinment_grid>1:
                                        berry_spacing=function_minimum_difference_pair(k_point_tmp,len(k_point_tmp))/2
                                    else:
                                        berry_spacing=1.0e-8
                                    
                                    k_points_path[0,:]=k_point_tmp+chosen_reciprocal_plane[0,:]*berry_spacing
                                    k_points_path[1,:]=k_point_tmp-chosen_reciprocal_plane[0,:]*berry_spacing
                                    k_points_path[2,:]=k_point_tmp+chosen_reciprocal_plane[1,:]*berry_spacing
                                    k_points_path[3,:]=k_point_tmp-chosen_reciprocal_plane[1,:]*berry_spacing                  
                                    
                                    for s in range(0,4):
                                        list_magnonic_branches,omega_k_tmp[s],u_k_tmp[s]=calculation_omega_and_u_checking_degeneracy(k_points_path[s],threshold_omega,**list_magnonic_branches)
                                    ##in the not-degenerate case the Berry curvature is straightforwardath
                                    
                                    if len(list_magnonic_branches)==self.N:
                                        phase=np.asarray(list(map(lambda x,y,z,w: np.angle(np.asarray(x))+np.angle(np.asarray(y))-np.angle(np.asarray(z))-np.angle(np.asarray(w)), function_product_of_rows_eigenvector_matrix(u_k_tmp[0],u_k_tmp[1]), \
                                            function_product_of_rows_eigenvector_matrix(u_k_tmp[1],u_k_tmp[2]), function_product_of_rows_eigenvector_matrix(u_k_tmp[3],u_k_tmp[2]), function_product_of_rows_eigenvector_matrix(u_k_tmp[3],u_k_tmp[0]))))
                                    ##in the degenerate case the Berry curvature is calculated using the Non-Abelian formulation
                                    else:
                                        print("degeneracy detected")
                                        phase_degenerate={}
                                        for key,values in list_magnonic_branches.items():
                                            phase_degenerate[key]=np.asarray(function_phase_eigenspace_indicated(u_k_tmp[0],u_k_tmp[1],values))+np.asarray(function_phase_eigenspace_indicated(u_k_tmp[1],u_k_tmp[2],values))-np.asarray(\
                                                    function_phase_eigenspace_indicated(u_k_tmp[3],u_k_tmp[2],values))-np.asarray(function_phase_eigenspace_indicated(u_k_tmp[3],u_k_tmp[0],values))
                                            phase_degenerate[key]=phase_degenerate[key]*k_points_list_tmp[r,3]
                                            for value in values:
                                                phase[value]=phase_degenerate[key]
                                    check_tmp=check_tmp+k_points_list_tmp[r,3]
                                    total_phase=total_phase+phase
                                    local_phase=local_phase+phase*k_points_list_tmp[r,3]
                                    phases_list.append(phase*k_points_list_tmp[r,3])
                                    omegas_list.append(list(map(lambda x,y,z,t: (x+y+z+t)/4,omega_k_tmp[0],omega_k_tmp[1],omega_k_tmp[2],omega_k_tmp[3])))
                                    kpoints_list.append(k_points_list_tmp[r,:3])
                                ##we need at least a point to evaluate if the phase is increasing significantly, if the 00 point has some problems, you can consider a translation in the 2d k grid
                                if i!=0 or j!=0:
                                    ##norm of a vector of type 2 (:max element)
                                    max_difference_phase=(np.abs(np.asarray(total_phase)-np.asarray(old_total_phase))).max()
                                    max_phase=(np.abs(np.asarray(old_total_phase))).max()
                                    variation=max_difference_phase/max_phase
                                    #print(variation)
                                    if variation < phase_threshold:
                                        chern_number_tmp[i][j]=local_phase
                                        check_tot=check_tot+check_tmp
                                        flag_local_evaluation=number_iterations_dynamical_refinment+1
                                    else:
                                    ##a dynmical refinment is considered if the phase has a percentage increase of more than phase_threshold
                                        print(f"dynamical refinment {flag_local_evaluation}/{number_iterations_dynamical_refinment-1}")
                                        k_points_grid_2d[i][j].dynamical_refinment(chosen_reciprocal_plane,refinment_spacing,symmetry,threshold_k_grid,added_refinment_iteration,refinment_iteration+added_refinment_iteration*flag_local_evaluation)
                                        flag_local_evaluation=flag_local_evaluation+1
                                else:
                                    chern_number_tmp[i][j]=local_phase
                                    check_tot=check_tot+check_tmp
                                    flag_local_evaluation=number_iterations_dynamical_refinment+1
                            print("\n")
                            old_total_phase=total_phase
                            if flag_local_evaluation==number_iterations_dynamical_refinment:
                                chern_number_tmp[i][j]=local_phase
                                check_tot=check_tot+check_tmp
                            chern_number=chern_number+chern_number_tmp[i][j]
                            ###printing
                            k_point_coordinates=np.asarray(k_points_grid_2d[i][j].coordinates)
                            chern_number_local=np.asarray(chern_number_tmp[i][j])
                            file_1.write(f"{k_point_coordinates} {chern_number_local}   \n")
                            kpoints_list=np.reshape(kpoints_list,(len(kpoints_list),3))
                            phases_list=np.reshape(phases_list,(len(phases_list),self.N))
                            omegas_list=np.reshape(omegas_list,(len(omegas_list),self.N))
                            for s in range(len(kpoints_list)):
                                file_2.write(f"{kpoints_list[s]} {omegas_list[s]}   \n")
                            for t in range(len(kpoints_list)):
                                file_3.write(f"{kpoints_list[t]} {phases_list[t]}   \n")
                            ###to check if the weight of the different k points is summing to one
                            print("check ", (check_tot/1), " over 1\n")
                file_3.close()
            file_2.close()
        file_1.close()
        
        ###printing of the new k point list after the dynamical refinment (to check)
        with open("k_points_dynamically_refined.txt",'w') as file:
            for i in range(0,k0):
                for j in range(0,k1):
                    file.write(f"{i} {j} {k_points_grid_2d[i][j].refined_grid}\n")
        file.close()

        return chern_number   

    def __call__(self, *args, **kwargs):
        return self.omegas(*args, **kwargs)
##
##if __name__ == "__main__":
##    import timeit
##    import os
##    import matplotlib.pyplot as plt
##    import numpy as np
##    from termcolor import cprint
##    from radtools.io.internal import load_template
##    from radtools.io.tb2j import load_tb2j_model
##    from radtools.magnons.dispersion import MagnonDispersion
##    from radtools.decorate.stats import logo
##    from radtools.spinham.constants import TXT_FLAGS
##    from radtools.decorate.array import print_2d_array
##    from radtools.decorate.axes import plot_hlines
##
##    spin={"Cr1":[0,0,1],"Cr2":[0,0,1]}
##    input_filename='/home/marcomarino/TestBerry_curvature/exchange.2.out'
##    spinham = load_tb2j_model(input_filename)
##    for key,values in spin.items():
##        atom_name=key
##        atom=spinham.get_atom(atom_name)
##        atom_spin=list(values)
##        atom.spin_vector=atom_spin
##    ##kp = spinham.kpoints
##    ##fig, ax = plt.subplots()
##    ##
##   # dispersion = MagnonDispersion(spinham,nodmi=True,noaniso=True)
##    dispersion = MagnonDispersion(spinham,nodmi=False,noaniso=False)
##    ##omegas = dispersion(kp)
##    ##ax.set_xticks(kp.coordinates(), kp.labels, fontsize=15)
##    ##ax.set_ylabel("E, meV", fontsize=15)
##    ##for omega in omegas:
##    ##    ax.plot(kp.flatten_points(), omega)
##    ##plt.show()
##    plane_2d=[1,1,0]
##    #threshold for understanding if there is a degeneracy
##    threshold_k_grid=0.0000001
##    shift_in_plane=[0,0]
##    shift_in_space=[0,0,0]
##    symmetry=[[0,0,0]]
##    #symmetry=[[0,0,0]]
##    grid_spacing=0.01
##    refinment_iteration=1
##    refinment_spacing=0.005
##    #chosen_reciprocal_plane,k0,k1,total_number_k_points,k_points_grid_2d=dispersion.k_points_grid_2d_refinment_and_symmetry(plane_2d,grid_spacing,shift_in_plane,shift_in_space,symmetry,refinment_spacing,refinment_iteration,threshold_k_grid)
##    ##print(k0,k1,total_number_k_points)
##    #import json
##    #file_k_points=json.load('/home/marcomarino/rad-tools/k_points.txt')
##    #print(file_k_points)
##    threshold_omega=0.6
##    added_refinment_iteration=1
##    ###how many times the dynamical refinement is tried
##    number_iterations_dynamical_refinment=2
##    chern_number=dispersion.berry_curvature(plane_2d,grid_spacing,shift_in_plane,shift_in_space,symmetry,refinment_spacing,refinment_iteration,threshold_k_grid,threshold_omega,added_refinment_iteration,number_iterations_dynamical_refinment)
##    print(chern_number)
##    #sed 's/[][",'"'"']//g' berry_points.txt | awk 'NF>1{print}' > berry_points.cleaned.txt
##    #sed 's/[][",'"'"']//g' magnonic_surfaces.txt | awk 'NF>1{print}' > magnonic_surfaces.cleaned.txt
##     #sed 's/[][",'"'"']//g' k_points.txt | awk 'NF>1{print}' > k_points.cleaned.txt
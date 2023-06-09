from radtools.crystal.crystal import Crystal
from radtools.crystal.atom import Atom

import numpy as np

from tqdm import tqdm
from math import pi

from multiprocessing import Pool, shared_memory


def compute_chunk(data):
    i, pool_size, buffer, shape, dtype = data
    buffer = shared_memory.SharedMemory(name=buffer)
    array = np.ndarray(shape, dtype=dtype, buffer=buffer.buf)[i * pool_size :]
    end_index = min(len(array), pool_size)
    energy = 0
    for i in range(0, end_index):
        m_i = array[i, 0]
        r_i = array[i, 1]
        r_ij = array[i + 1 :, 1] - r_i
        r_ij_norm = np.linalg.norm(r_ij, axis=1)
        r_ij_norm_5 = (r_ij.T / r_ij_norm**5).T
        first_term = np.sum(
            r_ij_norm_5.T * np.diag(array[i + 1 :, 0] @ r_ij.T),
            axis=1,
        )
        second_term = np.sum((array[i + 1 :, 0].T) / r_ij_norm**3, axis=1)
        energy += m_i @ (-3 * first_term + second_term)
    return energy


def dipole_dipole_energy(crystal: Crystal, na, nb, nc, nproc=1):
    r"""
    Computes dipole-dipole energy.

    Formula:

    .. math::

        E = \sum_{i\nej}\left(3(\vec{m_i} \cdot \ver{r_ij})(\vec{m_j} \cdot \ver{r_ij}) - (\vec{m_i}\cdot\vec{m_j})\right)
    """

    constant = 1.25663706212 * 9.2740100783**2 * 6.241509074 / 1000 / 4 / pi

    translations = (
        np.transpose(np.indices((na, nb, nc)), (1, 2, 3, 0)).reshape((na * nb * nc, 3))
        @ crystal.cell
    )
    n_t = len(translations)

    magnetic_atoms = []
    for atom in crystal:
        try:
            tmp = atom.magmom
            magnetic_atoms.append(atom)
        except ValueError:
            pass

    n_mag = len(magnetic_atoms)
    magnetic_centres = np.zeros((n_t * len(crystal), 2, 3), dtype=float)
    for a_i, atom in enumerate(magnetic_atoms):
        magnetic_centres[a_i * n_t : (a_i + 1) * n_t, 0] = np.tile(
            atom.magmom, (n_t, 1)
        )
        magnetic_centres[a_i * n_t : (a_i + 1) * n_t, 1] = translations + atom.position

    if nproc == 1 or len(magnetic_centres) // nproc == 0:
        energy = 0
        for i in tqdm(range(len(magnetic_centres) - 1)):
            m_i = magnetic_centres[i, 0]
            r_i = magnetic_centres[i, 1]
            r_ij = magnetic_centres[i + 1 :, 1] - r_i
            r_ij_norm = np.linalg.norm(r_ij, axis=1)
            r_ij_norm_5 = (r_ij.T / r_ij_norm**5).T
            first_term = np.sum(
                r_ij_norm_5.T * np.diag(magnetic_centres[i + 1 :, 0] @ r_ij.T),
                axis=1,
            )
            second_term = np.sum(
                (magnetic_centres[i + 1 :, 0].T) / r_ij_norm**3, axis=1
            )
            energy += (m_i @ (-3 * first_term + second_term)) / na / nb / nc / n_mag
    else:
        pool_size = len(magnetic_centres) // nproc + 1
        print(magnetic_centres.shape)
        shm = shared_memory.SharedMemory(create=True, size=magnetic_centres.nbytes)
        b = np.ndarray(
            magnetic_centres.shape, dtype=magnetic_centres.dtype, buffer=shm.buf
        )
        b[:] = magnetic_centres[:]
        print(b.shape, b.dtype)
        buffer_name = shm.name

        with Pool(nproc) as p:
            energy = (
                np.sum(
                    list(
                        p.map(
                            compute_chunk,
                            [
                                (i, pool_size, buffer_name, b.shape, b.dtype)
                                for i in range(nproc)
                            ],
                        )
                    )
                )
                / na
                / nb
                / nc
                / n_mag
            )

    return energy * constant


if __name__ == "__main__":
    from radtools.crystal.atom import Atom

    crystal = Crystal()
    crystal.cell = [
        [3.513012, 0.000000000, 0.000000000],
        [0.000000000, 4.752699, 0.000000000],
        [0.000000000, 0.000000000, 23.5428674],
    ]

    Cr1 = Atom(
        "Cr1",
        np.array((0.2500000000, 0.7500000000, 0.1616620796)) @ crystal.cell,
    )
    Cr2 = Atom(
        "Cr2",
        np.array((0.7500000000, 0.2500000000, 0.0737774367)) @ crystal.cell,
    )

    print(Cr1.position, Cr2.position, sep="\n")

    crystal.add_atom(Cr1)
    crystal.add_atom(Cr2)

    crystal.Cr1.magmom = [0, 0, 3]
    crystal.Cr2.magmom = [0, 0, 3]
    from time import time

    start_time = time()
    z = dipole_dipole_energy(crystal, 20, 20, 1)
    print(z, time() - start_time, sep="\n")
    z = dipole_dipole_energy(crystal, 20, 20, 1, 10)
    print(z, time() - start_time, sep="\n")
    # z = dipole_dipole_energy(crystal, 10, 10, 1)
    # crystal.Cr1.magmom = [0, 3, 0]
    # crystal.Cr2.magmom = [0, 3, 0]
    # y = dipole_dipole_energy(crystal, 10, 10, 1)
    # crystal.Cr1.magmom = [3, 0, 0]
    # crystal.Cr2.magmom = [3, 0, 0]
    # x = dipole_dipole_energy(crystal, 10, 10, 1)

    # print("x = ", x)
    # print("y = ", y)
    # print("z = ", z)
    # print("x - y = ", x - y)
    # print("z - y = ", z - y)
    # print("rel = ", (x - y) / (z - y))

    # print("Previous:")
    # print("x = -0.026982064971")
    # print("y = -0.016608788855")
    # print("z = 0.043590853826")
    # print("x - y = -0.010373276115")
    # print("z - y = 0.060199642681")
    # print("rel = ", -0.010373276115 / 0.060199642681)
    # print("Previous 10:")
    # print("x = -0.020418832737")
    # print("y = -0.012159728056")
    # print("z = 0.032578560794")
    # print("x - y = -0.008259104681")
    # print("z - y = 0.04473828885")
    # print("rel = ", -0.008259104681 / 0.04473828885)

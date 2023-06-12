import numpy as np

from tqdm import tqdm
from math import pi

__all__ = ["dipole_dipole_energy", "dipole_dipole_interaction"]

CONSTANT = 1.25663706212 * 9.2740100783**2 * 6.241509074 / 1000 / 4 / pi


def dipole_dipole_energy(magnetic_centres, progress_bar=True, normalize=True):
    r"""
    Computes magnetic dipole-dipole energy.

    This function computes the magnetic dipole-dipole energy of the set of magnetic centres:

    .. math::

        E = \frac{\mu_0}{4\pi}\sum_{i > j}\frac{1}{\vert r_{ij}\vert^3}\left(3(\vec{m_i} \cdot \vec{r_{ij}})(\vec{m_j} \cdot \vec{r_{ij}}) - (\vec{m_i}\cdot\vec{m_j})\right)

    Parameters
    ----------
    magnetic_centres: (N, 2, 3) |array_like|_
        List of N magnetic centres.
        First element along second axis is magnetic moment (in Bohr magnetons),
        Second element along second axis if position (in Angstrom).
    progress_bar : bool, default True
        Whether to show progressbar.
    normalize : bool, default True
        Whether to normalize energy to the number of magnetic centres.

    Returns
    -------
    energy : float
        Dipole-dipole energy of the system of ``magnetic_centres``.
        Normalized to the number of magnetic centres N.

    See Also
    --------
    dipole_dipole_interaction : between two sets of magnetic centres.
    """

    energy = 0
    if progress_bar:
        for_range = tqdm(range(len(magnetic_centres) - 1))
    else:
        for_range = range(len(magnetic_centres) - 1)

    N = len(magnetic_centres)

    for i in for_range:
        m_i = magnetic_centres[i, 0]
        r_i = magnetic_centres[i, 1]
        r_ij = magnetic_centres[i + 1 :, 1] - r_i
        r_ij_norm = np.linalg.norm(r_ij, axis=1)
        r_ij_norm_5 = (r_ij.T / r_ij_norm**5).T
        first_term = np.sum(
            r_ij_norm_5.T * np.diag(magnetic_centres[i + 1 :, 0] @ r_ij.T),
            axis=1,
        )
        second_term = np.sum((magnetic_centres[i + 1 :, 0].T) / r_ij_norm**3, axis=1)
        energy += m_i @ (-3 * first_term + second_term)

    if normalize:
        energy /= N

    return energy * CONSTANT


def dipole_dipole_interaction(
    magnetic_centres1, magnetic_centres2, progress_bar=True, normalize=True
):
    r"""
    Computes magnetic dipole-dipole interaction.

    This function computes dipole-dipole interaction between two sets of magnetic
    centres (:math:`ms_1` and :math:`ms_2`):

    .. math::

        E = \frac{\mu_0}{4\pi}\sum_{i \in mc_1, j \in mc_2}\frac{1}{\vert r_{ij}\vert^3}\left(3(\vec{m_i} \cdot \vec{r_{ij}})(\vec{m_j} \cdot \vec{r_{ij}}) - (\vec{m_i}\cdot\vec{m_j})\right)

    Parameters
    ----------
    magnetic_centres1: (N, 2, 3) |array_like|_
        List of N magnetic centres of the first set.
        First element along second axis is magnetic moment (in Bohr magnetons),
        Second element along second axis if position (in Angstrom).
    magnetic_centres2: (M, 2, 3) |array_like|_
        List of M magnetic centres of the second set.
        First element along second axis is magnetic moment (in Bohr magnetons),
        Second element along second axis if position (in Angstrom).
        Python cycle is running along this set.
    progress_bar : bool, default True
        Whether to show progressbar.
    normalize : bool, default True
        Whether to normalize energy to the number of magnetic centres.

    Returns
    -------
    energy : float
        Magnetic dipole-dipole interaction between ``magnetic_centres1`` and ``magnetic_centres2``.
        Normalized to :math:`N\cdot M`.

    See Also
    --------
    dipole_dipole_energy : within one set of magnetic centres.
    """

    energy = 0
    if progress_bar:
        for_range = tqdm(range(len(magnetic_centres2) - 1))
    else:
        for_range = range(len(magnetic_centres2) - 1)

    NM = len(magnetic_centres1) + len(magnetic_centres2)

    for i in for_range:
        m_i = magnetic_centres2[i, 0]
        r_i = magnetic_centres2[i, 1]
        r_ij = magnetic_centres1[:, 1] - r_i
        r_ij_norm = np.linalg.norm(r_ij, axis=1)
        r_ij_norm_5 = (r_ij.T / r_ij_norm**5).T
        first_term = np.sum(
            r_ij_norm_5.T * np.diag(magnetic_centres1[:, 0] @ r_ij.T),
            axis=1,
        )
        second_term = np.sum((magnetic_centres1[:, 0].T) / r_ij_norm**3, axis=1)
        energy += m_i @ (-3 * first_term + second_term)

    if normalize:
        energy /= NM

    return energy * CONSTANT


if __name__ == "__main__":
    from radtools.crystal.atom import Atom
    from radtools.crystal.crystal import Crystal

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

    crystal.add_atom(Cr1)
    crystal.add_atom(Cr2)

    crystal.Cr1.magmom = [0, 0, 3]
    crystal.Cr2.magmom = [0, 0, 3]

    z10 = crystal.mag_dipdip_energy(10, 10, 1)
    z15 = crystal.mag_dipdip_energy(15, 15, 1)
    energies = crystal.converge_mag_dipdip_energy((10, 10, 1), (5, 5, 0), eps=0.1)
    for e in energies:
        print(f"{e[0]} ({e[1]})<------")

    print(z10, z15, sep="\n")
    # crystal.Cr1.magmom = [0, 3, 0]
    # crystal.Cr2.magmom = [0, 3, 0]
    # y = crystal.mag_dipdip_energy(10, 10, 1)
    # crystal.Cr1.magmom = [3, 0, 0]
    # crystal.Cr2.magmom = [3, 0, 0]
    # x = crystal.mag_dipdip_energy(10, 10, 1)

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

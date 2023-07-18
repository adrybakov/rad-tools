from math import pi

import numpy as np
from tqdm import tqdm

__all__ = ["dipole_dipole_energy", "dipole_dipole_interaction"]

CONSTANT = 1.25663706212 * 9.2740100783**2 * 6.241509074 / 1000 / 4 / pi


def dipole_dipole_energy(magnetic_centres, progress_bar=True, normalize=True):
    r"""
    Computes magnetic dipole-dipole energy.

    This function computes the magnetic dipole-dipole energy of the set of magnetic centres:

    .. math::

        E = -\frac{\mu_0}{4\pi}\sum_{i > j}\frac{1}{\vert r_{ij}\vert^3}\left(3(\vec{m_i} \cdot \vec{r_{ij}})(\vec{m_j} \cdot \vec{r_{ij}}) - (\vec{m_i}\cdot\vec{m_j})\right)

    Parameters
    ----------
    magnetic_centres: (N, 2, 3) |array_like|_
        List of N magnetic centres.
        First element along second axis is magnetic moment (in Bohr magnetons).
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

    magnetic_centres = np.array(magnetic_centres)

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
    centres (:math:`mc_1` and :math:`mc_2`):

    .. math::

        E = -\frac{\mu_0}{4\pi}\sum_{i \in mc_1, j \in mc_2}\frac{1}{\vert r_{ij}\vert^3}\left(3(\vec{m_i} \cdot \vec{r_{ij}})(\vec{m_j} \cdot \vec{r_{ij}}) - (\vec{m_i}\cdot\vec{m_j})\right)

    Parameters
    ----------
    magnetic_centres1: (N, 2, 3) |array_like|_
        List of N magnetic centres of the first set.
        First element along second axis is magnetic moment (in Bohr magnetons).
        Second element along second axis if position (in Angstrom).
    magnetic_centres2: (M, 2, 3) |array_like|_
        List of M magnetic centres of the second set.
        First element along second axis is magnetic moment (in Bohr magnetons).
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

    magnetic_centres1 = np.array(magnetic_centres1)
    magnetic_centres2 = np.array(magnetic_centres2)
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

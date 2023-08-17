from math import cos, sin, tan

import numpy as np

from radtools.constants import TORADIANS

__all__ = [
    "CUB_hs_points",
    "FCC_hs_points",
    "BCC_hs_points",
    "TET_hs_points",
    "BCT_hs_points",
    "ORC_hs_points",
    "ORCF_hs_points",
    "ORCI_hs_points",
    "ORCC_hs_points",
    "HEX_hs_points",
    "RHL_hs_points",
    "MCL_hs_points",
    "MCLC_hs_points",
    "TRI_hs_points",
]


def CUB_hs_points():
    r"""
    Get high-symmetry points for the CUB lattice.

    See :ref:`guide_cub` for the details.

    Returns
    -------
    kpoints : dict
        High-symmetry points.
    """

    return {
        "G": np.array([0, 0, 0]),
        "M": np.array([1 / 2, 1 / 2, 0]),
        "R": np.array([1 / 2, 1 / 2, 1 / 2]),
        "X": np.array([0, 1 / 2, 0]),
    }


def FCC_hs_points():
    r"""
    Get high-symmetry points for the FCC lattice.

    See :ref:`guide_fcc` for the details.

    Returns
    -------
    kpoints : dict
        High-symmetry points.
    """

    return {
        "G": np.array([0, 0, 0]),
        "K": np.array([3 / 8, 3 / 8, 3 / 4]),
        "L": np.array([1 / 2, 1 / 2, 1 / 2]),
        "U": np.array([5 / 8, 1 / 4, 5 / 8]),
        "W": np.array([1 / 2, 1 / 4, 3 / 4]),
        "X": np.array([1 / 2, 0, 1 / 2]),
    }


def BCC_hs_points():
    r"""
    Get high-symmetry points for the CUB lattice.

    See :ref:`guide_bcc` for the details.

    Returns
    -------
    kpoints : dict
        High-symmetry points.
    """

    return {
        "G": np.array([0, 0, 0]),
        "H": np.array([1 / 2, -1 / 2, 1 / 2]),
        "P": np.array([1 / 4, 1 / 4, 1 / 4]),
        "N": np.array([0, 0, 1 / 2]),
    }


def TET_hs_points():
    r"""
    Get high-symmetry points for the TET lattice.

    See :ref:`guide_tet` for the details.

    Returns
    -------
    kpoints : dict
        High-symmetry points.
    """
    return {
        "G": np.array([0, 0, 0]),
        "A": np.array([1 / 2, 1 / 2, 1 / 2]),
        "M": np.array([1 / 2, 1 / 2, 0]),
        "R": np.array([0, 1 / 2, 1 / 2]),
        "X": np.array([0, 1 / 2, 0]),
        "Z": np.array([0, 0, 1 / 2]),
    }


def BCT_hs_points(variation, conv_a, conv_c):
    r"""
    Get high-symmetry points for the BCT lattice.

    See :ref:`guide_bct` for the details.

    Parameters
    ----------
    variation : str
        BCT variation.
    conv_a : float
        Length of the lattice vector of the conventional cell.
    conv_c : float
        Length of the lattice vector of the conventional cell.

    Returns
    -------
    kpoints : dict
        High-symmetry points.
    """

    if variation == "BCT1":
        eta = (1 + conv_c**2 / conv_a**2) / 4
        kpoints = {
            "G": np.array([0, 0, 0]),
            "M": np.array([-1 / 2, 1 / 2, 1 / 2]),
            "N": np.array([0, 1 / 2, 0]),
            "P": np.array([1 / 4, 1 / 4, 1 / 4]),
            "X": np.array([0, 0, 1 / 2]),
            "Z": np.array([eta, eta, -eta]),
            "Z1": np.array([-eta, 1 - eta, eta]),
        }

    elif variation == "BCT2":
        eta = (1 + conv_a**2 / conv_c**2) / 4
        zeta = conv_a**2 / (2 * conv_c**2)
        kpoints = {
            "G": np.array([0, 0, 0]),
            "N": np.array([0, 1 / 2, 0]),
            "P": np.array([1 / 4, 1 / 4, 1 / 4]),
            "S": np.array([-eta, eta, eta]),
            "S1": np.array([eta, 1 - eta, -eta]),
            "X": np.array([0, 0, 1 / 2]),
            "Y": np.array([-zeta, zeta, 1 / 2]),
            "Y1": np.array([1 / 2, 1 / 2, -zeta]),
            "Z": np.array([1 / 2, 1 / 2, -1 / 2]),
        }
    return kpoints


def ORC_hs_points():
    r"""
    Get high-symmetry points for the ORC lattice.

    See :ref:`guide_orc` for the details.

    Returns
    -------
    kpoints : dict
        High-symmetry points.
    """
    return {
        "G": np.array([0, 0, 0]),
        "R": np.array([1 / 2, 1 / 2, 1 / 2]),
        "S": np.array([1 / 2, 1 / 2, 0]),
        "T": np.array([0, 1 / 2, 1 / 2]),
        "U": np.array([1 / 2, 0, 1 / 2]),
        "X": np.array([1 / 2, 0, 0]),
        "Y": np.array([0, 1 / 2, 0]),
        "Z": np.array([0, 0, 1 / 2]),
    }


def ORCF_hs_points(variation, conv_a, conv_b, conv_c):
    r"""
    Get high-symmetry points for the ORCF lattice.

    See :ref:`guide_orcf` for the details.

    Parameters
    ----------
    variation : str
        ORCF variation.
    conv_a : float
        Length of the lattice vector of the conventional cell.
    conv_b : float
        Length of the lattice vector of the conventional cell.
    conv_c : float
        Length of the lattice vector of the conventional cell.

    Returns
    -------
    kpoints : dict
        High-symmetry points.
    """

    if variation == "ORCF1":
        eta = (1 + conv_a**2 / conv_b**2 + conv_a**2 / conv_c**2) / 4
        zeta = (1 + conv_a**2 / conv_b**2 - conv_a**2 / conv_c**2) / 4

        kpoints = {
            "G": np.array([0, 0, 0]),
            "A": np.array([1 / 2, 1 / 2 + zeta, zeta]),
            "A1": np.array([1 / 2, 1 / 2 - zeta, 1 - zeta]),
            "L": np.array([1 / 2, 1 / 2, 1 / 2]),
            "T": np.array([1, 1 / 2, 1 / 2]),
            "X": np.array([0, eta, eta]),
            "X1": np.array([1, 1 - eta, 1 - eta]),
            "Y": np.array([1 / 2, 0, 1 / 2]),
            "Z": np.array([1 / 2, 1 / 2, 0]),
        }
    elif variation == "ORCF2":
        eta = (1 + conv_a**2 / conv_b**2 - conv_a**2 / conv_c**2) / 4
        delta = (1 + conv_b**2 / conv_a**2 - conv_b**2 / conv_c**2) / 4
        phi = (1 + conv_c**2 / conv_b**2 - conv_c**2 / conv_a**2) / 4

        kpoints = {
            "G": np.array([0, 0, 0]),
            "C": np.array([1 / 2, 1 / 2 - eta, 1 - eta]),
            "C1": np.array([1 / 2, 1 / 2 + eta, eta]),
            "D": np.array([1 / 2 - delta, 1 / 2, 1 - delta]),
            "D1": np.array([1 / 2 + delta, 1 / 2, delta]),
            "L": np.array([1 / 2, 1 / 2, 1 / 2]),
            "H": np.array([1 - phi, 1 / 2 - phi, 1 / 2]),
            "H1": np.array([phi, 1 / 2 + phi, 1 / 2]),
            "X": np.array([0, 1 / 2, 1 / 2]),
            "Y": np.array([1 / 2, 0, 1 / 2]),
            "Z": np.array([1 / 2, 1 / 2, 0]),
        }

    elif variation == "ORCF3":
        eta = (1 + conv_a**2 / conv_b**2 + conv_a**2 / conv_c**2) / 4
        zeta = (1 + conv_a**2 / conv_b**2 - conv_a**2 / conv_c**2) / 4

        kpoints = {
            "G": np.array([0, 0, 0]),
            "A": np.array([1 / 2, 1 / 2 + zeta, zeta]),
            "A1": np.array([1 / 2, 1 / 2 - zeta, 1 - zeta]),
            "L": np.array([1 / 2, 1 / 2, 1 / 2]),
            "T": np.array([1, 1 / 2, 1 / 2]),
            "X": np.array([0, eta, eta]),
            "Y": np.array([1 / 2, 0, 1 / 2]),
            "Z": np.array([1 / 2, 1 / 2, 0]),
        }
    return kpoints


def ORCI_hs_points(conv_a, conv_b, conv_c):
    r"""
    Get high-symmetry points for the ORCI lattice.

    See :ref:`guide_orci` for the details.

    Parameters
    ----------
    conv_a : float
        Length of the lattice vector of the conventional cell.
    conv_b : float
        Length of the lattice vector of the conventional cell.
    conv_c : float
        Length of the lattice vector of the conventional cell.

    Returns
    -------
    kpoints : dict
        High-symmetry points.
    """

    zeta = (1 + conv_a**2 / conv_c**2) / 4
    eta = (1 + conv_b**2 / conv_c**2) / 4
    delta = (conv_b**2 - conv_a**2) / (4 * conv_c**2)
    mu = (conv_a**2 + conv_b**2) / (4 * conv_c**2)

    return {
        "G": np.array([0, 0, 0]),
        "L": np.array([-mu, mu, 1 / 2 - delta]),
        "L1": np.array([mu, -mu, 1 / 2 + delta]),
        "L2": np.array([1 / 2 - delta, 1 / 2 + delta, -mu]),
        "R": np.array([0, 1 / 2, 0]),
        "S": np.array([1 / 2, 0, 0]),
        "T": np.array([0, 0, 1 / 2]),
        "W": np.array([1 / 4, 1 / 4, 1 / 4]),
        "X": np.array([-zeta, zeta, zeta]),
        "X1": np.array([zeta, 1 - zeta, -zeta]),
        "Y": np.array([eta, -eta, eta]),
        "Y1": np.array([1 - eta, eta, -eta]),
        "Z": np.array([1 / 2, 1 / 2, -1 / 2]),
    }


def ORCC_hs_points(conv_a, conv_b):
    r"""
    Get high-symmetry points for the ORCC lattice.

    See :ref:`guide_orcc` for the details.

    Parameters
    ----------
    conv_a : float
        Length of the lattice vector of the conventional cell.
    conv_b : float
        Length of the lattice vector of the conventional cell.

    Returns
    -------
    kpoints : dict
        High-symmetry points.
    """

    zeta = (1 + conv_a**2 / conv_b**2) / 4

    return {
        "G": np.array([0, 0, 0]),
        "A": np.array([zeta, zeta, 1 / 2]),
        "A1": np.array([-zeta, 1 - zeta, 1 / 2]),
        "R": np.array([0, 1 / 2, 1 / 2]),
        "S": np.array([0, 1 / 2, 0]),
        "T": np.array([-1 / 2, 1 / 2, 1 / 2]),
        "X": np.array([zeta, zeta, 0]),
        "X1": np.array([-zeta, 1 - zeta, 0]),
        "Y": np.array([-1 / 2, 1 / 2, 0]),
        "Z": np.array([0, 0, 1 / 2]),
    }


def HEX_hs_points():
    r"""
    Get high-symmetry points for the HEX lattice.

    See :ref:`guide_hex` for the details.

    Returns
    -------
    kpoints : dict
        High-symmetry points.
    """

    return {
        "G": np.array([0, 0, 0]),
        "A": np.array([0, 0, 1 / 2]),
        "H": np.array([1 / 3, 1 / 3, 1 / 2]),
        "K": np.array([1 / 3, 1 / 3, 0]),
        "L": np.array([1 / 2, 0, 1 / 2]),
        "M": np.array([1 / 2, 0, 0]),
    }


def RHL_hs_points(variation, conv_alpha):
    r"""
    Get high-symmetry points for the RHL lattice.

    See :ref:`guide_rhl` for the details.

    Parameters
    ----------
    variation : str
        RHL variation.
    alpha : float
        Angle between the lattice vectors.

    Returns
    -------
    kpoints : dict
        High-symmetry points.
    """
    conv_alpha *= TORADIANS

    if variation == "RHL1":
        eta = (1 + 4 * cos(conv_alpha)) / (2 + 4 * cos(conv_alpha))
        nu = 3 / 4 - eta / 2

        return {
            "G": np.array([0, 0, 0]),
            "B": np.array([eta, 1 / 2, 1 - eta]),
            "B1": np.array([1 / 2, 1 - eta, eta - 1]),
            "F": np.array([1 / 2, 1 / 2, 0]),
            "L": np.array([1 / 2, 0, 0]),
            "L1": np.array([0, 0, -1 / 2]),
            "P": np.array([eta, nu, nu]),
            "P1": np.array([1 - nu, 1 - nu, 1 - eta]),
            "P2": np.array([nu, nu, eta - 1]),
            "Q": np.array([1 - nu, nu, 0]),
            "X": np.array([nu, 0, -nu]),
            "Z": np.array([1 / 2, 1 / 2, 1 / 2]),
        }

    elif variation == "RHL2":
        eta = 1 / (2 * tan(conv_alpha / 2) ** 2)
        nu = 3 / 4 - eta / 2

        return {
            "G": np.array([0, 0, 0]),
            "F": np.array([1 / 2, -1 / 2, 0]),
            "L": np.array([1 / 2, 0, 0]),
            "P": np.array([1 - nu, -nu, 1 - nu]),
            "P1": np.array([nu, nu - 1, nu - 1]),
            "Q": np.array([eta, eta, eta]),
            "Q1": np.array([1 - eta, -eta, -eta]),
            "Z": np.array([1 / 2, -1 / 2, 1 / 2]),
        }


def MCL_hs_points(conv_b, conv_c, conv_alpha):
    r"""
    Get high-symmetry points for the MCL lattice.

    See :ref:`guide_mcl` for the details.

    Parameters
    ----------
    conv_b : float
        Length of the lattice vector of the conventional cell.
    conv_c : float
        Length of the lattice vector of the conventional cell.
    conv_alpha : float
        Angle between the lattice vectors.

    Returns
    -------
    kpoints : dict
        High-symmetry points.
    """
    conv_alpha *= TORADIANS

    eta = (1 - conv_b * cos(conv_alpha) / conv_c) / (2 * sin(conv_alpha) ** 2)
    nu = 1 / 2 - eta * conv_c * cos(conv_alpha) / conv_b

    return {
        "G": np.array([0, 0, 0]),
        "A": np.array([1 / 2, 1 / 2, 0]),
        "C": np.array([0, 1 / 2, 1 / 2]),
        "D": np.array([1 / 2, 0, 1 / 2]),
        "D1": np.array([1 / 2, 0, -1 / 2]),
        "E": np.array([1 / 2, 1 / 2, 1 / 2]),
        "H": np.array([0, eta, 1 - nu]),
        "H1": np.array([0, 1 - eta, nu]),
        "H2": np.array([0, eta, -nu]),
        "M": np.array([1 / 2, eta, 1 - nu]),
        "M1": np.array([1 / 2, 1 - eta, nu]),
        "M2": np.array([1 / 2, eta, -nu]),
        "X": np.array([0, 1 / 2, 0]),
        "Y": np.array([0, 0, 1 / 2]),
        "Y1": np.array([0, 0, -1 / 2]),
        "Z": np.array([1 / 2, 0, 0]),
    }


def MCLC_hs_points(variation, conv_a, conv_b, conv_c, conv_alpha):
    r"""
    Get high-symmetry points for the MCLC lattice.

    See :ref:`guide_mclc` for the details.

    Parameters
    ----------
    variation : str
        MCLC variation.
    conv_a : float
        Length of the lattice vector of the conventional cell.
    conv_b : float
        Length of the lattice vector of the conventional cell.
    conv_c : float
        Length of the lattice vector of the conventional cell.
    conv_alpha : float
        Angle between the lattice vectors.

    Returns
    -------
    kpoints : dict
        High-symmetry points.
    """
    conv_alpha *= TORADIANS
    # Parameters
    if variation in ["MCLC1", "MCLC2"]:
        zeta = (2 - conv_b * cos(conv_alpha) / conv_c) / (4 * sin(conv_alpha) ** 2)
        eta = 1 / 2 + 2 * zeta * conv_c * cos(conv_alpha) / conv_b
        psi = 3 / 4 - conv_a**2 / (4 * conv_b**2 * sin(conv_alpha) ** 2)
        phi = psi + (3 / 4 - psi) * conv_b * cos(conv_alpha) / conv_c
    elif variation in ["MCLC3", "MCLC4"]:
        mu = (1 + conv_b**2 / conv_a**2) / 4
        delta = conv_b * conv_c * cos(conv_alpha) / (2 * conv_a**2)
        zeta = (
            mu
            - 1 / 4
            + (1 - conv_b * cos(conv_alpha) / conv_c) / (4 * sin(conv_alpha) ** 2)
        )
        eta = 1 / 2 + 2 * zeta * conv_c * cos(conv_alpha) / conv_b
        phi = 1 + zeta - 2 * mu
        psi = eta - 2 * delta
    elif variation == "MCLC5":
        zeta = (
            conv_b**2 / conv_a**2
            + (1 - conv_b * cos(conv_alpha) / conv_c) / sin(conv_alpha) ** 2
        ) / 4
        eta = 1 / 2 + 2 * zeta * conv_c * cos(conv_alpha) / conv_b
        mu = (
            eta / 2
            + conv_b**2 / (4 * conv_a**2)
            - conv_b * conv_c * cos(conv_alpha) / (2 * conv_a**2)
        )
        nu = 2 * mu - zeta
        rho = 1 - zeta * conv_a**2 / conv_b**2
        omega = (
            (4 * nu - 1 - conv_b**2 * sin(conv_alpha) ** 2 / conv_a**2)
            * conv_c
            / (2 * conv_b * cos(conv_alpha))
        )
        delta = zeta * conv_c * cos(conv_alpha) / conv_b + omega / 2 - 1 / 4

    # Path
    if variation == "MCLC1":
        return {
            "G": np.array([0, 0, 0]),
            "N": np.array([1 / 2, 0, 0]),
            "N1": np.array([0, -1 / 2, 0]),
            "F": np.array([1 - zeta, 1 - zeta, 1 - eta]),
            "F1": np.array([zeta, zeta, eta]),
            "F2": np.array([-zeta, -zeta, 1 - eta]),
            "I": np.array([phi, 1 - phi, 1 / 2]),
            "I1": np.array([1 - phi, phi - 1, 1 / 2]),
            "L": np.array([1 / 2, 1 / 2, 1 / 2]),
            "M": np.array([1 / 2, 0, 1 / 2]),
            "X": np.array([1 - psi, psi - 1, 0]),
            "X1": np.array([psi, 1 - psi, 0]),
            "X2": np.array([psi - 1, -psi, 0]),
            "Y": np.array([1 / 2, 1 / 2, 0]),
            "Y1": np.array([-1 / 2, -1 / 2, 0]),
            "Z": np.array([0, 0, 1 / 2]),
        }
    elif variation == "MCLC2":
        return {
            "G": np.array([0, 0, 0]),
            "N": np.array([1 / 2, 0, 0]),
            "N1": np.array([0, -1 / 2, 0]),
            "F": np.array([1 - zeta, 1 - zeta, 1 - eta]),
            "F1": np.array([zeta, zeta, eta]),
            "F2": np.array([-zeta, -zeta, 1 - eta]),
            "F3": np.array([1 - zeta, -zeta, 1 - eta]),
            "I": np.array([phi, 1 - phi, 1 / 2]),
            "I1": np.array([1 - phi, phi - 1, 1 / 2]),
            "L": np.array([1 / 2, 1 / 2, 1 / 2]),
            "M": np.array([1 / 2, 0, 1 / 2]),
            "X": np.array([1 - psi, psi - 1, 0]),
            "Y": np.array([1 / 2, 1 / 2, 0]),
            "Y1": np.array([-1 / 2, -1 / 2, 0]),
            "Z": np.array([0, 0, 1 / 2]),
        }
    elif variation == "MCLC3":
        return {
            "G": np.array([0, 0, 0]),
            "F": np.array([1 - phi, 1 - phi, 1 - psi]),
            "F1": np.array([phi, phi - 1, psi]),
            "F2": np.array([1 - phi, -phi, 1 - psi]),
            "H": np.array([zeta, zeta, eta]),
            "H1": np.array([1 - zeta, -zeta, 1 - eta]),
            "H2": np.array([-zeta, -zeta, 1 - eta]),
            "I": np.array([1 / 2, -1 / 2, 1 / 2]),
            "M": np.array([1 / 2, 0, 1 / 2]),
            "N": np.array([1 / 2, 0, 0]),
            "N1": np.array([0, -1 / 2, 0]),
            "X": np.array([1 / 2, -1 / 2, 0]),
            "Y": np.array([mu, mu, delta]),
            "Y1": np.array([1 - mu, -mu, -delta]),
            "Y2": np.array([-mu, -mu, -delta]),
            "Y3": np.array([mu, mu - 1, delta]),
            "Z": np.array([0, 0, 1 / 2]),
        }
    elif variation == "MCLC4":
        return {
            "G": np.array([0, 0, 0]),
            "F": np.array([1 - phi, 1 - phi, 1 - psi]),
            "H": np.array([zeta, zeta, eta]),
            "H1": np.array([1 - zeta, -zeta, 1 - eta]),
            "H2": np.array([-zeta, -zeta, 1 - eta]),
            "I": np.array([1 / 2, -1 / 2, 1 / 2]),
            "M": np.array([1 / 2, 0, 1 / 2]),
            "N": np.array([1 / 2, 0, 0]),
            "N1": np.array([0, -1 / 2, 0]),
            "X": np.array([1 / 2, -1 / 2, 0]),
            "Y": np.array([mu, mu, delta]),
            "Y1": np.array([1 - mu, -mu, -delta]),
            "Y2": np.array([-mu, -mu, -delta]),
            "Y3": np.array([mu, mu - 1, delta]),
            "Z": np.array([0, 0, 1 / 2]),
        }
    elif variation == "MCLC5":
        return {
            "G": np.array([0, 0, 0]),
            "F": np.array([nu, nu, omega]),
            "F1": np.array([1 - nu, 1 - nu, 1 - omega]),
            "F2": np.array([nu, nu - 1, omega]),
            "H": np.array([zeta, zeta, eta]),
            "H1": np.array([1 - zeta, -zeta, 1 - eta]),
            "H2": np.array([-zeta, -zeta, 1 - eta]),
            "I": np.array([rho, 1 - rho, 1 / 2]),
            "I1": np.array([1 - rho, rho - 1, 1 / 2]),
            "L": np.array([1 / 2, 1 / 2, 1 / 2]),
            "M": np.array([1 / 2, 0, 1 / 2]),
            "N": np.array([1 / 2, 0, 0]),
            "N1": np.array([0, -1 / 2, 0]),
            "X": np.array([1 / 2, -1 / 2, 0]),
            "Y": np.array([mu, mu, delta]),
            "Y1": np.array([1 - mu, -mu, -delta]),
            "Y2": np.array([-mu, -mu, -delta]),
            "Y3": np.array([mu, mu - 1, delta]),
            "Z": np.array([0, 0, 1 / 2]),
        }


def TRI_hs_points(variation):
    r"""
    Get high-symmetry points for the TRI lattice.

    See :ref:`guide_tri` for the details.

    Parameters
    ----------
    variation : str
        TRI variation.

    Returns
    -------
    kpoints : dict
        High-symmetry points.
    """

    if variation in ["TRI1a", "TRI2a"]:
        return {
            "G": np.array([0, 0, 0]),
            "L": np.array([1 / 2, 1 / 2, 0]),
            "M": np.array([0, 1 / 2, 1 / 2]),
            "N": np.array([1 / 2, 0, 1 / 2]),
            "R": np.array([1 / 2, 1 / 2, 1 / 2]),
            "X": np.array([1 / 2, 0, 0]),
            "Y": np.array([0, 1 / 2, 0]),
            "Z": np.array([0, 0, 1 / 2]),
        }

    elif variation in ["TRI1b", "TRI2b"]:
        return {
            "G": np.array([0, 0, 0]),
            "L": np.array([1 / 2, -1 / 2, 0]),
            "M": np.array([0, 0, 1 / 2]),
            "N": np.array([-1 / 2, -1 / 2, 1 / 2]),
            "R": np.array([0, -1 / 2, 1 / 2]),
            "X": np.array([0, -1 / 2, 0]),
            "Y": np.array([1 / 2, 0, 0]),
            "Z": np.array([-1 / 2, 0, 1 / 2]),
        }

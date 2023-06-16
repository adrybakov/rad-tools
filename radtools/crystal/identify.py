r"""
Crystal/lattice identification.
"""

from math import acos, cos, floor, log10, sqrt

import numpy as np
from termcolor import cprint

from radtools.routines import (
    cell_from_param,
    get_permutation,
    print_2d_array,
    reciprocal_cell,
    todegrees,
    toradians,
    volume,
)

__all__ = ["niggli", "lepage"]


def niggli(
    a=1,
    b=1,
    c=1,
    alpha=90,
    beta=90,
    gamma=90,
    eps_rel=1e-5,
    verbose=False,
    return_cell=False,
):
    r"""
    Computes Niggli matrix form.

    Parameters
    ----------
    a : float, default 1
        Length of the :math:`a_1` vector.
    b : float, default 1
        Length of the :math:`a_2` vector.
    c : float, default 1
        Length of the :math:`a_3` vector.
    alpha : float, default 90
        Angle between vectors :math:`a_2` and :math:`a_3`. In degrees.
    beta : float, default 90
        Angle between vectors :math:`a_1` and :math:`a_3`. In degrees.
    gamma : float, default 90
        Angle between vectors :math:`a_1` and :math:`a_2`. In degrees.
    eps_rel : float, default 1e-5
        Relative epsilon as defined in [2]_.
    verbose : bool, default False
        Whether to print the steps of an algorithm.
    return_cell : bool, default False
        Whether to return cell parameters instead of Niggli matrix form.

    Returns
    -------
    result : (3,2) :numpy:`ndarray`
        Niggli matrix form as defined in [1]_:

        .. math::

            \begin{pmatrix}
                A & B & C \\
                \xi/2 & \eta/2 & \zeta/2
            \end{pmatrix}
        
        If return_cell == True, then return Niggli cell: (a, b, c, alpha, beta, gamma).


    References
    ----------
    .. [1] Křivý, I. and Gruber, B., 1976.
        A unified algorithm for determining the reduced (Niggli) cell.
        Acta Crystallographica Section A: Crystal Physics, Diffraction,
        Theoretical and General Crystallography,
        32(2), pp.297-298.
    .. [2] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.

    Examples
    --------
    Example from [1]_
    (parameters are reproducing :math:`A=9`, :math:`B=27`, :math:`C=4`, 
    :math:`\xi` = -5, :math:`\eta` = -4, :math:`\zeta = -22`):

    .. doctest::

        >>> import radtools as rad
        >>> from radtools.routines import todegrees
        >>> from math import acos, sqrt
        >>> a = 3
        >>> b = sqrt(27)
        >>> c = 2
        >>> print(f"{a} {b:.3f} {c}")
        3 5.196 2
        >>> alpha = acos(-5 / 2 / b / c) * todegrees
        >>> beta = acos(-4 / 2 / a / c) * todegrees
        >>> gamma = acos(-22 / 2 / a / b) * todegrees
        >>> print(f"{alpha:.2f} {beta:.2f} {gamma:.2f}")
        103.92 109.47 134.88
        >>> niggli_matrix_form = rad.niggli(a, b, c, alpha, beta, gamma, verbose=True)
                       A         B         C        xi        eta      zeta   
        start:       9.00000  27.00000   4.00000  -5.00000  -4.00000 -22.00000
        2 appl. to   9.00000  27.00000   4.00000  -5.00000  -4.00000 -22.00000
        1 appl. to   9.00000   4.00000  27.00000  -5.00000 -22.00000  -4.00000
        4 appl. to   4.00000   9.00000  27.00000 -22.00000  -5.00000  -4.00000
        5 appl. to   4.00000   9.00000  27.00000 -22.00000  -5.00000  -4.00000
        4 appl. to   4.00000   9.00000  14.00000  -4.00000  -9.00000  -4.00000
        6 appl. to   4.00000   9.00000  14.00000  -4.00000  -9.00000  -4.00000
        4 appl. to   4.00000   9.00000   9.00000  -8.00000  -1.00000  -4.00000
        7 appl. to   4.00000   9.00000   9.00000  -8.00000  -1.00000  -4.00000
        3 appl. to   4.00000   9.00000   9.00000  -9.00000  -1.00000   4.00000
        5 appl. to   4.00000   9.00000   9.00000   9.00000   1.00000   4.00000
        3 appl. to   4.00000   9.00000   9.00000  -9.00000  -3.00000   4.00000
        result:      4.00000   9.00000   9.00000   9.00000   3.00000   4.00000
        >>> niggli_matrix_form
        array([[4. , 9. , 9. ],
               [4.5, 1.5, 2. ]])

    """

    eps = eps_rel * volume(a, b, c, alpha, beta, gamma) ** (1 / 3.0)
    n = abs(floor(log10(abs(eps))))

    # 0
    A = a**2
    B = b**2
    C = c**2
    xi = 2 * b * c * cos(alpha * toradians)
    eta = 2 * a * c * cos(beta * toradians)
    zeta = 2 * a * b * cos(gamma * toradians)
    N = (
        max(
            len(str(A).split(".")[0]),
            len(str(B).split(".")[0]),
            len(str(C).split(".")[0]),
            len(str(xi).split(".")[0]),
            len(str(eta).split(".")[0]),
            len(str(zeta).split(".")[0]),
        )
        + 1
        + n
    )

    def compare(x, condition, y):
        if condition == "<":
            return x < y - eps
        if condition == ">":
            return y < x - eps
        if condition == "<=":
            return not y < x - eps
        if condition == ">=":
            return not x < y - eps
        if condition == "==":
            return not (x < y - eps or y < x - eps)

    def summary_line():
        return (
            f"{A:{N}.{n}f} {B:{N}.{n}f} {C:{N}.{n}f} "
            + f"{xi:{N}.{n}f} {eta:{N}.{n}f} {zeta:{N}.{n}f}"
        )

    if verbose:
        print(
            f"           {'A':^{N}} {'B':^{N}} {'C':^{N}} "
            + f"{'xi':^{N}} {'eta':^{N}} {'zeta':^{N}}"
        )
        cprint(
            f"start:     {summary_line()}",
            color="yellow",
        )
        phrase = "appl. to"
    while True:
        # 1
        if compare(A, ">", B) or (
            compare(A, "==", B) and compare(abs(xi), ">", abs(eta))
        ):
            if verbose:
                print(f"1 {phrase} {summary_line()}")
            A, xi, B, eta = B, eta, A, xi
        # 2
        if compare(B, ">", C) or (
            compare(B, "==", C) and compare(abs(eta), ">", abs(zeta))
        ):
            if verbose:
                print(f"2 {phrase} {summary_line()}")
            B, eta, C, zeta = C, zeta, B, eta
            continue
        # 3
        if compare(xi * eta * zeta, ">", 0):
            if verbose:
                print(f"3 {phrase} {summary_line()}")
            xi, eta, zeta = abs(xi), abs(eta), abs(zeta)
        # 4
        if compare(xi * eta * zeta, "<=", 0):
            if verbose:
                print(f"4 {phrase} {summary_line()}")
            xi, eta, zeta = -abs(xi), -abs(eta), -abs(zeta)
        # 5
        if (
            compare(abs(xi), ">", B)
            or (compare(xi, "==", B) and compare(2 * eta, "<", zeta))
            or (compare(xi, "==", -B) and compare(zeta, "<", 0))
        ):
            if verbose:
                print(f"5 {phrase} {summary_line()}")
            C = B + C - xi * np.sign(xi)
            eta = eta - zeta * np.sign(xi)
            xi = xi - 2 * B * np.sign(xi)
            continue
        # 6
        if (
            compare(abs(eta), ">", A)
            or (compare(eta, "==", A) and compare(2 * xi, "<", zeta))
            or (compare(eta, "==", -A) and compare(zeta, "<", 0))
        ):
            if verbose:
                print(f"6 {phrase} {summary_line()}")
            C = A + C - eta * np.sign(eta)
            xi = xi - zeta * np.sign(eta)
            eta = eta - 2 * A * np.sign(eta)
            continue
        # 7
        if (
            compare(abs(zeta), ">", A)
            or (compare(zeta, "==", A) and compare(2 * xi, "<", eta))
            or (compare(zeta, "==", -A) and compare(eta, "<", 0))
        ):
            if verbose:
                print(f"7 {phrase} {summary_line()}")
            B = A + B - zeta * np.sign(zeta)
            xi = xi - eta * np.sign(zeta)
            zeta = zeta - 2 * A * np.sign(zeta)
            continue
        # 8
        if compare(xi + eta + zeta + A + B, "<", 0) or (
            compare(xi + eta + zeta + A + B, "==", 0)
            and compare(2 * (A + eta) + zeta, ">", 0)
        ):
            if verbose:
                print(f"8 {phrase} {summary_line()}")
            C = A + B + C + xi + eta + zeta
            xi = 2 * B + xi + zeta
            eta = 2 * A + eta + zeta
            continue
        break
    if verbose:
        cprint(
            f"result:    {summary_line()}",
            color="green",
        )
    if return_cell:
        a = sqrt(A)
        b = sqrt(B)
        c = sqrt(C)
        alpha = acos(xi / 2 / b / c) * todegrees
        beta = acos(eta / 2 / a / c) * todegrees
        gamma = acos(zeta / 2 / a / b) * todegrees
        return a, b, c, alpha, beta, gamma
    return np.array([[A, B, C], [xi / 2, eta / 2, zeta / 2]])


def check_cub(angles: np.ndarray, axes: np.ndarray, eps):
    target_angles = np.array(
        [
            [0, 45, 45, 45, 45, 90, 90, 90, 90],
            [0, 45, 45, 45, 45, 90, 90, 90, 90],
            [0, 45, 45, 45, 45, 90, 90, 90, 90],
            [0, 45, 45, 60, 60, 60, 60, 90, 90],
            [0, 45, 45, 60, 60, 60, 60, 90, 90],
            [0, 45, 45, 60, 60, 60, 60, 90, 90],
            [0, 45, 45, 60, 60, 60, 60, 90, 90],
            [0, 45, 45, 60, 60, 60, 60, 90, 90],
            [0, 45, 45, 60, 60, 60, 60, 90, 90],
        ]
    )

    conventional_axis = np.array([0, 45, 45, 45, 45, 90, 90, 90, 90])

    axes = np.array([i[0] for i in axes])

    if 9 <= angles.shape[0]:
        indices = get_permutation(angles.shape[0], 9)
        for index in indices:
            sub_angles = angles[np.ix_(index, index)]
            sub_axes = axes[index]
            if (
                np.abs(np.sort(sub_angles.flatten()) - np.sort(target_angles.flatten()))
                < eps
            ).all():
                xyz = []
                for i in range(9):
                    if (np.abs(np.sort(sub_angles[i]) - conventional_axis) < eps).all():
                        xyz.append(sub_axes[i])
                det = np.abs(np.linalg.det(xyz))
                if det == 1:
                    result = "CUB"
                elif det == 4:
                    result = "FCC"
                elif det == 2:
                    result = "BCC"
                return result, False
        return None, True
    return None, True


def check_hex(angles: np.ndarray, eps):
    target_angles = np.array(
        [
            [0, 90, 90, 90, 90, 90, 90],
            [0, 30, 30, 60, 60, 90, 90],
            [0, 30, 30, 60, 60, 90, 90],
            [0, 30, 30, 60, 60, 90, 90],
            [0, 30, 30, 60, 60, 90, 90],
            [0, 30, 30, 60, 60, 90, 90],
            [0, 30, 30, 60, 60, 90, 90],
        ]
    )
    angles = angles[:7]
    if 7 <= angles.shape[0]:
        indices = get_permutation(angles.shape[0], 7)
        for index in indices:
            sub_angles = angles[np.ix_(index, index)]
            if (
                np.abs(np.sort(sub_angles.flatten()) - np.sort(target_angles.flatten()))
                < eps
            ).all():
                return "HEX", False
        return None, True
    return None, True


def check_tet(angles: np.ndarray, axes: np.ndarray, eps, cell):
    target_angles = np.array(
        [
            [0, 90, 90, 90, 90],
            [0, 45, 45, 90, 90],
            [0, 45, 45, 90, 90],
            [0, 45, 45, 90, 90],
            [0, 45, 45, 90, 90],
        ]
    )

    conventional_axis = np.array([0, 90, 90, 90, 90])

    axes = np.array([i[0] for i in axes])

    angles = angles[:5]
    if 5 <= angles.shape[0]:
        indices = get_permutation(angles.shape[0], 5)
        for index in indices:
            sub_angles = angles[np.ix_(index, index)]
            sub_axes = axes[index]
            if (
                np.abs(np.sort(sub_angles.flatten()) - np.sort(target_angles.flatten()))
                < eps
            ).all():
                xy = []
                for i in range(5):
                    if (np.abs(np.sort(sub_angles[i]) - conventional_axis) < eps).all():
                        z = sub_axes[i]
                    else:
                        xy.append(sub_axes[i])
                xy.sort(key=lambda x: np.linalg.norm(x @ cell))

                xyz = [z, xy[0], xy[1]]

                det = np.abs(np.linalg.det(xyz))
                if det == 1:
                    result = "TET"
                elif det == 2:
                    result = "BCT"
                return result, False
        return None, True
    return None, True


def check_rhl(angles: np.ndarray, eps):
    target_angles = np.array(
        [
            [0, 60, 60],
            [0, 60, 60],
            [0, 60, 60],
        ]
    )

    angles = angles[:3]
    if 3 <= angles.shape[0]:
        indices = get_permutation(angles.shape[0], 3)
        for index in indices:
            sub_angles = angles[np.ix_(index, index)]
            if (
                np.abs(np.sort(sub_angles.flatten()) - np.sort(target_angles.flatten()))
                < eps
            ).all():
                return "RHL", False
        return None, True
    return None, True


def check_orc(angles: np.ndarray, axes: np.ndarray, eps):
    target_angles = np.array(
        [
            [0, 90, 90],
            [0, 90, 90],
            [0, 90, 90],
        ]
    )

    angles = angles[:3]
    axes = np.array([i[0] for i in axes])
    if 3 <= angles.shape[0]:
        indices = get_permutation(angles.shape[0], 3)
        for index in indices:
            sub_angles = angles[np.ix_(index, index)]
            sub_axes = axes[index]
            if (
                np.abs(np.sort(sub_angles.flatten()) - np.sort(target_angles.flatten()))
                < eps
            ).all():
                C = np.array(sub_axes, dtype=float).T
                det = np.abs(np.linalg.det(C))
                if det == 1:
                    result = "ORC"
                if det == 4:
                    result = "ORCF"
                if det == 2:
                    v = C @ [1, 1, 1]

                    def gcd(p, q):
                        while q != 0:
                            p, q = q, p % q
                        return p

                    if (
                        gcd(abs(v[0]), abs(v[1])) > 1
                        and gcd(abs(v[0]), abs(v[2])) > 1
                        and gcd(abs(v[1]), abs(v[2])) > 1
                    ):
                        result = "ORCI"
                    else:
                        result = "ORCC"
                return result, False
        return None, True
    return None, True


def get_perpendicular_shortest(v, cell, eps):
    perp_axes = []

    miller_indices = (np.indices((3, 3, 3)) - 1).transpose((1, 2, 3, 0)).reshape(27, 3)

    for index in miller_indices:
        if (index != [0, 0, 0]).any():
            if abs((index @ cell) @ (v @ cell)) < eps:
                perp_axes.append(index)

    perp_axes.sort(key=lambda x: np.linalg.norm(x @ cell))

    # indices 0 and 2 (not 0 and 1), since v and -v are present in miller_indices
    return perp_axes[0], perp_axes[2]


def check_mcl(angles: np.ndarray, axes: np.ndarray, eps, cell):
    axes = np.array([i[0] for i in axes])
    angles = angles[:1]
    if 1 <= angles.shape[0]:
        indices = get_permutation(angles.shape[0], 1)
        for index in indices:
            sub_axes = axes[index]
            b = sub_axes[0]

            a, c = get_perpendicular_shortest(b, cell, eps)

            C = np.array(
                [
                    a,
                    b,
                    c,
                ],
                dtype=float,
            ).T
            det = np.abs(np.linalg.det(C))
            if det == 1:
                return "MCL", False
            if det == 2:
                return "MCLC", False
        return None, True
    return None, True


def lepage(
    a=1,
    b=1,
    c=1,
    alpha=90,
    beta=90,
    gamma=90,
    eps_rel=1e-4,
    verbose=False,
    very_verbose=False,
    give_all_results=False,
    delta_max=0.01,
):
    r"""
    Le Page algorithm [1]_.

    Parameters
    ----------
    a : float, default 1
        Length of the :math:`a_1` vector.
    b : float, default 1
        Length of the :math:`a_2` vector.
    c : float, default 1
        Length of the :math:`a_3` vector.
    alpha : float, default 90
        Angle between vectors :math:`a_2` and :math:`a_3`. In degrees.
    beta : float, default 90
        Angle between vectors :math:`a_1` and :math:`a_3`. In degrees.
    gamma : float, default 90
        Angle between vectors :math:`a_1` and :math:`a_2`. In degrees.
    eps_rel : float, default 1e-4
        Relative epsilon.
    verbose : bool, default False
        Whether to print the steps of an algorithm.
    very_verbose : bool, default False
        Whether to print the detailed steps of an algorithm.
    give_all_results : bool, default False
        Whether to give the whole list of consecutive results.
    delta_max : float, default 0.01
        Maximum angle tolerance, in degrees.

    References
    ----------
    .. [1] Le Page, Y., 1982.
        The derivation of the axes of the conventional unit cell from
        the dimensions of the Buerger-reduced cell.
        Journal of Applied Crystallography, 15(3), pp.255-259.
    """

    if very_verbose:
        verbose = True

    limit = max(1.5, delta_max * 1.1)

    eps_volumetric = eps_rel * volume(a, b, c, alpha, beta, gamma) ** (1 / 3.0)
    decimals = abs(floor(log10(abs(eps_volumetric))))
    if delta_max is None:
        delta_max = eps

    # Niggli reduction
    a, b, c, alpha, beta, gamma = niggli(a, b, c, alpha, beta, gamma, return_cell=True)
    cell = cell_from_param(a, b, c, alpha, beta, gamma)
    rcell = reciprocal_cell(cell)
    if very_verbose:
        print("Cell:")
        print_2d_array(cell, fmt=f"{4+decimals}.{1+decimals}f")
        print("Reciprocal cell:")
        print_2d_array(rcell, fmt=f"{4+decimals}.{1+decimals}f")

    # Find all axes with twins
    miller_indices = (np.indices((5, 5, 5)) - 2).transpose((1, 2, 3, 0)).reshape(125, 3)
    axes = []
    for U in miller_indices:
        for h in miller_indices:
            if abs(U @ h) == 2:
                t = U @ cell
                tau = h @ rcell
                delta = (
                    np.arctan(np.linalg.norm(np.cross(t, tau)) / abs(t @ tau))
                    * todegrees
                )
                if delta < limit:
                    axes.append([U, t / np.linalg.norm(t), abs(U @ h), delta])

    # Sort and filter
    axes.sort(key=lambda x: x[-1])
    keep_index = np.ones(len(axes))
    for i in range(len(axes)):
        if keep_index[i]:
            for j in range(i + 1, len(axes)):
                if (
                    (axes[i][0] == axes[j][0]).all()
                    or (axes[i][0] == -axes[j][0]).all()
                    or (axes[i][0] == 2 * axes[j][0]).all()
                    or (axes[i][0] == 0.5 * axes[j][0]).all()
                ):
                    keep_index[i] = 0
                    break
    new_axes = []
    for i in range(len(axes)):
        if keep_index[i]:
            if set(axes[i][0]) == {0, 2}:
                axes[i][0] = axes[i][0] / 2
            new_axes.append(axes[i])
    axes = new_axes

    if very_verbose:
        print(f"Axes with delta < {limit}:")
        print(f"       U     {'delta':>{4+decimals}}")
        for ax in axes:
            print(
                f"  ({ax[0][0]:2.0f} "
                + f"{ax[0][1]:2.0f} "
                + f"{ax[0][2]:2.0f}) "
                + f"{ax[-1]:{4+decimals}.{1+decimals}f}"
            )

    # Compute angles matrix
    n = len(axes)
    angles = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i, n):
            angles[i][j] = (
                acos(abs(np.clip(np.array(axes[i][1]) @ np.array(axes[j][1]), -1, 1)))
                * todegrees
            )
    angles += angles.T

    # Main check cycle
    delta = None
    separator = lambda x: "=" * 20 + f" Cycle {x} " + "=" * 20
    cycle = 0
    if give_all_results:
        results = []
    while delta is None or delta > delta_max:
        if verbose:
            cycle += 1
            print(separator(cycle))

        try:
            delta = max(axes, key=lambda x: x[-1])[-1]
        except ValueError:
            delta = 0
        eps = max(eps_volumetric, delta)
        decimals = abs(floor(log10(abs(eps))))
        if very_verbose:
            decimals = abs(floor(log10(abs(eps))))
            print(
                f"Epsilon = {eps:{4+decimals}.{1+decimals}f}\n"
                + f"delta   = {delta:{4+decimals}.{1+decimals}f}"
            )
            print("Axes:")
            print(f"       U     {'delta':>{4+decimals}}")
            for ax in axes:
                print(
                    f"  ({ax[0][0]:2.0f} "
                    + f"{ax[0][1]:2.0f} "
                    + f"{ax[0][2]:2.0f}) "
                    + f"{ax[-1]:{4+decimals}.{1+decimals}f}"
                )
            print("Angles:")
            print_2d_array(angles, fmt=f"{4+decimals}.{1+decimals}f")

        continue_search = True
        n = len(axes)
        result = None

        # CUB
        result, continue_search = check_cub(angles, axes, eps)

        # HEX
        if continue_search:
            result, continue_search = check_hex(angles, eps)

        # TET
        if continue_search:
            result, continue_search = check_tet(angles, axes, eps, cell)

        # RHL
        if continue_search:
            result, continue_search = check_rhl(angles, eps)

        # ORC
        if continue_search:
            result, continue_search = check_orc(angles, axes, eps)

        # MCL
        if continue_search:
            result, continue_search = check_mcl(angles, axes, eps, cell)

        # TRI
        if continue_search:
            result = "TRI"

        if verbose:
            print(
                f"System {result} with the worst delta = {delta:{4+decimals}.{1+decimals}f}"
            )

        if len(axes) > 0:
            # remove worst axes
            while len(axes) >= 2 and axes[-1][-1] == axes[-2][-1]:
                axes = axes[:-1]
                angles = np.delete(angles, -1, -1)[:-1]
            axes = axes[:-1]
            angles = np.delete(angles, -1, -1)[:-1]

        if give_all_results:
            results.append((result, delta))

    if give_all_results:
        return results

    return result


if __name__ == "__main__":
    print(
        lepage(
            6.137,
            6.136925044352425,
            20.718,
            90.0,
            90.0,
            119.99501387877993,
            very_verbose=True,
        )
    )
    # from radtools.crystal.bravais_lattice import lattice_example

    # for i, name in enumerate(
    #     [
    #         "CUB",
    #         "FCC",
    #         "BCC",
    #         "HEX",
    #         "TET",
    #         "BCT1",
    #         "BCT2",
    #         "RHL1",
    #         "RHL2",
    #         "ORC",
    #         "ORCF1",
    #         "ORCF2",
    #         "ORCF3",
    #         "ORCI",
    #         "ORCC",
    #         "MCL",
    #         "MCLC1",
    #         "MCLC2",
    #         "MCLC3",
    #         "MCLC4",
    #         "MCLC5",
    #         "TRI1a",
    #         "TRI2a",
    #         "TRI1b",
    #         "TRI2b",
    #     ]
    # ):
    #     lattice = lattice_example(name)

    #     print(
    #         name,
    #         lepage(
    #             lattice.a,
    #             lattice.b,
    #             lattice.c,
    #             lattice.alpha,
    #             lattice.beta,
    #             lattice.gamma,
    #             very_verbose=True,
    #         ),
    #     )

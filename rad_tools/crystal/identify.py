r"""
Crystal/lattice identification.
"""

from math import acos, cos, floor, log10, sqrt

import numpy as np

from termcolor import cprint

from rad_tools.routines import (
    _todegrees,
    _toradians,
    volume,
    cell_from_param,
    reciprocal_cell,
    print_2D_array,
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

        >>> import rad_tools as rad
        >>> from rad_tools.routines import _todegrees
        >>> from math import acos, sqrt
        >>> a = 3
        >>> b = sqrt(27)
        >>> c = 2
        >>> print(f"{a} {b:.3f} {c}")
        3 5.196 2
        >>> alpha = acos(-5 / 2 / b / c) * _todegrees
        >>> beta = acos(-4 / 2 / a / c) * _todegrees
        >>> gamma = acos(-22 / 2 / a / b) * _todegrees
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
    xi = 2 * b * c * cos(alpha * _toradians)
    eta = 2 * a * c * cos(beta * _toradians)
    zeta = 2 * a * b * cos(gamma * _toradians)
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
        a = round(sqrt(A), n)
        b = round(sqrt(B), n)
        c = round(sqrt(C), n)
        alpha = round(acos(xi / 2 / b / c) * _todegrees, n)
        beta = round(acos(eta / 2 / a / c) * _todegrees, n)
        gamma = round(acos(zeta / 2 / a / b) * _todegrees, n)
        return a, b, c, alpha, beta, gamma
    return np.around(np.array([[A, B, C], [xi / 2, eta / 2, zeta / 2]]), decimals=n)


def lepage(
    a=1,
    b=1,
    c=1,
    alpha=90,
    beta=90,
    gamma=90,
    eps_rel=1e-4,
    limit=1.5,
    verbose=False,
    very_verbose=False,
    delta_max=None,
):
    r"""
    Le Page algorithm [1_].

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
        Relative epsilon as defined in [2]_.
    limit : float, default 1.5
        Limit for the Miller index search.
    verbose : bool, default False
        Whether to print the steps of an algorithm.
    very_verbose : bool, default False
        Whether to print the detailed steps of an algorithm.

    References
    ----------
    .. [1] Le Page, Y., 1982.
        The derivation of the axes of the conventional unit cell from
        the dimensions of the Buerger-reduced cell.
        Journal of Applied Crystallography, 15(3), pp.255-259.
    """

    if very_verbose:
        verbose = True

    eps = eps_rel * volume(a, b, c, alpha, beta, gamma) ** (1 / 3.0)
    decimals = 3
    if delta_max is None:
        delta_max = eps

    target_angles = {
        "CUB": np.concatenate(
            (np.zeros(9), 45 * np.ones(24), 60 * np.ones(24), 90 * np.ones(24))
        ),
        "HEX": np.concatenate(
            (np.zeros(7), 30 * np.ones(12), 60 * np.ones(12), 90 * np.ones(18))
        ),
        "TET": np.concatenate((np.zeros(5), 45 * np.ones(8), 90 * np.ones(12))),
        "RHL": np.concatenate((np.zeros(3), 60 * np.ones(6))),
        "ORC": np.concatenate((np.zeros(3), 90 * np.ones(6))),
    }

    conventional_axis = {
        "CUB": np.concatenate((np.zeros(1), 45 * np.ones(4), 90 * np.ones(4))),
        "TET": {
            "z": np.concatenate((np.zeros(1), 90 * np.ones(4))),
            "xy": np.concatenate((np.zeros(1), 45 * np.ones(2), 90 * np.ones(2))),
        },
    }

    # Niggli reduction
    a, b, c, alpha, beta, gamma = niggli(a, b, c, alpha, beta, gamma, return_cell=True)
    cell = cell_from_param(a, b, c, alpha, beta, gamma)
    rcell = reciprocal_cell(cell)

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
                    * _todegrees
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
        cprint(f"Axes with delta < {limit}:", color="yellow")
        print(f"       U     {'delta':>{3+decimals}}")
        for ax in axes:
            print(
                f"  ({ax[0][0]:2.0f} "
                + f"{ax[0][1]:2.0f} "
                + f"{ax[0][0]:2.0f}) "
                + f"{ax[-1]:{3+decimals}.{decimals}f}"
            )

    # Compute angles matrix
    n = len(axes)
    angles = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            angles[i][j] = (
                acos(abs(np.clip(np.array(axes[i][1]) @ np.array(axes[j][1]), -1, 1)))
                * _todegrees
            )

    delta = None
    separator = lambda x: "=" * 20 + f" Cycle {x} " + "=" * 20
    i = 0
    while delta is None or delta >= delta_max:
        if verbose:
            i += 1
            print(separator(i))

        if very_verbose:
            print("Axes:")
            print_2D_array(list(map(lambda x: x[0], axes)), fmt="2.0f")
            print("Angles:")
            print_2D_array(angles, fmt=f"{3+decimals}.{decimals}f")

        try:
            delta = max(axes, key=lambda x: x[-1])[-1]
        except ValueError:
            delta = 0

        continue_search = True
        n = len(axes)
        result = None

        # CUB
        if (
            n**2 == target_angles["CUB"].shape[0]
            and (
                np.abs(np.sort(np.abs(angles.flatten())) - target_angles["CUB"]) < eps
            ).all()
        ):
            xyz = []
            for i in range(n):
                if (
                    np.abs(np.sort(np.abs(angles[i])) - conventional_axis["CUB"]) < eps
                ).all():
                    xyz.append(axes[i])
            det = np.abs(
                np.linalg.det(
                    [
                        xyz[0][0],
                        xyz[1][0],
                        xyz[2][0],
                    ]
                )
            )
            if det == 1:
                result = "CUB"
            elif det == 4:
                result = "FCC"
            elif det == 2:
                result = "BCC"
            continue_search = False

        # HEX
        if continue_search and (
            n**2 == target_angles["HEX"].shape[0]
            and (
                np.abs(np.sort(np.abs(angles.flatten())) - target_angles["HEX"]) < eps
            ).all()
        ):
            result = "HEX"
            continue_search = False

        # TET
        if continue_search and (
            n**2 == target_angles["TET"].shape[0]
            and (
                np.abs(np.sort(np.abs(angles.flatten())) - target_angles["TET"]) < eps
            ).all()
        ):
            x, y, z = None, None, None
            x_alt, y_alt = None, None
            for i in range(n):
                if (
                    np.abs(np.sort(np.abs(angles[i])) - conventional_axis["TET"]["z"])
                    < eps
                ).all():
                    z = axes[i]
                if (
                    np.abs(np.sort(np.abs(angles[i])) - conventional_axis["TET"]["xy"])
                    < eps
                ).all():
                    if x is None:
                        x = axes[i]
                    elif y is None and abs(axes[i][1] @ x[1]) < eps:
                        y = axes[i]
                    elif x_alt is None:
                        x_alt = axes[i]
                    else:
                        y_alt = axes[i]

            if np.linalg.norm(x_alt[0] @ cell) < np.linalg.norm(x[0] @ cell):
                x, y = x_alt, y_alt

            det = np.abs(
                np.linalg.det(
                    [
                        x[0],
                        y[0],
                        z[0],
                    ]
                )
            )
            if det == 1:
                result = "TET"
            elif det == 2:
                result = "BCT"
            continue_search = False

        # RHL
        if continue_search and (
            n**2 == target_angles["RHL"].shape[0]
            and (
                np.abs(np.sort(np.abs(angles.flatten())) - target_angles["RHL"]) < eps
            ).all()
        ):
            result = "RHL"
            continue_search = False

        # ORC
        if continue_search and (
            n**2 == target_angles["ORC"].shape[0]
            and (
                np.abs(np.sort(np.abs(angles.flatten())) - target_angles["ORC"]) < eps
            ).all()
        ):
            C = np.array(
                [
                    axes[0][0],
                    axes[1][0],
                    axes[2][0],
                ],
                dtype=float,
            ).T
            det = np.abs(np.linalg.det(C))
            if det == 1:
                result = "ORC"
            if det == 4:
                result = "ORCF"
            if det == 2:
                v1, v2, v3, v4 = (
                    C @ [0, 1, 1],
                    C @ [1, 0, 1],
                    C @ [1, 1, 0],
                    C @ [1, 1, 1],
                )

                def gcd(p, q):
                    while q != 0:
                        p, q = q, p % q
                    return p

                if (
                    gcd(abs(v4[0]), abs(v4[1])) > 1
                    and gcd(abs(v4[0]), abs(v4[2])) > 1
                    and gcd(abs(v4[1]), abs(v4[2])) > 1
                ):
                    result = "ORCI"
                else:
                    result = "ORCC"
            continue_search = False

        # MCL
        if continue_search and (n == 1):
            v = axes[0][0] @ cell
            a, b, c = cell

            ax = []
            test_ax = []

            if abs(a @ v) < eps:
                ax.append(np.array([1, 0, 0]))
            else:
                test_ax.append(np.array([1, 0, 0]))
            if abs(b @ v) < eps:
                ax.append(np.array([0, 1, 0]))
            else:
                test_ax.append(np.array([0, 1, 0]))
            if abs(c @ v) < eps:
                ax.append(np.array([0, 0, 1]))
            else:
                test_ax.append(np.array([0, 0, 1]))

            indices = [[1, 0], [0, 1], [1, 1], [1, -1]]
            if len(ax) == 2:
                a, b = ax
                ax = [a, b, a + b, a - b]
                ax.sort(key=lambda x: np.linalg.norm(x @ cell))
                a = ax[0]
                b = ax[1]
                c = axes[0][0]
            elif len(test_ax) == 2:
                tmp = ax
                a, b = test_ax
                ax = [a, b, a + b, a - b]
                new_ax = []
                for i in ax:
                    if abs((i @ cell) @ v) < eps:
                        new_ax.append(i)

                a, b = tmp[0], new_ax[0]
                ax = [a, b, a + b, a - b]
                ax.sort(key=lambda x: np.linalg.norm(x @ cell))
                a = ax[0]
                b = ax[1]
                c = axes[0][0]

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
                result = "MCL"
            if det == 2:
                result = "MCLC"
            continue_search = False

        # TRI
        if continue_search:
            result = "TRI"

        if verbose:
            print(
                f"System {result} with the worst delta = {delta:{3+decimals}.{decimals}f}"
            )

        if len(axes) > 0:
            # remove worst axes
            while len(axes) >= 2 and axes[-1][-1] == axes[-2][-1]:
                axes = axes[:-1]
                angles = np.delete(angles, -1, -1)[:-1]
            axes = axes[:-1]
            angles = np.delete(angles, -1, -1)[:-1]

    return result


if __name__ == "__main__":
    print(
        lepage(
            4,
            4.472,
            4.583,
            79.030,
            64.130,
            64.150,
            very_verbose=True,
            eps_rel=0.001,
            delta_max=0.006,
        )
    )
    print(
        lepage(
            4,
            4,
            4,
            90,
            90,
            90,
            very_verbose=True,
            eps_rel=0.001,
            delta_max=0.006,
        )
    )

    # from rad_tools.crystal.bravais_lattice import lattice_example

    # delta_max = [
    #     None,
    #     None,
    #     None,
    #     None,
    #     None,
    #     None,
    #     None,
    #     1e-5,
    #     None,
    #     None,
    #     None,
    #     None,
    #     None,
    #     None,
    #     None,
    #     None,
    #     None,
    #     None,
    #     None,
    #     None,
    #     1e-5,
    #     None,
    #     None,
    #     None,
    #     None,
    # ]
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
    #     print("\n" + "=" * 80)
    #     print(
    #         name,
    #         lepage(
    #             lattice.a,
    #             lattice.b,
    #             lattice.c,
    #             lattice.alpha,
    #             lattice.beta,
    #             lattice.gamma,
    #             verbose=True,
    #             delta_max=delta_max[i],
    #         ),
    #     )

    # print("\n" + "=" * 80)

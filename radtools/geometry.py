from math import cos, sqrt

import numpy as np
from scipy.spatial.transform import Rotation

from radtools.constants import TODEGREES, TORADIANS
from radtools.crystal.constants import ABS_TOL, ABS_TOL_ANGLE
from radtools.numerical import compare_numerically

__all__ = [
    "volume",
    "angle",
    "parallelepiped_check",
    "absolute_to_relative",
    "span_orthonormal_set",
]


def volume(*args):
    r"""
    Computes volume.

    .. versionadded:: 0.7

    Three type of arguments are expected:

    * One argument.
        Matrix ``cell``.
        Volume is computed as:

        .. math::
            V = v_1 \cdot (v_2 \times v_3)
    * Three arguments.
        Vectors ``v1``, ``v2``, ``v3``.
        Volume is computed as:

        .. math::
            V = v_1 \cdot (v_2 \times v_3)
    * Six arguments.
        Parameters ``a``, ``b``, ``c``, ``alpha``, ``beta``, ``gamma``.
        Volume is computed as:

        .. math::
            V = abc\sqrt{1+2\cos(\alpha)\cos(\beta)\cos(\gamma)-\cos^2(\alpha)-\cos^2(\beta)-\cos^2(\gamma)}

    Notes
    -----

    Volume can be negative if the three vectors are not right-handed.

    Parameters
    ----------
    v1 : (3,) |array_like|_
        First vector.
    v2 : (3,) |array_like|_
        Second vector.
    v3 : (3,) |array_like|_
        Third vector.
    cell : (3,3) |array_like|_
        Cell matrix, rows are interpreted as vectors.
    a : float, default 1
        Length of the :math:`v_1` vector.
    b : float, default 1
        Length of the :math:`v_2` vector.
    c : float, default 1
        Length of the :math:`v_3` vector.
    alpha : float, default 90
        Angle between vectors :math:`v_2` and :math:`v_3`. In degrees.
    beta : float, default 90
        Angle between vectors :math:`v_1` and :math:`v_3`. In degrees.
    gamma : float, default 90
        Angle between vectors :math:`v_1` and :math:`v_2`. In degrees.

    Returns
    -------
    volume : float
        Volume of corresponding region in space.
    """

    if len(args) == 1:
        cell = np.array(args[0])
    elif len(args) == 3:
        cell = np.array(args)
    elif len(args) == 6:
        a, b, c, alpha, beta, gamma = args
        alpha = alpha * TORADIANS
        beta = beta * TORADIANS
        gamma = gamma * TORADIANS
        sq_root = (
            1
            + 2 * cos(alpha) * cos(beta) * cos(gamma)
            - cos(alpha) ** 2
            - cos(beta) ** 2
            - cos(gamma) ** 2
        )
        return a * b * c * sqrt(sq_root)
    else:
        raise ValueError(
            "Unable to identify input. "
            + "Supported: one (3,3) array_like, or three (3,) array_like, or 6 floats."
        )

    return np.linalg.det(cell)


def angle(v1, v2, radians=False):
    r"""
    Angle between two vectors.

    .. versionadded:: 0.7

    .. math::

        \alpha = \dfrac{(\vec{v_1} \cdot \vec{v_2})}{\vert\vec{v_1}\vert\cdot\vert\vec{v_2}\vert}

    Parameters
    ----------
    v1 : (3,) |array_like|_
        First vector.
    v2 : (3,) |array_like|_
        Second vector.
    radians : bool, default False
        Whether to return value in radians. Return value in degrees by default.

    Returns
    -------
    angle: float
        Angle in degrees or radians (see ``radians``).

    Raises
    ------
    ValueError
        If one of the vectors is zero vector (or both). Norm is compared against
        :numpy:`finfo`\ (float).eps.
    """

    # Normalize vectors
    v1_norm = np.linalg.norm(v1)
    v2_norm = np.linalg.norm(v2)
    if abs(v1_norm) <= np.finfo(float).eps or abs(v2_norm) <= np.finfo(float).eps:
        raise ValueError("Angle is ill defined (zero vector).")

    v1 = np.array(v1) / v1_norm
    v2 = np.array(v2) / v2_norm

    alpha = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
    if radians:
        return alpha
    return alpha * TODEGREES


def parallelepiped_check(a, b, c, alpha, beta, gamma, raise_error=False):
    r"""
    Check if parallelepiped is valid.

    Parameters
    ----------
    a : float
        Length of the :math:`v_1` vector.
    b : float
        Length of the :math:`v_2` vector.
    c : float
        Length of the :math:`v_3` vector.
    alpha : float
        Angle between vectors :math:`v_2` and :math:`v_3`. In degrees.
    beta : float
        Angle between vectors :math:`v_1` and :math:`v_3`. In degrees.
    gamma : float
        Angle between vectors :math:`v_1` and :math:`v_2`. In degrees.
    raise_error : bool, default False
        Whether to raise error if parameters could not form a parallelepiped.

    Returns
    -------
    result: bool
        Whether the parameters could from a parallelepiped.

    Raises
    ------
    ValueError
        If parameters could not form a parallelepiped.
        Only raised if ``raise_error`` is ``True`` (it is ``False`` by default).
    """

    result = (
        compare_numerically(a, ">", 0.0, ABS_TOL)
        and compare_numerically(b, ">", 0.0, ABS_TOL)
        and compare_numerically(c, ">", 0.0, ABS_TOL)
        and compare_numerically(alpha, "<", 180.0, ABS_TOL_ANGLE)
        and compare_numerically(beta, "<", 180.0, ABS_TOL_ANGLE)
        and compare_numerically(gamma, "<", 180.0, ABS_TOL_ANGLE)
        and compare_numerically(alpha, ">", 0.0, ABS_TOL_ANGLE)
        and compare_numerically(beta, ">", 0.0, ABS_TOL_ANGLE)
        and compare_numerically(gamma, ">", 0.0, ABS_TOL_ANGLE)
        and compare_numerically(gamma, "<", alpha + beta, ABS_TOL_ANGLE)
        and compare_numerically(alpha + beta, "<", 360.0 - gamma, ABS_TOL_ANGLE)
        and compare_numerically(beta, "<", alpha + gamma, ABS_TOL_ANGLE)
        and compare_numerically(alpha + gamma, "<", 360.0 - beta, ABS_TOL_ANGLE)
        and compare_numerically(alpha, "<", beta + gamma, ABS_TOL_ANGLE)
        and compare_numerically(beta + gamma, "<", 360.0 - alpha, ABS_TOL_ANGLE)
    )

    if not result and raise_error:
        message = "Parameters could not form a parallelepiped:\n"
        message += f"a = {a}"
        if not compare_numerically(a, ">", 0.0, ABS_TOL):
            message += f" (a <= 0 with numerical tolerance: {ABS_TOL})"
        message += "\n"
        message += f"b = {b}"
        if not compare_numerically(b, ">", 0.0, ABS_TOL):
            message += f" (b <= 0 with numerical tolerance: {ABS_TOL})"
        message += "\n"
        message += f"c = {c}"
        if not compare_numerically(c, ">", 0.0, ABS_TOL):
            message += f" (c <= 0 with numerical tolerance: {ABS_TOL})"
        message += "\n"
        message += f"alpha = {alpha}\n"
        if not compare_numerically(alpha, "<", 180.0, ABS_TOL_ANGLE):
            message += f"  (alpha >= 180 with numerical tolerance: {ABS_TOL_ANGLE})\n"
        if not compare_numerically(alpha, ">", 0.0, ABS_TOL_ANGLE):
            message += f"  (alpha <= 0 with numerical tolerance: {ABS_TOL_ANGLE})\n"
        message += f"beta = {beta}\n"
        if not compare_numerically(beta, "<", 180.0, ABS_TOL_ANGLE):
            message += f"  (beta >= 180 with numerical tolerance: {ABS_TOL_ANGLE})\n"
        if not compare_numerically(beta, ">", 0.0, ABS_TOL_ANGLE):
            message += f"  (beta <= 0 with numerical tolerance: {ABS_TOL_ANGLE})\n"
        message += f"gamma = {gamma}\n"
        if not compare_numerically(gamma, "<", 180.0, ABS_TOL_ANGLE):
            message += f"  (gamma >= 180 with numerical tolerance: {ABS_TOL_ANGLE})\n"
        if not compare_numerically(gamma, ">", 0.0, ABS_TOL_ANGLE):
            message += f"  (gamma <= 0 with numerical tolerance: {ABS_TOL_ANGLE})\n"
        if not compare_numerically(gamma, "<", alpha + beta, ABS_TOL_ANGLE):
            message += f"Inequality gamma < alpha + beta is not satisfied with numerical tolerance: {ABS_TOL_ANGLE}\n"
        if not compare_numerically(alpha + beta, "<", 360.0 - gamma, ABS_TOL_ANGLE):
            message += f"Inequality alpha + beta < 360 - gamma is not satisfied with numerical tolerance: {ABS_TOL_ANGLE}\n"
        if not compare_numerically(beta, "<", alpha + gamma, ABS_TOL_ANGLE):
            message += f"Inequality beta < alpha + gamma is not satisfied with numerical tolerance: {ABS_TOL_ANGLE}\n"
        if not compare_numerically(alpha + gamma, "<", 360.0 - beta, ABS_TOL_ANGLE):
            message += f"Inequality alpha + gamma < 360 - beta is not satisfied with numerical tolerance: {ABS_TOL_ANGLE}\n"
        if not compare_numerically(alpha, "<", beta + gamma, ABS_TOL_ANGLE):
            message += f"Inequality alpha < beta + gamma is not satisfied with numerical tolerance: {ABS_TOL_ANGLE}\n"
        if not compare_numerically(beta + gamma, "<", 360.0 - alpha, ABS_TOL_ANGLE):
            message += f"Inequality beta + gamma < 360 - alpha is not satisfied with numerical tolerance: {ABS_TOL_ANGLE}\n"
        raise ValueError(message)

    return result


def absolute_to_relative(basis, vector):
    r"""
    Compute relative coordinates of the vector with respect to the basis.

    Parameters
    ----------
    basis : (3, 3) |array_like|_
        Lattice vectors.
    vector : (3,) |array_like|_
        Absolute coordinates.

    Returns
    -------
    relative : (3,) :numpy:`ndarray`
        Relative coordinates.
    """

    # Three vectors of the cell
    v1 = np.array(basis[0], dtype=float)
    v2 = np.array(basis[1], dtype=float)
    v3 = np.array(basis[2], dtype=float)

    v = np.array(vector, dtype=float)
    if (v == np.zeros(3)).all():
        return np.zeros(3)

    # Compose system of linear equations
    B = np.array([np.dot(v1, v), np.dot(v2, v), np.dot(v3, v)])
    A = np.array(
        [
            [np.dot(v1, v1), np.dot(v1, v2), np.dot(v1, v3)],
            [np.dot(v2, v1), np.dot(v2, v2), np.dot(v2, v3)],
            [np.dot(v3, v1), np.dot(v3, v2), np.dot(v3, v3)],
        ]
    )

    # Solve and return
    return np.linalg.solve(A, B)


def span_orthonormal_set(vec):
    r"""
    Span orthonormal set of vectors.

    Parameters
    ----------
    vec : (3,) |array_like|_
        Vector, which serves as :math:`e_3`

    Returns
    -------
    e1 : (3,) :numpy:`ndarray`
    e2 : (3,) :numpy:`ndarray`
    e3 : (3,) :numpy:`ndarray`
    """

    vec = np.array(vec) / np.linalg.norm(vec)

    if np.allclose(vec, [0, 0, 1]):
        return (
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
        )

    if np.allclose(vec, [0, 0, -1]):
        return (
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, -1.0, 0.0]),
            np.array([0.0, 0.0, -1.0]),
        )

    z_dir = [0, 0, 1]
    n = (
        np.cross(vec, z_dir)
        / np.linalg.norm(np.cross(vec, z_dir))
        * np.arccos(np.dot(vec, z_dir) / np.linalg.norm(vec))
    )
    rotation_matrix = Rotation.from_rotvec(n).as_matrix()

    return rotation_matrix

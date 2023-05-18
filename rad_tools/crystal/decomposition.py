r"""
The algorithm for defining the planes, edges and corners of the Brillouin zone (or Wigner-Seitz cell):

Planes
------
* Generate all closest lattice points.
    By varying the :math:`i`, :math:`j`, :math:`k` indexes in the range :math:`(-1,0,1)`. 
    The point (0,0,0) is excluded. 
    Each vector  of the lattice point is defined as :math:`\vec{v} = i\vec{v}_1 + j\vec{v}_2 + k\vec{v}_3`.
* For each lattice point check if the middle point :math:`\vec{v} / 2` is closer to :math:`\Gamma` point then to any of the other lattice points. 
* If yes, then the vector :math:`\vec{v}/2` defines the plane.

Corners
-------
* Define potential corners by computing intersection points for each group of three planes.
* For each point check if it is closer (or at the same distance) to :math:`\Gamma`, then to any of the other lattice points.
* If yes, then take the point.

Edges
-----
* For each pair of corners check if they share two planes.
* If yes, then take an edge.
"""


import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from matplotlib import rcParams

TOLERANCE = 1e-8
TOL_BASE = 8

__all__ = [
    "deduct_zone",
    "get_lattice_points",
    "define_planes",
    "define_corners",
    "define_edges",
]


def deduct_zone(cell=None, v1=None, v2=None, v3=None):
    r"""
    Construct Brillouin zone (or Wigner-Seitz cell), from the given set of basis vectors.

    Parameters
    ----------
    cell : (3,3) array_like, default None
        Matrix of the basis vectors, rows are interpreted as vectors,
        columns as cartesian coordinates:

        .. code-block:: python

            cell = [
                [v1_x, v1_y, v1_z],
                [v2_x, v2_y, v2_z],
                [v3_x, v3_y, v3_z]
                ]
    v1 : (3,) |array_like|_, default None
        First basis vector (first primitive lattice vector)
    v2 : (3,) |array_like|_, default None
        Second basis vector (second primitive lattice vector)
    v2 : (3,) |array_like|_, default None
        Third basis vector (third primitive lattice vector)

    Returns
    -------
    planes:  list
        List of N planes, constraining Brillouin zone (or Wigner-Seitz cell).
        Plane is defined by the vector :math:`v`, which is perpendicular to the plane
        and gives the coordinate of the point from the plane at the same time.
        The vector is given in absolute coordinates.
        Units are the same as in the input (``cell`` or ``v1``, ``v2``, ``v3``).

        .. code-block:: python

            planes = [v1, ...vN]

        where :math:`v = (v_x, v_y, v_z)`.
    corners : list
        List of M corners of the Brillouin zone (or Wigner-Seitz cell).
        Corner is defined by the vector :math:`v = v_x, v_y, v_z` in absolute coordinates.
        Units are the same as in the input (``cell`` or ``v1``, ``v2``, ``v3``).
    edges: list
        List of pair of point, each entry define one edge.

        .. code-block:: python

            edges = [(p1, p2), ...]

        where :math:`p = (p_x, p_y, p_z)` in absolute coordinates.
        Units are the same as in the input (``cell`` or ``v1``, ``v2``, ``v3``).

    See Also
    --------
    get_lattice_points
    define_planes
    define_corners
    define_edges
    """

    lattice_points, vectors = get_lattice_points(cell)
    planes = define_planes(lattice_points, vectors, relative=False)
    corners, plane_indices = define_corners(planes, vectors)
    edges = define_edges(corners, plane_indices)
    return planes, edges, corners


def get_lattice_points(cell=None, v1=None, v2=None, v3=None):
    r"""
    Get all lattice points and corresponding vectors. Does not include (0,0,0).

    Parameters
    ----------
    cell : (3,3) array_like, default None
        Matrix of the basis vectors, rows are interpreted as vectors,
        columns as cartesian coordinates:

        .. code-block:: python

            cell = [
                [v1_x, v1_y, v1_z],
                [v2_x, v2_y, v2_z],
                [v3_x, v3_y, v3_z]
                ]
    v1 : (3,) |array_like|_, default None
        First basis vector (first primitive lattice vector)
    v2 : (3,) |array_like|_, default None
        Second basis vector (second primitive lattice vector)
    v2 : (3,) |array_like|_, default None
        Third basis vector (third primitive lattice vector)

    Returns
    -------
    lattice_points : list
        List of all lattice points, constructed by the permutations of (i,j,k)
        in the range (-1,0,1). Without (0,0,0) point.

        .. code-block:: python

            lattice_points = [p_1, ..., p_26]

        where :math:`p = (i,j,k)` in relative coordinates with respect to
        ``cell`` or ``v1``, ``v2``, ``v3``.
    vectors : (26, 3) :numpy:`array`
        Array of the vectors corresponding to the lattice points.

        .. code-block:: python

            vectors = [v_1, ... v_26]

        where :math:`v = (v_x, v_y, v_z)` in absolute coordinates.
        Units are the same as in the input (``cell`` or ``v1``, ``v2``, ``v3``).
    """

    # Check if provided data are sufficient.
    if cell is None and (v1 is None or v2 is None or v3 is None):
        raise ValueError(
            "Either cell or three vectors (v1, v2, v3) have to be passed to the function."
        )
    if cell is None:
        cell = np.array([v1, v2, v3])
    else:
        cell = np.array(cell)

    # Check if the cell is not coplanar
    if np.dot(cell[0], np.cross(cell[1], cell[2])) == 0:
        raise ValueError(
            f"Three provided vectors are coplanar.\n"
            + f"  v1 = {cell[0]}\n  v2 = {cell[1]}\n  v3 = {cell[2]}\n"
        )

    # Compute plattice_points
    lattice_points = []
    vectors = np.zeros((26, 3))
    index = 0
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            for k in [-1, 0, 1]:
                if i**2 + j**2 + k**2 != 0:
                    lattice_points.append((i, j, k))
                    vectors[index] = np.matmul(np.array((i, j, k)), cell)
                    index += 1

    return lattice_points, vectors


def define_planes(lattice_points, vectors, relative=False):
    r"""
    Define planes, which borders the first brillouin zone (or Wigner-Seitz cell).

    Parameters
    ----------
    lattice_points : list
        List of all lattice points, constructed by the permutations of (i,j,k)
        in the range (-1,0,1). Without (0,0,0) point.

        .. code-block:: python

            lattice_points = [p_1, ..., p_26]

        where :math:`p = (i,j,k)`.
    vectors : (26, 3) :numpy:`array`
        Array of the vectors corresponding to the lattice points.

        .. code-block:: python

            vectors = [v_1, ... v_26]

        where :math:`v = (v_x, v_y, v_z)`.
    relative : bool, default False.
        Whether to return planes in relative or absolute coordinates.

    Returns
    -------
    planes:  list
        List of N planes, constraining Brillouin zone (or Wigner-Seitz cell).
        Plane is defined by the vector :math:`v`, which is perpendicular to the plane
        and gives the coordinate of the point from the plane. If ``relative = True``,
        then :math:`v = (i,j,k)` in relative coordinates. If ``relative = False``,
        then :math:`v = (v_x,v_y,v_z)` in absolute coordinates.
        Units are the same as in the input (``vectors``).

    See Also
    --------
    get_lattice_points
    define_corners
    define_edges
    """

    points = np.transpose(
        np.tile(vectors / 2, (len(vectors), 1)).reshape(len(vectors), len(vectors), 3),
        (1, 0, 2),
    )
    lpoints = np.tile(vectors, (len(vectors), 1)).reshape(len(vectors), len(vectors), 3)

    compare_matrix = np.around(
        np.linalg.norm(points, axis=2) - np.linalg.norm(points - lpoints, axis=2),
        decimals=TOL_BASE,
    )
    compare_matrix = np.sum((compare_matrix >= 0) * 1, axis=1)

    planes = []

    for i in range(compare_matrix.shape[0]):
        if compare_matrix[i] <= 1:
            if relative:
                planes.append(
                    (
                        lattice_points[i][0] / 2,
                        lattice_points[i][1] / 2,
                        lattice_points[i][2] / 2,
                    )
                )
            else:
                planes.append(vectors[i] / 2)

    return planes


def define_corners(planes, vectors):
    r"""
    Define corners of the first brillouin zone (or Wigner-Seitz cell).

    Parameters
    ----------
    planes:  list
        List of N planes, constraining Brillouin zone (or Wigner-Seitz cell).
        Plane is defined by the vector :math:`v = (v_x,v_y,v_z)` in absolute coordinates,
        which is perpendicular to the plane and gives the coordinate of the point from the plane.
    vectors : (26, 3) :numpy:`array`
        Array of the vectors corresponding to the lattice points.

        .. code-block:: python

            vectors = [v_1, ... v_26]

        where :math:`v = (v_x, v_y, v_z)`.

    Returns
    -------
    corners : list
        List of M corners of the Brillouin zone (or Wigner-Seitz cell).
        Corner is defined by the vector :math:`v = v_x, v_y, v_z` in absolute coordinates.
        Units are the same as in the input (``vectors`` and ``planes``).
    plane_indices : list
        Indices of the three planes, which intersection produced the corners.

    See Also
    --------
    get_lattice_points
    define_planes
    define_edges
    """

    # Compute all 3-plane intersections
    intersection_points = []
    planes_indices = []
    for f, fplane in enumerate(planes):
        nf = fplane
        for s in range(f + 1, len(planes)):
            splane = planes[s]
            ns = splane
            for t in range(s + 1, len(planes)):
                tplane = planes[t]
                nt = tplane
                # If not around - overflow error occurs
                A = np.around(np.array([nf, ns, nt]), decimals=TOL_BASE)
                b = np.around(
                    np.array(
                        [
                            np.linalg.norm(nf) ** 2,
                            np.linalg.norm(ns) ** 2,
                            np.linalg.norm(nt) ** 2,
                        ]
                    ),
                    decimals=TOL_BASE,
                )
                try:
                    x = np.linalg.solve(A, b)
                    xprime = sp.linalg.solve(A, b)
                    intersection_points.append(x)
                    planes_indices.append({f, s, t})
                except:
                    pass

    corners = []

    # Check if intersection point is closer (or at the same distance)
    # to the Gamma point then to any other point of lattice.
    ipoints = np.transpose(
        np.tile(intersection_points, (len(vectors), 1)).reshape(
            len(vectors), len(intersection_points), 3
        ),
        (1, 0, 2),
    )
    lpoints = np.tile(vectors, (len(intersection_points), 1)).reshape(
        len(intersection_points), len(vectors), 3
    )

    compare_matrix = np.around(
        np.linalg.norm(ipoints, axis=2) - np.linalg.norm(ipoints - lpoints, axis=2),
        decimals=TOL_BASE,
    )

    compare_matrix = np.sum((compare_matrix > 0) * 1, axis=1)

    corners = []
    tmp = []
    for i in range(compare_matrix.shape[0]):
        if compare_matrix[i] == 0:
            corners.append(intersection_points[i])
            tmp.append(planes_indices[i])

    # Filter non-unique entries
    unique_corners = []
    planes_indices = []
    for i in range(len(corners)):
        takeit = True
        for j in range(len(unique_corners)):
            if (np.linalg.norm(corners[i] - unique_corners[j]) < TOLERANCE).all():
                takeit = False
                break
        if takeit:
            unique_corners.append(corners[i])
            planes_indices.append(tmp[i])
        else:
            planes_indices[j] = planes_indices[j].union(tmp[i])

    return unique_corners, planes_indices


def define_edges(corners, plane_indices):
    r"""
    Define edges of the first brillouin zone (or Wigner-Seitz cell).

    Parameters
    ----------
    corners : list
        List of M corners of the Brillouin zone (or Wigner-Seitz cell).
        Corner is defined by the vector :math:`v = v_x, v_y, v_z` in absolute coordinates.
    plane_indices : list
        Indices of the three planes, which intersection produced the corners.

    Returns
    -------
    edges: list
        List of pair of point, each entry define one edge.

        .. code-block:: python

            edges = [(p1, p2), ...]

        where :math:`p = (p_x, p_y, p_z)` in absolute coordintes.
        Units are the same as in the input (``corners``).

    See Also
    --------
    get_lattice_points
    define_planes
    define_corners
    """

    edges = []
    for f in range(len(corners)):
        for s in range(f + 1, len(corners)):
            if len(set(plane_indices[f]).intersection(set(plane_indices[s]))) == 2:
                edges.append([corners[f], corners[s]])

    return edges


if __name__ == "__main__":
    from rad_tools import orcf2

    cell = orcf2.reciprocal_cell
    lattice_points, vectors = get_lattice_points(cell)
    planes = define_planes(lattice_points, vectors)
    corners, planes_indices = define_corners(planes, vectors)
    edges = define_edges(corners, planes_indices)
    print(len(edges))

r"""
The algorithm for defining the planes, edges and corners of the Brillouin zone:

Planes
------
* Generate all closest lattice points.
    By varying the :math:`i`, :math:`j`, :math:`k` indexes in the range :math:`(-1,0,1)`. 
    The point (0,0,0) is excluded. 
    Each vector  of the lattice point is defined as :math:`\vec{v} = i\vec{v}_1 + j\vec{v}_2 + k\vec{v}_3`.
* For each lattice point check if the middle point :math:`\vec{v} / 2` is closer 
    to :math:`\Gamma` point then to any of the other lattice points. 
* If yes, then the vector :math:`\vec{v}/2` defines the plane.

Corners
=======
* Define potential corners by computing intersection points for each group of three planes.
* For each point check if it is closer to :math:`\Gamma`, then to any of the other lattice points.
* If yes, then take the point.

Edges
=====
* For each pair of 
"""


import numpy as np

TOLERANCE = 1e-8
TOL_BASE = 8


def deduct_zone(cell):
    r"""
    Construct Brillouin zone (or Wigner-Seitz cell), from the given set of basis vectors.

    Parameters
    ----------
    cell : (3,3) array_like
        Matrix of the basis vectors, rows are interpreted as vectors,
        columns as cartesian coordinates:

        .. code-block:: python
            cell = [
                [b1_x, b1_y, b1_z],
                [b2_x, b2_y, b2_z],
                [b3_x, b3_y, b3_z]
                ]
    Returns
    -------
    planes: (N, 3) array
        Array of N planes, constraining Brillouin zone (or Wigner-Seitz cell).
        Plane is defined by the vector :math:`v`, which is perpendicular to the plane
        and gives the coordinate of the point from the plane.
        The vector is given in relative coordinates, with respect to the basis vectors
    edges : (K, 2, 3) array
        Array of K edges of the Brillouin zone (or Wigner-Seitz cell).
        Define by the two points.
    corners: (M, 3) array
        Array of M corners of the Brillouin zone (or Wigner-Seitz cell).
    """

    planes, other_points = define_planes(cell)

    # Compute corners
    corners_candidates = []
    for f, fplane in enumerate(planes):
        nf = np.matmul(fplane, cell) / 2
        for s in range(f + 1, len(planes)):
            splane = planes[s]
            ns = np.matmul(splane, cell) / 2
            for t in range(s + 1, len(planes)):
                tplane = planes[t]
                nt = np.matmul(tplane, cell) / 2
                A = np.array([nf, ns, nt])
                b = np.array(
                    [
                        np.linalg.norm(nf) ** 2,
                        np.linalg.norm(ns) ** 2,
                        np.linalg.norm(nt) ** 2,
                    ]
                )
                try:
                    x = np.linalg.solve(A, b)
                    corners_candidates.append(x)
                except:
                    pass

    corners = []

    # Check if corner is closer to the Gamma point then to any other point of lattice.
    for corner in corners_candidates:
        takeit = True
        for other_point in other_points:
            if (
                np.linalg.norm(corner) - np.linalg.norm(corner - other_point)
            ) > TOLERANCE:
                takeit = False
                break
        if takeit:
            corners.append(corner)

    # Filter equal entries
    new_corners = []
    for i in range(len(corners)):
        takeit = True
        for j in range(len(new_corners)):
            if (np.abs(corners[i] - new_corners[j]) < TOLERANCE).all():
                takeit = False
                break
        if takeit:
            new_corners.append(corners[i])
    corners = new_corners

    # Compute corners
    corners_planes = []
    for corner in corners:
        corners_planes.append([])
        for plane in planes:
            A, B, C = tuple(np.matmul(plane, cell) / 2)
            D = -np.linalg.norm(np.matmul(plane, cell) / 2) ** 2
            if abs(A * corner[0] + B * corner[1] + C * corner[2] + D) < TOLERANCE:
                corners_planes[-1].append(plane)
    edges = []
    for fc, fcorner in enumerate(corners):
        for sc in range(fc + 1, len(corners)):
            scorner = corners[sc]
            n = 0
            for fplane in corners_planes[fc]:
                for splane in corners_planes[sc]:
                    if np.linalg.norm(fplane - splane) < TOLERANCE:
                        n += 1
            if n == 2:
                edges.append([fcorner, scorner])
    return planes, edges, corners


def get_lattice_points_vectors(cell=None, v1=None, v2=None, v3=None):
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
        where :math:`p = (i,j,k)`.
    vectors : (26, 3) :numpy:`array`
        Array of the vectors corresponding to the lattice points.

        .. code-block:: python

            vectors = [v_1, ... v_26]
        where :math:`v = (v_x, v_y, v_z)`.
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


def define_planes(cell=None, v1=None, v2=None, v3=None):
    r"""
    Define planes, which borders the first brillouin zone (or Wigner-Seitz cell).

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
        List of N planes (i,j,k), constraining Brillouin zone (or Wigner-Seitz cell).
        Plane is defined by the vector :math:`v`, which is perpendicular to the plane
        and gives the coordinate of the point from the plane.
        The vector is given in relative coordinates, with respect to the basis vectors
    """

    lattice_points, vectors = get_lattice_points_vectors(cell=cell, v1=v1, v2=v2, v3=v3)

    points = np.transpose(
        np.tile(vectors / 2, (len(vectors), 1)).reshape(len(vectors), len(vectors), 3),
        (1, 0, 2),
    )
    lpoints = np.tile(vectors, (len(vectors), 1)).reshape(len(vectors), len(vectors), 3)

    compare_matrix = np.around(
        np.linalg.norm(points, axis=2) - np.linalg.norm(points - lpoints, axis=2),
        decimals=TOL_BASE,
    )

    planes = []

    for i in range(compare_matrix.shape[0]):
        if compare_matrix[i] <= 1:
            planes.append(lattice_points[i])

    return planes


if __name__ == "__main__":
    cell = [[-2, 1.6, 1.2], [2, -1.6, 1.2], [2, 1.6, -1.2]]
    planes = define_planes(cell=cell)
    for i in planes:
        print(i)

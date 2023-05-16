r"""
Routines defining lattice type.
"""


import numpy as np

TOLERANCE = 10e-8


def deduct_zone(
    cell,
):
    r"""
    Construct Brillouin zone (or Wigner-Seitz cell), from the given set of basis vectors.

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

    cell = np.array(cell)

    # Compute planes
    planes_candidates = []
    other_points = []
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            for k in [-1, 0, 1]:
                if i**2 + j**2 + k**2 != 0:
                    rel_vector = np.array((i, j, k))
                    vector = np.matmul(rel_vector, cell)
                    distance_gamma = np.linalg.norm(vector) / 2
                    planes_candidates.append((rel_vector, vector, distance_gamma))
                    other_points.append(vector)

    planes_candidates.sort(key=lambda x: x[2])
    planes = []
    for plane in planes_candidates:
        takeit = True
        for other_point in other_points:
            if (
                not (other_point == plane[1]).all()
                and (plane[2] - np.linalg.norm(plane[1] / 2 - other_point)) > -TOLERANCE
            ):
                takeit = False
                break
        if takeit:
            planes.append(plane[0])

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

import pytest
from radtools.crystal.kpoints import Kpoints
import numpy as np


points = {
    "G": [0, 0, 0],
    "K": [0.5, 0.5, 0],
    "R": [0.5, 0.5, 0.5],
    "X": [0.5, 0, 0],
    "Y": [0, 0.5, 0],
    "Z": [0, 0, 0.5],
    "E": [0.5, 0.5, 1],
    "F": [0.5, 0, 0.5],
    "A": [0, 0.5, 0.5],
    "Q": [0.5, 0.5, 0.5],
}

labels = {
    "G": R"$\Gamma$",
    "K": "K",
    "R": "R",
    "X": "X",
    "Y": "Y",
    "Z": "Z",
    "E": "E",
    "F": "F",
    "A": "A",
    "Q": "Q",
}
paths = [
    [["G", "K", "R"]],
    [["X", "Y", "Z"], ["R", "E"]],
    [["X", "Y", "Z"], ["R", "E"], ["F", "A", "Q"]],
]
path_strings = ["G-K-R", "X-Y-Z|R-E", "X-Y-Z|R-E|F-A-Q"]
correct_labels = [
    [R"$\Gamma$", "K", "R"],
    ["X", "Y", "Z|R", "E"],
    ["X", "Y", "Z|R", "E|F", "A", "Q"],
]
b1 = [1, 1, 0]
b2 = [0, 1, 1]
b3 = [0, 0, 1]
correct_coordinates = [
    np.array([0, 0.7071067811865476, 1.2071067811865475]),
    np.array([0, 0.7071067811865476, 1.4142135623730951, 1.9142135623730951]),
    np.array(
        [
            0,
            0.7071067811865476,
            1.4142135623730951,
            1.9142135623730951,
            2.621320343559643,
            3.121320343559643,
        ]
    ),
]

correct_points = [
    np.array(
        [
            [0, 0, 0, 0],  # G
            [0.1, 0.1, 0, 0.14142135623730953],
            [0.2, 0.2, 0, 0.28284271247461906],
            [0.3, 0.3, 0, 0.42426406871192857],
            [0.4, 0.4, 0, 0.5656854249492381],
            [0.5, 0.5, 0, 0.7071067811865477],  # K
            [0.5, 0.5, 0, 0.7071067811865477],  # K
            [0.5, 0.5, 0.1, 0.8071067811865477],
            [0.5, 0.5, 0.2, 0.90710678118654777],
            [0.5, 0.5, 0.3, 1.0071067811865477],
            [0.5, 0.5, 0.4, 1.1071067811865477],
            [0.5, 0.5, 0.5, 1.2071067811865477],  # R
        ]
    ),
    np.array(
        [
            [0.5, 0, 0, 0],  # X
            [0.4, 0.1, 0, 0.14142135623730953],
            [0.3, 0.2, 0, 0.28284271247461906],
            [0.2, 0.3, 0, 0.42426406871192857],
            [0.1, 0.4, 0, 0.5656854249492381],
            [0, 0.5, 0, 0.7071067811865477],  # Y
            [0, 0.5, 0, 0.7071067811865477],  # Y
            [0, 0.4, 0.1, 0.8485281374238571],
            [0, 0.3, 0.2, 0.9899494936611667],
            [0, 0.2, 0.3, 1.1313708498984762],
            [0, 0.1, 0.4, 1.2727922061357857],
            [0, 0, 0.5, 1.4142135623730954],  # Z
            [0.5, 0.5, 0.5, 1.4142135623730954],  # R
            [0.5, 0.5, 0.6, 1.5142135623730954],
            [0.5, 0.5, 0.7, 1.6142135623730954],
            [0.5, 0.5, 0.8, 1.7142135623730954],
            [0.5, 0.5, 0.9, 1.8142135623730954],
            [0.5, 0.5, 1, 1.9142135623730954],  # E
        ]
    ),
    np.array(
        [
            [0.5, 0, 0, 0],  # X
            [0.4, 0.1, 0, 0.14142135623730953],
            [0.3, 0.2, 0, 0.28284271247461906],
            [0.2, 0.3, 0, 0.42426406871192857],
            [0.1, 0.4, 0, 0.5656854249492381],
            [0, 0.5, 0, 0.7071067811865477],  # Y
            [0, 0.5, 0, 0.7071067811865477],  # Y
            [0, 0.4, 0.1, 0.8485281374238571],
            [0, 0.3, 0.2, 0.9899494936611667],
            [0, 0.2, 0.3, 1.1313708498984762],
            [0, 0.1, 0.4, 1.2727922061357857],
            [0, 0, 0.5, 1.4142135623730954],  # Z
            [0.5, 0.5, 0.5, 1.4142135623730954],  # R
            [0.5, 0.5, 0.6, 1.5142135623730954],
            [0.5, 0.5, 0.7, 1.6142135623730954],
            [0.5, 0.5, 0.8, 1.7142135623730954],
            [0.5, 0.5, 0.9, 1.8142135623730954],
            [0.5, 0.5, 1, 1.9142135623730954],  # E
            [0.5, 0, 0.5, 1.9142135623730954],  # F
            [0.4, 0.1, 0.5, 2.055634918610405],
            [0.3, 0.2, 0.5, 2.1970562748477143],
            [0.2, 0.3, 0.5, 2.338477631085024],
            [0.1, 0.4, 0.5, 2.4798989873223336],
            [0, 0.5, 0.5, 2.621320343559643],  # A
            [0, 0.5, 0.5, 2.621320343559643],  # A
            [0.1, 0.5, 0.5, 2.721320343559643],
            [0.2, 0.5, 0.5, 2.821320343559643],
            [0.3, 0.5, 0.5, 2.921320343559643],
            [0.4, 0.5, 0.5, 3.021320343559643],
            [0.5, 0.5, 0.5, 3.121320343559643],  # Q
        ]
    ),
]

path_input = []
label_input = []
coord_input = []
point_input = []
flat_point_input = []
for i in range(len(paths)):
    path_input.append((paths[i], path_strings[i]))
    label_input.append((paths[i], correct_labels[i]))
    coord_input.append((paths[i], correct_coordinates[i]))
    point_input.append((paths[i], correct_points[i][:, :3]))
    flat_point_input.append((paths[i], correct_points[i][:, 3]))


@pytest.mark.parametrize("path, corr_path", path_input)
def test_path(path, corr_path):
    kp = Kpoints(
        b1,
        b2,
        b3,
        [points[i] for i in points],
        names=[i for i in points],
        labels=[labels[i] for i in points],
        n=4,
        path=path,
    )
    assert kp.path_string == corr_path


@pytest.mark.parametrize("path, corr_lab", label_input)
def test_labels(path, corr_lab):
    kp = Kpoints(
        b1,
        b2,
        b3,
        [points[i] for i in points],
        names=[i for i in points],
        labels=[labels[i] for i in points],
        n=4,
        path=path,
    )
    assert kp.labels == corr_lab


@pytest.mark.parametrize("path, corr_coord", coord_input)
def test_coordinates(path, corr_coord):
    kp = Kpoints(
        b1,
        b2,
        b3,
        [points[i] for i in points],
        names=[i for i in points],
        n=4,
        path=path,
    )
    assert (np.abs(kp.coordinates(relative=True) - corr_coord) < 1e-5).all()


@pytest.mark.parametrize("path, corr_points", point_input)
def test_points(path, corr_points):
    kp = Kpoints(
        b1,
        b2,
        b3,
        [points[i] for i in points],
        names=[i for i in points],
        n=4,
        path=path,
    )
    assert (np.abs(kp.points(relative=True) - corr_points) < 1e-5).all()


@pytest.mark.parametrize("path, corr_flat_points", flat_point_input)
def test_flatten_points(path, corr_flat_points):
    kp = Kpoints(
        b1,
        b2,
        b3,
        [points[i] for i in points],
        names=[i for i in points],
        n=4,
        path=path,
    )
    assert (np.abs(kp.flatten_points(relative=True) - corr_flat_points) < 1e-5).all()

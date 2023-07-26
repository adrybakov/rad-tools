from math import acos, sqrt

import numpy as np
import pytest

from radtools.crystal.bravais_lattice import Lattice, lattice_example
from radtools.crystal.identify import lepage, niggli
from radtools.routines import todegrees


def test_niggli():
    a = 3
    b = sqrt(27)
    c = 2
    alpha = acos(-5 / 2 / sqrt(27) / 2) * todegrees
    beta = acos(-4 / 2 / 3 / 2) * todegrees
    gamma = acos(-22 / 2 / 3 / sqrt(27)) * todegrees
    assert (
        np.array([[4, 9, 9], [9 / 2, 3 / 2, 2]]) - niggli(a, b, c, alpha, beta, gamma)
        < 1e-5
    ).all()


def test_niggli_run():
    alpha = 79.030
    beta = 64.130
    gamma = 64.150
    a = 4
    b = 4.472
    c = 4.583
    ap, bp, cp, alphap, betap, gammap = niggli(
        a, b, c, alpha, beta, gamma, return_cell=True
    )

    assert a - ap < 1e-3
    assert b - bp < 1e-3
    assert c - cp < 1e-3
    assert alpha - alphap < 1e-3
    assert beta - betap < 1e-3
    assert gamma - gammap < 1e-3


lattices = lattice_example()


@pytest.mark.parametrize("name", lattices, ids=lattices)
def test_lepage(name):
    lattice = lattice_example(name)
    type_name = ""
    for i in name:
        if i in "12345":
            break
        type_name += i

    assert (
        lepage(
            lattice.a,
            lattice.b,
            lattice.c,
            lattice.alpha,
            lattice.beta,
            lattice.gamma,
        )
        == type_name
    )


@pytest.mark.parametrize(
    "cell, name, eps_rel",
    [
        (
            [[3.588, 0.000, 0.000], [0.000, 4.807, 0.000], [0.000, 0.000, 23.571]],
            "ORC",
            1e-5,
        ),
        (
            [[6.137, 0.000, 0.000], [-3.068, 5.315, 0.000], [0.000, 0.000, 20.718]],
            "HEX",
            1e-3,
        )
    ],
    ids=["crsbr", "nii2", "custom"],
)
def test_custom_lepage(cell, name, eps_rel):
    lattice = Lattice(cell)
    assert (
        lepage(
            lattice.a,
            lattice.b,
            lattice.c,
            lattice.alpha,
            lattice.beta,
            lattice.gamma,
            eps_rel=eps_rel,
        )
        == name
    )


def test_lepage_paper():
    results = lepage(
        4,
        4.472,
        4.583,
        79.030,
        64.130,
        64.150,
        give_all_results=True,
        eps_rel=0.001,
        delta_max=0.006,
    )
    assert results[0][0] == "BCT"
    assert results[0][1] - 1.482 < 0.001
    assert results[1][0] == "ORCF"
    assert results[1][1] - 1.480 < 0.001
    assert results[2][0] == "ORCF"
    assert results[2][1] - 0.714 < 0.001
    assert results[3][0] == "MCLC"
    assert results[3][1] - 0.714 < 0.001
    assert results[4][0] == "MCLC"
    assert results[4][1] - 0.005 < 0.001


@pytest.mark.parametrize(
    "a, b, c, alpha, beta, gamma, ltype",
    [
        (1, 1, 2, 60, 90, 90, "TET"),
        (1, 1, 1, 30, 90, 90, "ORCC"),
        (1, 1, 1, 45, 90, 90, "ORCC"),
        (3.2, 1, 1, 45, 90, 90, "ORCC"),
        (1, 1, 1, 60, 90, 90, "HEX"),
    ],
)
def test_lepage_custom(a, b, c, alpha, beta, gamma, ltype):
    assert ltype == lepage(
        a,
        b,
        c,
        alpha,
        beta,
        gamma,
        eps_rel=0.001,
    )

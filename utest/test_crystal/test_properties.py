import numpy as np
import pytest

from radtools.crystal.atom import Atom
from radtools.crystal.crystal import Crystal
from radtools.crystal.properties import CONSTANT, dipole_dipole_energy


def test_custom():
    crystal = Crystal()
    crystal.cell = [
        [3.513012, 0.000000000, 0.000000000],
        [0.000000000, 4.752699, 0.000000000],
        [0.000000000, 0.000000000, 23.5428674],
    ]

    Cr1 = Atom(
        "Cr1",
        np.array((0.2500000000, 0.7500000000, 0.1616620796)),
    )
    Cr2 = Atom(
        "Cr2",
        np.array((0.7500000000, 0.2500000000, 0.0737774367)),
    )

    crystal.add_atom(Cr1)
    crystal.add_atom(Cr2)

    crystal.Cr1.magmom = [0, 0, 3]
    crystal.Cr2.magmom = [0, 0, 3]
    z10 = crystal.mag_dipdip_energy(10, 10, 1, progress_bar=False)
    z15 = crystal.mag_dipdip_energy(15, 15, 1, progress_bar=False)
    z20 = crystal.mag_dipdip_energy(20, 20, 1, progress_bar=False)
    assert abs(z10 - 0.032578560794) < 1e-4
    assert abs(z15 - 0.036420741905) < 1e-4
    assert abs(z20 - 0.038660829368) < 1e-4

    crystal.Cr1.magmom = [0, 3, 0]
    crystal.Cr2.magmom = [0, 3, 0]
    y10 = crystal.mag_dipdip_energy(10, 10, 1, progress_bar=False)
    y15 = crystal.mag_dipdip_energy(15, 15, 1, progress_bar=False)
    y20 = crystal.mag_dipdip_energy(20, 20, 1, progress_bar=False)
    assert abs(y10 - -0.012159728056) < 1e-4
    assert abs(y15 - -0.01369598845) < 1e-4
    assert abs(y20 - -0.014599888887) < 1e-4

    crystal.Cr1.magmom = [3, 0, 0]
    crystal.Cr2.magmom = [3, 0, 0]
    x10 = crystal.mag_dipdip_energy(10, 10, 1, progress_bar=False)
    x15 = crystal.mag_dipdip_energy(15, 15, 1, progress_bar=False)
    x20 = crystal.mag_dipdip_energy(20, 20, 1, progress_bar=False)
    assert abs(x10 - -0.020418832737) < 1e-4
    assert abs(x15 - -0.022724753454) < 1e-4
    assert abs(x20 - -0.024060940481) < 1e-4


def test_simple():
    e = dipole_dipole_energy(
        [[[0, 0, 1], [0, 0, 0]], [[0, 0, 1], [1, 0, 0]]], progress_bar=False
    )
    assert abs(e) - CONSTANT / 2 < 1e-5
    e = dipole_dipole_energy(
        [[[0, 0, 1], [0, 0, 0]], [[0, 0, 1], [1, 0, 0]]], progress_bar=False
    )
    assert abs(e) - CONSTANT / 2 < 1e-5
    e = dipole_dipole_energy(
        [[[0, 1, 0], [0, 0, 0]], [[0, 1, 0], [1, 0, 0]]], progress_bar=False
    )
    assert abs(e) - CONSTANT / 2 < 1e-5
    e = dipole_dipole_energy(
        [[[1, 0, 0], [0, 0, 0]], [[1, 0, 0], [1, 0, 0]]], progress_bar=False
    )
    assert abs(e) - 2 * CONSTANT / 2 < 1e-5
    e = dipole_dipole_energy(
        [[[1, 0, 0], [0, 0, 0]], [[0, 1, 0], [1, 0, 0]]], progress_bar=False
    )
    assert abs(e) < 1e-5

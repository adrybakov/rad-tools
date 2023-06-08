import numpy as np
import pytest

from radtools import Atom


def test_atom():
    atom = Atom()
    assert (atom.position == np.zeros(3)).all()
    assert atom.name == "X"
    atom.name = "Cr1"
    assert atom.type == "Cr"
    atom.name = "CR1"
    assert atom.type == "Cr"
    atom.name = "cr1"
    assert atom.type == "Cr"
    with pytest.raises(ValueError):
        a = atom.spin
    with pytest.raises(ValueError):
        a = atom.spin_vector
    with pytest.raises(ValueError):
        a = atom.magmom
    with pytest.raises(ValueError):
        atom.spin_vector = 23
    with pytest.raises(ValueError):
        atom.spin_vector = "adsfasdfs"
    with pytest.raises(ValueError):
        atom.spin_vector = (4, 5)
    with pytest.raises(ValueError):
        atom.spin_vector = (4, 5, 4, 5)
    atom.spin_vector = (0, 0, 1 / 2)
    assert atom.spin - 0.5 < 1e-10
    assert (atom.spin_vector - np.array([0, 0, 0.5]) < 1e-10).all()
    assert (atom.spin_direction - np.array([0, 0, 1]) < 1e-10).all()
    atom.spin = 2
    assert (atom.spin_vector - np.array([0, 0, 2]) < 1e-10).all()
    assert (atom.spin_direction - np.array([0, 0, 1]) < 1e-10).all()
    assert atom.spin - 2 < 1e-10

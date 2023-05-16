import numpy as np
import pytest

from rad_tools import Atom


def test_atom():
    atom = Atom()
    assert (atom.position == np.zeros(3)).all()
    assert atom.literal == "X"
    with pytest.raises(ValueError):
        a = atom.spin
    with pytest.raises(ValueError):
        a = atom.magmom
    with pytest.raises(ValueError):
        atom.spin = 23
    with pytest.raises(ValueError):
        atom.spin = "adsfasdfs"
    with pytest.raises(ValueError):
        atom.spin = (4, 5)
    with pytest.raises(ValueError):
        atom.spin = (4, 5, 4, 5)
    atom.spin = (0, 0, 1 / 2)
    assert (atom.spin == np.array([0, 0, 1 / 2])).all()
    atom.spin = [0, 0, 1 / 2]
    assert (atom.spin == np.array([0, 0, 1 / 2])).all()

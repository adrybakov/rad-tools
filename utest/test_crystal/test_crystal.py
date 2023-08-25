import pytest

from radtools.crystal.crystal import Crystal
from radtools.crystal.atom import Atom
from hypothesis import given, strategies as st, settings

import numpy as np


def test_deepcopy():
    from copy import deepcopy

    c = Crystal()
    a = deepcopy(c)


@given(
    st.text(min_size=1, max_size=10),
    st.lists(st.floats(min_value=0, max_value=1), min_size=3, max_size=3),
)
def test_add_atom(name, position):
    c = Crystal(cell=[[1, 0, 0], [0, 2, 0], [0, 0, 3]], standardize=False)
    c.add_atom(name=name, position=position)
    assert c.atoms[0].name == name
    assert np.allclose(c.atoms[0].position, position)
    c.add_atom(Atom(name, position))
    assert c.atoms[1].name == name
    assert np.allclose(c.atoms[1].position, position)
    c.add_atom(name, position=position)
    assert c.atoms[2].name == name
    assert np.allclose(c.atoms[2].position, position)
    c.add_atom(name, position=position, relative=False)
    assert c.atoms[3].name == name
    assert np.allclose(
        c.atoms[3].position, [position[0], position[1] / 2.0, position[2] / 3.0]
    )


def test_add_atom_raises():
    c = Crystal()

    with pytest.raises(TypeError):
        c.add_atom(nam=4)

    with pytest.raises(TypeError):
        c.add_atom(1)


@given(
    st.text(min_size=1, max_size=10),
    st.lists(st.floats(min_value=0, max_value=1), min_size=3, max_size=3),
)
def test_get_atom(name, position):
    if not name.startswith("__") and not name.endswith("__"):
        position = np.array(position)
        c = Crystal(cell=[[1, 0, 0], [0, 2, 0], [0, 0, 3]], standardize=False)
        c.add_atom(name=name, position=position)
        c.add_atom(name, position=position * 2.0)
        atom = c.get_atom(name, 2)
        assert np.allclose(atom.position, position * 2.0)
        with pytest.raises(ValueError):
            c.get_atom(name)

        atoms = c.get_atom(name, return_all=True)
        assert len(atoms) == 2


def test_remove_atom():
    c = Crystal()
    c.add_atom(Atom("H", (0, 0, 0)))
    c.add_atom(Atom("H", (1, 0, 0)))
    c.add_atom(Atom("O", (2, 0, 0)))
    c.remove_atom("H", index=1)
    assert c.atoms[0].name == "H"
    assert c.atoms[0].index == 2
    assert np.allclose(c.atoms[0].position, (1, 0, 0))
    c = Crystal()
    c.add_atom(Atom("H", (0, 0, 0)))
    c.add_atom(Atom("H", (1, 0, 0)))
    c.add_atom(Atom("O", (2, 0, 0)))
    O = Atom("O", (0, 4, 0))
    c.add_atom(O)
    c.remove_atom("H")
    assert c.atoms[0].name == "O"
    assert c.atoms[0].index == 3
    assert np.allclose(c.atoms[0].position, (2, 0, 0))
    c.remove_atom(O)
    assert len(c.atoms) == 1

    with pytest.raises(ValueError):
        c.remove_atom("H", index=1)


@given(
    st.text(min_size=1, max_size=10),
    st.lists(st.floats(min_value=0, max_value=1), min_size=3, max_size=3),
    st.lists(st.integers(min_value=-1e9, max_value=1e9), min_size=3, max_size=3),
)
def test_get_atom_coordinates(name, position, R):
    position = np.array(position)
    c = Crystal(cell=[[1, 0, 0], [0, 2, 0], [0, 0, 3]], standardize=False)
    c.add_atom(name=name, position=position)
    assert np.allclose(c.get_atom_coordinates(name, R=R), position + np.array(R))
    assert np.allclose(
        c.get_atom_coordinates(name, R=R, relative=False),
        (position + np.array(R)) @ c.cell,
    )


@given(
    st.text(min_size=1, max_size=10),
    st.lists(st.floats(min_value=0, max_value=1), min_size=3, max_size=3),
    st.text(min_size=1, max_size=10),
    st.lists(st.floats(min_value=0, max_value=1), min_size=3, max_size=3),
    st.lists(st.integers(min_value=-1e9, max_value=1e9), min_size=3, max_size=3),
)
def test_get_vector(name1, position1, name2, position2, R):
    position1 = np.array(position1)
    position2 = np.array(position2)
    c = Crystal(cell=[[1, 0, 0], [0, 2, 0], [0, 0, 3]], standardize=False)
    c.add_atom(name=name1, position=position1)
    c.add_atom(name=name2, position=position2)
    assert np.allclose(
        c.get_vector(name1, name2, R=R, index1=1, index2=2),
        (position2 - position1 + np.array(R)) @ c.cell,
    )
    assert np.allclose(
        c.get_vector(name1, name2, R=R, index1=1, index2=2, relative=True),
        position2 - position1 + np.array(R),
    )


@given(
    st.text(min_size=1, max_size=10),
    st.lists(st.floats(min_value=0, max_value=1), min_size=3, max_size=3),
    st.text(min_size=1, max_size=10),
    st.lists(st.floats(min_value=0, max_value=1), min_size=3, max_size=3),
    st.lists(st.integers(min_value=-1e9, max_value=1e9), min_size=3, max_size=3),
)
def test_get_distance(name1, position1, name2, position2, R):
    position1 = np.array(position1)
    position2 = np.array(position2)
    c = Crystal(cell=[[1, 0, 0], [0, 2, 0], [0, 0, 3]], standardize=False)
    c.add_atom(name=name1, position=position1)
    c.add_atom(name=name2, position=position2)
    assert np.allclose(
        c.get_distance(name1, name2, R=R, index1=1, index2=2),
        np.linalg.norm((position2 - position1 + np.array(R)) @ c.cell),
    )
    assert np.allclose(
        c.get_distance(name1, name2, R=R, index1=1, index2=2, relative=True),
        np.linalg.norm(position2 - position1 + np.array(R)),
    )

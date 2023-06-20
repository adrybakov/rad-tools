import pytest

from radtools.crystal.crystal import Crystal
from radtools.crystal.atom import Atom


def test_deepcopy():
    from copy import deepcopy

    c = Crystal()
    a = deepcopy(c)


def test_get_atom():
    c = Crystal()
    c.add_atom(Atom("H", (0, 0, 0)))
    assert c.get_atom("H").name == "H"
    c.add_atom(Atom("H", (1, 0, 0)))
    assert c.get_atom("H", 1).name == "H"
    assert c.get_atom("H", 2).name == "H"
    assert c.get_atom("H__2").name == "H"

import pytest

from radtools.crystal.crystal import Crystal


def test_deepcopy():
    from copy import deepcopy

    c = Crystal()
    a = deepcopy(c)

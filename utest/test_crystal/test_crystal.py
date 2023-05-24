import pytest

from rad_tools.crystal.crystal import Crystal


def test_deepcopy():
    from copy import deepcopy

    c = Crystal()
    a = deepcopy(c)

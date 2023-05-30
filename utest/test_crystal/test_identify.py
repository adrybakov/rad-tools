from math import sqrt

import numpy as np
import pytest
import numpy as np
from math import acos, sqrt
from rad_tools.crystal.identify import niggli, lepage
from rad_tools.routines import _todegrees


def test_niggli():
    a = 3
    b = sqrt(27)
    c = 2
    alpha = acos(-5 / 2 / sqrt(27) / 2) * _todegrees
    beta = acos(-4 / 2 / 3 / 2) * _todegrees
    gamma = acos(-22 / 2 / 3 / sqrt(27)) * _todegrees
    assert (
        np.array([[4, 9, 9], [9 / 2, 3 / 2, 2]]) == niggli(a, b, c, alpha, beta, gamma)
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

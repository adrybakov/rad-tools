from math import sqrt

import pytest

from rad_tools.routines import (
    absolute_to_relative,
    atom_mark_to_latex,
    rot_angle,
    angle,
)


@pytest.mark.parametrize(
    "mark, new_mark",
    [
        ("Cr1", "$Cr_{1}$"),
        ("Cr11", "$Cr_{11}$"),
    ],
)
def test_atom_mark_to_latex(mark, new_mark):
    assert atom_mark_to_latex(mark) == new_mark


class TestRotAngle:
    @pytest.mark.parametrize(
        "x, y, angle",
        [
            (1, 0, 0),
            (1, 1, 45),
            (0, 1, 90),
            (-1, 1, 135),
            (-1, 0, 180),
            (-1, -1, 225),
            (0, -1, 270),
            (1, -1, 315),
        ],
    )
    def test_dummy(self, x, y, angle):
        assert round(rot_angle(x, y, dummy=True), 4) == round(angle, 4)

    @pytest.mark.parametrize(
        "x, y, angle",
        [
            (1, 0, 0),
            (1, 1, 45),
            (0, 1, 90),
            (-1, 1, -45),
            (-1, 0, 0),
            (-1, -1, 45),
            (0, -1, 90),
            (1, -1, -45),
        ],
    )
    def test_not_dummy(self, x, y, angle):
        assert round(rot_angle(x, y), 4) == round(angle, 4)

    def test_ill_case(self):
        with pytest.raises(ValueError):
            rot_angle(0, 0)


@pytest.mark.parametrize(
    "cell, absolute, relative",
    [
        ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0], [0, 0, 0]),
        ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 1], [0, 0, 1]),
        ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 1, 0], [0, 1, 0]),
        ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [1, 0, 0], [1, 0, 0]),
        ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.5, 0.5, 0], [0.5, 0.5, 0]),
        ([[1, 1, 0], [0, 1, 0], [0, 0, 1]], [0.5, 1, 0], [0.5, 0.5, 0]),
        ([[2, 1, 0], [1, 1, 0], [0, 0, 1]], [0.9, 0.7, 0.4], [0.2, 0.5, 0.4]),
    ],
)
def test_absolute_to_relative(cell, absolute, relative):
    new_relative = absolute_to_relative(cell, *tuple(absolute))
    assert (new_relative == relative).all()


def test_angle():
    tol_angle = 1e-5
    assert abs(angle([1, 0, 0], [0, 1, 0]) - 90) < tol_angle
    assert abs(angle([1, 0, 0], [0, 0, 1]) - 90) < tol_angle
    assert abs(angle([0, 1, 0], [1, 0, 0]) - 90) < tol_angle
    assert abs(angle([0, 1, 0], [0, 0, 1]) - 90) < tol_angle
    assert abs(angle([0, 0, 1], [1, 0, 0]) - 90) < tol_angle
    assert abs(angle([0, 0, 1], [0, 1, 0]) - 90) < tol_angle
    assert abs(angle([1, 1, 0], [0, 0, 1]) - 90) < tol_angle
    assert abs(angle([1, 0, 0], [1, 1, 0]) - 45) < tol_angle
    assert abs(angle([1, 0, 0], [1, sqrt(3), 0]) - 60) < tol_angle
    assert abs(angle([1, 0, 0], [sqrt(3), 1, 0]) - 30) < tol_angle
    assert abs(angle([1, 1, 0], [0, 1, 1]) - 60) < tol_angle
    assert abs(angle([1, 1, 0], [3, 1, 1]) - 31.48215) < tol_angle
    assert abs(angle([3, 1, 1], [-3, -1, -1]) - 180) < tol_angle
    assert abs(angle([3, 1, 1], [3, 1, 1])) < tol_angle

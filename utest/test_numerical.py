import pytest

from radtools.numerical import compare_numerically


@pytest.mark.parametrize(
    "x, sign, y, result",
    [
        (-4, "==", -4.0, True),
        (-4.0, "==", -4, True),
        (-4.0, "==", -4.0, True),
        (-4.0, "==", 4.0, False),
    ],
)
def test_compare_numerically(x, sign, y, result):
    assert compare_numerically(x, sign, y) == result

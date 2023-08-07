from math import sqrt

import pytest

from radtools.spinham.parameter import ExchangeParameter

import numpy as np

from hypothesis import given, strategies as st, example
from hypothesis.extra.numpy import arrays as harrays

MAX_MODULUS = 1e8


@given(
    harrays(
        np.float64,
        (3, 3),
        elements=st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
    )
)
def test_init_from_matrix(matrix):
    parameter = ExchangeParameter(matrix=matrix)
    assert np.allclose(parameter.matrix, matrix)


@given(st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS))
def test_init_from_iso(iso):
    parameter = ExchangeParameter(iso=iso)
    assert np.allclose(parameter.matrix, parameter.iso * np.eye(3))


@given(
    harrays(
        np.float64,
        (3,),
        elements=st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
    )
)
def test_init_from_dmi(dmi):
    parameter = ExchangeParameter(dmi=dmi)
    assert np.allclose(
        parameter.matrix,
        parameter.dmi_matrix,
    )


@given(
    st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
    st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
    st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
)
def test_init_from_aniso(a, b, c):
    aniso = np.array([[0, a, b], [a, 0, c], [b, c, 0]])
    parameter = ExchangeParameter(aniso=aniso)
    assert np.allclose(parameter.matrix, parameter.aniso)


@given(
    harrays(
        np.float64,
        (3, 3),
        elements=st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
    )
)
def test_T(matrix):
    parameter = ExchangeParameter(matrix=matrix)
    assert np.allclose(parameter.T.matrix, matrix.T)


@given(
    harrays(
        np.float64,
        (3, 3),
        elements=st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
    )
)
def test_symm_matrix(matrix):
    parameter = ExchangeParameter(matrix=matrix)
    assert np.allclose(parameter.symm_matrix, (matrix + matrix.T) / 2)


@given(
    harrays(
        np.float64,
        (3, 3),
        elements=st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
    )
)
def test_asymm_matrix(matrix):
    parameter = ExchangeParameter(matrix=matrix)
    assert np.allclose(parameter.asymm_matrix, (matrix - matrix.T) / 2)


@given(
    harrays(
        np.float64,
        (3, 3),
        elements=st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
    ),
    st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
)
def test_iso(matrix, new_iso):
    parameter = ExchangeParameter(matrix=matrix)
    assert np.allclose(parameter.iso, np.trace(matrix) / 3.0)
    assert np.allclose(
        parameter.iso_matrix, np.eye(3, dtype=float) * np.trace(matrix) / 3.0
    )

    # Test setter
    parameter.iso = new_iso
    assert np.allclose(parameter.iso, new_iso, atol=1e-7)


@given(
    harrays(
        np.float64,
        (3, 3),
        elements=st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
    ),
    st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
    st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
    st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
)
def test_aniso(matrix, a, b, c):
    parameter = ExchangeParameter(matrix=matrix)
    aniso = (matrix + matrix.T) / 2.0 - np.eye(3, dtype=float) * np.trace(matrix) / 3.0
    assert np.allclose(parameter.aniso, aniso)
    assert np.allclose(parameter.aniso_diagonal, np.diag(aniso))
    assert np.allclose(parameter.aniso_diagonal_matrix, np.diag(np.diag(aniso)))

    # Test setter
    new_aniso = np.array([[0, a, b], [a, 0, c], [b, c, 0]])
    parameter.aniso = new_aniso
    assert np.allclose(parameter.aniso, new_aniso)


@given(
    harrays(
        np.float64,
        (3, 3),
        elements=st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
    ),
    harrays(
        np.float64,
        (3,),
        elements=st.floats(min_value=-MAX_MODULUS, max_value=MAX_MODULUS),
    ),
)
def test_dmi(matrix, new_dmi):
    parameter = ExchangeParameter(matrix=matrix)
    dmi_matrix = (matrix - matrix.T) / 2.0
    dmi = [dmi_matrix[1, 2], dmi_matrix[2, 0], dmi_matrix[0, 1]]
    assert np.allclose(parameter.dmi, dmi)
    assert np.allclose(parameter.dmi_matrix, dmi_matrix)
    assert np.allclose(parameter.dmi_module, np.linalg.norm(dmi))
    if abs(parameter.iso) > 1e-10:
        assert np.allclose(parameter.rel_dmi, parameter.dmi_module / abs(parameter.iso))

    # Test setter
    if (
        max(np.linalg.norm(new_dmi), np.linalg.norm(dmi)) > 1e-8
        and min(np.linalg.norm(new_dmi), np.linalg.norm(dmi))
        / max(np.linalg.norm(new_dmi), np.linalg.norm(dmi))
        > 1e-8
    ):
        parameter.dmi = new_dmi
        assert np.allclose(parameter.dmi, new_dmi)


# Legacy tests
class TestExchangeParameter:
    parameter = ExchangeParameter([[1, 5, 2], [5, 8, 4], [2, 6, 3]])

    def test_init(self):
        parameter = ExchangeParameter(iso=23)
        assert parameter.iso == 23
        assert (parameter.dmi == np.zeros(3, dtype=float)).all()
        assert (parameter.aniso == np.zeros((3, 3), dtype=float)).all()
        parameter.dmi = (1, 1, 1)
        assert parameter.iso == 23
        assert (parameter.dmi == np.ones(3, dtype=float)).all()
        assert (parameter.aniso == np.zeros((3, 3), dtype=float)).all()
        parameter.aniso = [[1, 1, 0], [1, -0.5, 0], [0, 0, -0.5]]
        assert parameter.iso == 23
        assert (parameter.dmi == np.ones(3, dtype=float)).all()
        assert (
            parameter.aniso == np.array([[1, 1, 0], [1, -0.5, 0], [0, 0, -0.5]])
        ).all()
        parameter = ExchangeParameter(iso=23, matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        assert parameter.iso == 1

    def test_matrix(self):
        parameter = ExchangeParameter()
        with pytest.raises(ValueError):
            parameter.matrix = 1
        parameter.matrix = [[1, 2, 0], [1, 1, 0], [0, 0, 1]]
        assert (parameter == np.array([[1, 2, 0], [1, 1, 0], [0, 0, 1]])).all()
        assert (parameter.matrix == np.array([[1, 2, 0], [1, 1, 0], [0, 0, 1]])).all()
        assert (self.parameter == np.array([[1, 5, 2], [5, 8, 4], [2, 6, 3]])).all()

    def test_symm_assym_matrix(self):
        parameter = ExchangeParameter(matrix=[[1, 2, 0], [1, 1, 0], [0, 0, 1]])
        assert (
            parameter.symm_matrix == np.array([[1, 1.5, 0], [1.5, 1, 0], [0, 0, 1]])
        ).all()

        assert (
            parameter.asymm_matrix == np.array([[0, 0.5, 0], [-0.5, 0, 0], [0, 0, 0]])
        ).all()

        assert (
            self.parameter.symm_matrix == np.array([[1, 5, 2], [5, 8, 5], [2, 5, 3]])
        ).all()

        assert (
            self.parameter.asymm_matrix == np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]])
        ).all()

    def test_iso(self):
        parameter = ExchangeParameter()
        parameter.iso = 23
        assert parameter.iso == 23
        parameter.iso = 0
        assert parameter.iso == 0
        assert self.parameter.iso == 4

    def test_iso_matrix(self):
        parameter = ExchangeParameter()
        parameter.iso = 23
        assert (
            parameter.iso_matrix == np.array([[23, 0, 0], [0, 23, 0], [0, 0, 23]])
        ).all()
        parameter.iso = 0
        assert (
            parameter.iso_matrix == np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
        ).all()
        assert (
            self.parameter.iso_matrix == np.array([[4, 0, 0], [0, 4, 0], [0, 0, 4]])
        ).all()

    def test_aniso(self):
        parameter = ExchangeParameter()
        parameter.aniso = [[1, 1, 0], [1, -0.5, 0], [0, 0, -0.5]]
        assert (
            parameter.aniso == np.array([[1, 1, 0], [1, -0.5, 0], [0, 0, -0.5]])
        ).all()
        parameter.iso = 23
        assert (
            parameter.aniso == np.array([[1, 1, 0], [1, -0.5, 0], [0, 0, -0.5]])
        ).all()
        parameter.aniso = np.zeros((3, 3))
        assert (parameter.aniso == np.zeros((3, 3))).all()
        with pytest.raises(ValueError):
            parameter.aniso = 1

        assert (
            self.parameter.aniso == np.array([[-3, 5, 2], [5, 4, 5], [2, 5, -1]])
        ).all()

    def test_aniso_diagonal(self):
        assert (self.parameter.aniso_diagonal == np.array([-3, 4, -1])).all()

    def test_aniso_diagonal_matrix(self):
        assert (
            self.parameter.aniso_diagonal_matrix
            == np.array([[-3, 0, 0], [0, 4, 0], [0, 0, -1]])
        ).all()

    def test_dmi(self):
        parameter = ExchangeParameter()
        parameter.dmi = (1, 2, 3)
        assert (parameter.dmi == np.array([1, 2, 3])).all()
        parameter.iso = 23
        parameter.aniso = [[1, 1, 0], [1, -0.5, 0], [0, 0, -0.5]]
        assert (parameter.dmi == np.array([1, 2, 3])).all()
        parameter.dmi = [0, 0, 0]
        assert (parameter.dmi == np.zeros(3)).all()
        with pytest.raises(ValueError):
            parameter.dmi = 1
        assert (self.parameter.dmi == np.array([-1, 0, 0])).all()

    def test_dmi_matrix(self):
        parameter = ExchangeParameter()
        parameter.dmi = (1, 2, 3)
        assert (
            parameter.dmi_matrix == np.array([[0, 3, -2], [-3, 0, 1], [2, -1, 0]])
        ).all()
        parameter.iso = 23
        parameter.aniso = [[1, 1, 0], [1, -0.5, 0], [0, 0, -0.5]]
        assert (
            parameter.dmi_matrix == np.array([[0, 3, -2], [-3, 0, 1], [2, -1, 0]])
        ).all()
        parameter.dmi = [0, 0, 0]
        assert (
            parameter.dmi_matrix == np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
        ).all()
        assert (
            self.parameter.dmi_matrix == np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]])
        ).all()

    def test_dmi_module(self):
        parameter = ExchangeParameter()
        parameter.dmi = (1, 2, 3)
        assert parameter.dmi_module == sqrt(14)
        parameter.iso = 23
        parameter.aniso = [[1, 1, 0], [1, -0.5, 0], [0, 0, -0.5]]
        assert parameter.dmi_module == sqrt(14)
        parameter.dmi = [0, 0, 0]
        assert parameter.dmi_module == 0
        assert self.parameter.dmi_module == 1

    def test_add(self):
        parameter1 = ExchangeParameter(iso=1)
        parameter2 = ExchangeParameter(iso=2)
        parameter3 = ExchangeParameter(iso=3)
        with pytest.raises(TypeError):
            parameter = parameter1 + 1
        with pytest.raises(TypeError):
            parameter = 1 + parameter1
        parameter = parameter1 + parameter3
        assert parameter.iso == 4
        parameter += parameter2
        assert parameter.iso == 6

    def test_sub(self):
        parameter1 = ExchangeParameter(iso=1)
        parameter2 = ExchangeParameter(iso=2)
        parameter3 = ExchangeParameter(iso=3)
        with pytest.raises(TypeError):
            parameter = parameter1 - 1
        with pytest.raises(TypeError):
            parameter = 1 - parameter1
        parameter = parameter1 - parameter3
        assert parameter.iso == -2
        parameter -= parameter2
        assert parameter.iso == -4

    def test_mul(self):
        parameter1 = ExchangeParameter(iso=1)
        parameter = parameter1 * 5
        assert round(parameter.iso, 10) == 5
        parameter = 1.5 * parameter
        assert round(parameter.iso, 10) == 7.5
        parameter = 5 * parameter1
        assert round(parameter.iso, 10) == 5

    def test_truediv(self):
        parameter1 = ExchangeParameter(iso=1)
        parameter = parameter1 / 5
        assert round(parameter.iso, 10) == 0.2
        parameter = parameter / 2
        assert round(parameter.iso, 10) == 0.1
        with pytest.raises(TypeError):
            parameter = 5 / parameter1

    def test_matmul(self):
        parameter = ExchangeParameter(iso=3, dmi=[1, 0, 0])
        assert (np.array([3, 3, 11]) == [1, 2, 3] @ parameter).all()
        assert (np.array([3, 9, 7]) == parameter @ [1, 2, 3]).all()
        assert (np.array([3, 3, 11]) == np.array([1, 2, 3]) @ parameter).all()
        assert (np.array([3, 9, 7]) == parameter @ np.array([1, 2, 3])).all()
        assert isinstance(parameter @ [1, 2, 3], np.ndarray)
        assert isinstance([1, 2, 3] @ parameter, np.ndarray)
        assert isinstance(parameter @ np.array([1, 2, 3]), np.ndarray)
        assert isinstance(np.array([1, 2, 3]) @ parameter, np.ndarray)

    def test_floordiv(self):
        parameter1 = ExchangeParameter(iso=23.4345)
        parameter = parameter1 // 5
        assert round(parameter.iso, 10) == 4
        parameter = parameter // 2
        assert round(parameter.iso, 10) == 2
        with pytest.raises(TypeError):
            parameter = 5 // parameter1

    def test_mod(self):
        parameter1 = ExchangeParameter(iso=23.4345)
        parameter = parameter1 % 5
        assert round(parameter.iso, 10) == 3.4345
        parameter = parameter % 2
        assert round(parameter.iso, 10) == 1.4345
        with pytest.raises(TypeError):
            parameter = 5 % parameter1

    def test_neg_pos(self):
        parameter1 = ExchangeParameter(iso=23.4345)
        assert round((-parameter1).iso, 10) == -23.4345
        assert round((+parameter1).iso, 10) == 23.4345
        assert round((--parameter1).iso, 10) == 23.4345
        assert parameter1 == --parameter1

    def test_eq(self):
        parameter1 = ExchangeParameter(iso=1)
        parameter2 = ExchangeParameter(iso=2)
        parameter3 = ExchangeParameter(iso=3)
        assert parameter1 != parameter2
        assert parameter1 == parameter2 / 2
        assert parameter2 != parameter3
        assert parameter3 == parameter1 + parameter2

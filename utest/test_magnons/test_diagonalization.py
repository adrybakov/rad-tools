import pytest
import numpy as np

from radtools.magnons.diagonalization import solve_via_colpa, ColpaFailed


@pytest.mark.parametrize(
    "h",
    [
        ([[1, 0], [0, 1]]),
        (np.eye(20)),
        (np.diag(np.concatenate((np.linspace(1, 10, 10), np.linspace(10, 1, 10))))),
        ([[2, 1], [1, 2]]),
        ([[2, 0.1, 1, 0], [0.1, 2, 0, 1], [1, 0, 2, 0.1], [0, 1, 0.1, 2]]),
    ],
)
def test_solve_via_colpa(h):
    N = len(h) // 2
    E, G = solve_via_colpa(h)
    assert np.allclose(
        np.diag(E), np.linalg.inv(np.conjugate(G).T) @ h @ np.linalg.inv(G), rtol=1e-5
    )
    assert np.allclose(E[:N], E[: N - 1 : -1], rtol=1e-5)
    assert (E > 0).all()
    assert np.allclose(
        np.sum(G[:N, :N] ** 2, axis=1) - np.sum(G[:N, N:] ** 2, axis=1),
        np.ones(N),
        rtol=1e-5,
    )
    assert np.allclose(
        np.sum(G[N:, :N] ** 2, axis=1) - np.sum(G[N:, N:] ** 2, axis=1),
        -np.ones(N),
        rtol=1e-5,
    )


@pytest.mark.parametrize("h", [[[1, 0], [0, 0]], [[1, -1], [-1, -1]]])
def test_fail_via_colpa(h):
    with pytest.raises(ColpaFailed):
        solve_via_colpa(h)

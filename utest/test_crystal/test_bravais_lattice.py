import pytest

from rad_tools.crystal.bravais_lattice import *

lattice_names = lattice_example()
lattices = [lattice_example(i) for i in lattice_names]
n_edges = [
    12,
    36,
    24,
    12,
    28,
    36,
    12,
    28,
    36,
    24,
    36,
    18,
    18,
    36,
    24,
    18,
    36,
    28,
    28,
    24,
    36,
    36,
    28,
    36,
    28,
]

n_vertices = [
    8,
    24,
    14,
    8,
    18,
    24,
    8,
    18,
    24,
    14,
    24,
    12,
    12,
    24,
    14,
    12,
    24,
    18,
    18,
    14,
    24,
    24,
    18,
    24,
    18,
]


@pytest.mark.parametrize(
    "lattice, variation", list(zip(lattices, lattice_names)), ids=lattice_example()
)
def test_variants(lattice, variation):
    assert lattice.variation == variation


@pytest.mark.parametrize(
    "lattice, n_edge",
    list(zip(lattices, n_edges)),
    ids=lattice_example(),
)
def test_edges(lattice, n_edge):
    edges, vertices = lattice.voronoi_cell(reciprocal=True)
    assert len(edges) == n_edge


@pytest.mark.parametrize(
    "lattice, n_vertex",
    list(zip(lattices, n_vertices)),
    ids=lattice_example(),
)
def test_vertices(lattice, n_vertex):
    edges, vertices = lattice.voronoi_cell(reciprocal=True)
    assert len(vertices) == n_vertex

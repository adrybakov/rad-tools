import pytest

from rad_tools.crystal.bravais_lattice import *


@pytest.mark.parametrize(
    "lattice, variation",
    [
        (cub, "CUB"),
        (fcc, "FCC"),
        (bcc, "BCC"),
        (tet, "TET"),
        (bct1, "BCT1"),
        (bct2, "BCT2"),
        (orc, "ORC"),
        (orcf1, "ORCF1"),
        (orcf2, "ORCF2"),
        (orcf3, "ORCF3"),
        (orci, "ORCI"),
        (orcc, "ORCC"),
        (hex, "HEX"),
        (rhl1, "RHL1"),
        (rhl2, "RHL2"),
        (mcl, "MCL"),
        (mclc1, "MCLC1"),
        (mclc2, "MCLC2"),
        (mclc3, "MCLC3"),
        (mclc4, "MCLC4"),
        (mclc5, "MCLC5"),
        (tri1a, "TRI1a"),
        (tri1b, "TRI1b"),
        (tri2a, "TRI2a"),
        (tri2b, "TRI2b"),
    ],
)
def test_variants(lattice, variation):
    assert lattice.variation == variation

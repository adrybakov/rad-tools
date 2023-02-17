from math import sqrt, pi

import pytest
import numpy as np

from rad_tools.kpoints.kpoints import KPoints


class TestKPoints:

    def test_labels(self):
        tmp = KPoints()
        tmp.hex()
        assert tmp.labels == ["$\\Gamma$", "$M$", "$K$", "$\\Gamma$",
                              "$A$", "$L$", "$H$", "$A$$\\vert$$L$",
                              "$M$$\\vert$$K$", "$H$"]
        assert (tmp.labels_coordinates() - np.array([0,
                                                     0.5,
                                                     0.8726779962499649,
                                                     1.3440825170409965,
                                                     1.8440825170409965,
                                                     2.3440825170409965,
                                                     2.7167605132909616,
                                                     3.1881650340819934,
                                                     3.6881650340819934,
                                                     4.188165034081994]) == 0).all()

    def test_coordinatets(self):
        tmp = KPoints()
        tmp.hex()
        assert tmp.coordinates(n=1).shape == (18, 3)
        assert tmp.coordinates(n=10).shape == (99, 3)
        assert tmp.coordinates(n=1, linear=True).shape == (18,)
        assert (tmp.coordinates(n=1) - np.array([[0., 0., 0.],
                                                 [1/2, 0., 0.],
                                                 [1/2, 0., 0.],
                                                 [1/3, 1/3, 0.],
                                                 [1/3, 1/3, 0.],
                                                 [0., 0., 0.],
                                                 [0., 0., 0.],
                                                 [0., 0., 1/2],
                                                 [0., 0., 1/2],
                                                 [1/2, 0., 1/2],
                                                 [1/2, 0., 1/2],
                                                 [1/3, 1/3, 1/2],
                                                 [1/3, 1/3, 1/2],
                                                 [0., 0., 1/2],
                                                 [1/2, 0., 1/2],
                                                 [1/2, 0., 0.],
                                                 [1/3, 1/3, 0.],
                                                 [1/3, 1/3, 1/2]]) == 0).all()
        assert (tmp.coordinates(n=2) - np.array([[0., 0., 0.],
                                                 [1/4, 0., 0.],
                                                 [1/2, 0., 0.],
                                                 [1/2, 0., 0.],
                                                 [(1/3+1/2)/2, 1/6, 0.],
                                                 [1/3, 1/3, 0.],
                                                 [1/3, 1/3, 0.],
                                                 [1/6, 1/6, 0.],
                                                 [0., 0., 0.],
                                                 [0., 0., 0.],
                                                 [0., 0., 1/4],
                                                 [0., 0., 1/2],
                                                 [0., 0., 1/2],
                                                 [1/4, 0., 1/2],
                                                 [1/2, 0., 1/2],
                                                 [1/2, 0., 1/2],
                                                 [(1/3+1/2)/2, 1/6, 1/2],
                                                 [1/3, 1/3, 1/2],
                                                 [1/3, 1/3, 1/2],
                                                 [1/6, 1/6, 1/2],
                                                 [0., 0., 1/2],
                                                 [1/2, 0., 1/2],
                                                 [1/2, 0., 1/4],
                                                 [1/2, 0., 0.],
                                                 [1/3, 1/3, 0.],
                                                 [1/3, 1/3, 1/4],
                                                 [1/3, 1/3, 1/2]]) == 0).all()

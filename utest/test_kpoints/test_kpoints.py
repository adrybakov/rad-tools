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

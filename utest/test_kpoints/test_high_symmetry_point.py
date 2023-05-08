import numpy as np
import pytest

from rad_tools.kpoints import HighSymmetryPoints


class TestHighSymmetryPoints:
    def test_init(self):
        tmp = HighSymmetryPoints(
            points={"R": [0, 0, 0], "G": [1, 2, 3]}, path=[["R", "G"]]
        )
        assert tmp.path == [["R", "G"]]
        assert tmp._points == {"R": [0, 0, 0], "G": [1, 2, 3]}

    def test_path(self):
        tmp = HighSymmetryPoints()
        tmp.path = None
        assert tmp.path == []
        with pytest.raises(ValueError):
            tmp.path = "G"

        with pytest.raises(ValueError):
            tmp.path = ["G"]

        with pytest.raises(ValueError):
            tmp.path = [["G"]]

        with pytest.raises(ValueError):
            tmp.path = "A-K-|M|"

        with pytest.raises(ValueError):
            tmp.path = "A|-K-M|"

        tmp = HighSymmetryPoints().hex()
        assert tmp.path == [
            ["Gamma", "M", "K", "Gamma", "A", "L", "H", "A"],
            ["L", "M"],
            ["K", "H"],
        ]
        tmp.path = "A-K-M|H-L"
        assert tmp.path == [["A", "K", "M"], ["H", "L"]]

        tmp.path = "A-K-M|"
        assert tmp.path == [["A", "K", "M"]]

        tmp.path = "|A-K-M"
        assert tmp.path == [["A", "K", "M"]]

        tmp.path = "A"
        assert tmp.path == [["A"]]

        tmp.path = [["A", "K", "M", "H"]]
        assert tmp.path == [["A", "K", "M", "H"]]

        tmp.path = []
        assert tmp.path == []

        tmp.path = [[]]
        assert tmp.path == []

        tmp.path = ["A", "K", "M", "H"]
        assert tmp.path == [["A", "K", "M", "H"]]

    def test_add_kpoints(self):
        tmp = HighSymmetryPoints()
        tmp.add_kpoint("F", [4, 2, 1])
        assert tmp.points["F"].tolist() == [4, 2, 1]
        assert tmp._PLOT_LITERALS["F"] == "$F$"

        tmp.add_kpoint("R13", [4, 5, 6])
        assert tmp._PLOT_LITERALS["R13"] == "R13"
        assert tmp.points["R13"].tolist() == [4, 5, 6]

        tmp = HighSymmetryPoints()
        tmp.add_kpoint("R13", [4, 7, 6], plot_name="$R_{13}$")
        assert tmp._PLOT_LITERALS["R13"] == "$R_{13}$"
        assert tmp.points["R13"].tolist() == [4, 7, 6]

        tmp = HighSymmetryPoints()
        tmp.add_kpoint("R", [4, 7, 6], plot_name="Rabbit")
        assert tmp._PLOT_LITERALS["R"] == "Rabbit"
        assert tmp.points["R"].tolist() == [4, 7, 6]

        tmp = HighSymmetryPoints()
        tmp.add_kpoint("R1", [4, 7, 6], plot_name="Rabbit1")
        tmp.add_kpoint("R2", [3, 7, 1], plot_name="Rabbit2")
        print(type(tmp.points["R1"]))
        assert tmp._PLOT_LITERALS["R1"] == "Rabbit1"
        assert tmp.points["R1"].tolist() == [4, 7, 6]
        assert tmp._PLOT_LITERALS["R2"] == "Rabbit2"
        assert tmp.points["R2"].tolist() == [3, 7, 1]
        assert (tmp.points["R1"] - tmp.points["R2"] == np.array([1, 0, 5])).all()

# RAD-tools - Sandbox (mainly condense matter plotting).
# Copyright (C) 2022-2024  Andrey Rybakov
#
# e-mail: anry@uv.es, web: rad-tools.org
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os

import pytest

from radtools.io.internal import *


class ReadTemplate:
    template = load_template(
        os.path.join("utest", "test_io", "resources", "exchange.out")
    )

    def test_empty_filename(self):
        with pytest.raises(TypeError):
            template = load_template(None)

    def test_wrong_filename(self):
        with pytest.raises(FileNotFoundError):
            template = load_template("Ah, music. " + "A magic beyond all we do here!")

    def count_entries(self, dictionary):
        i = 0
        for name in dictionary:
            i += len(dictionary[name])
        return i

    def test_read_neighbors(self):
        assert self.count_entries(self.template.names) == 28
        assert self.template.latex_names == {
            "J1": "$J_1$",
            "J1(1)": "$J_1(1)$",
            "J1(2)": "$J_1(2)$",
            "J2": "J2",
            "J2(1)": "$J_2(1)$",
            "J2(2)": "$J_2(2)$",
            "J3": "$J_3$",
        }

    def test_template(self):
        assert self.template.names == {
            "J1": [
                ("Cr1", "Cr1", (1, 0, 0)),
                ("Cr1", "Cr1", (-1, 0, 0)),
                ("Cr2", "Cr2", (1, 0, 0)),
                ("Cr2", "Cr2", (-1, 0, 0)),
            ],
            "J1(1)": [("Cr1", "Cr1", (1, 0, 0)), ("Cr2", "Cr2", (-1, 0, 0))],
            "J1(2)": [("Cr1", "Cr1", (-1, 0, 0)), ("Cr2", "Cr2", (1, 0, 0))],
            "J2": [
                ("Cr1", "Cr2", (0, 0, 0)),
                ("Cr1", "Cr2", (1, 0, 0)),
                ("Cr1", "Cr2", (0, -1, 0)),
                ("Cr1", "Cr2", (1, -1, 0)),
                ("Cr2", "Cr1", (0, 0, 0)),
                ("Cr2", "Cr1", (-1, 0, 0)),
                ("Cr2", "Cr1", (0, 1, 0)),
                ("Cr2", "Cr1", (-1, 1, 0)),
            ],
            "J2(1)": [
                ("Cr1", "Cr2", (1, 0, 0)),
                ("Cr1", "Cr2", (1, -1, 0)),
                ("Cr2", "Cr1", (-1, 0, 0)),
                ("Cr2", "Cr1", (-1, 1, 0)),
            ],
            "J2(2)": [
                ("Cr1", "Cr2", (0, 0, 0)),
                ("Cr1", "Cr2", (0, -1, 0)),
                ("Cr2", "Cr1", (0, 0, 0)),
                ("Cr2", "Cr1", (0, 1, 0)),
            ],
            "J3": [
                ("Cr1", "Cr1", (0, 1, 0)),
                ("Cr1", "Cr1", (0, -1, 0)),
                ("Cr2", "Cr2", (0, 1, 0)),
                ("Cr2", "Cr2", (0, -1, 0)),
            ],
        }
        assert len(self.template.get_list()) == 28

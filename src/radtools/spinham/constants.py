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

__all__ = [
    "PREDEFINED_NOTATIONS",
    "TXT_FLAGS",
]

PREDEFINED_NOTATIONS = {
    "standard": (True, False, -1),
    "tb2j": (True, True, -1),
    "spinw": (True, False, 1),
    "vampire": (True, True, -0.5),
}

TXT_FLAGS = {
    "cell": "Primitive unit cell:",
    "atoms": "Atoms:",
    "magmoms": "Magnetic moments of magnetic atoms:",
    "spins": "Spin of magnetic atoms:",
    "spinham": "Spin Hamiltonian:",
    "iso": "Isotropic exchange:",
    "matrix": "Full exchange matrix:",
    "dmi": "DMI vector:",
    "dmis": "DMI vectors:",
    "dmi_module": "DMI module:",
    "dmi_relative": "|DMI| / |J_iso|:",
    "aniso": "Symmetric anisotropic exchange:",
    "spinmodel": "Spin model:",
    "kpath": "K-path:",
    "kpoints": "Points (in relative coordinates):",
    "klabels": "Labels:",
    "dispersion": "Magnon dispersion:",
    "separate": "Data are in the file:",
    "notation": "Notation of Hamiltonian:",
    "bonds": "Bonds:",
}

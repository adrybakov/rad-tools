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

__all__ = ["ColpaFailed", "NotationError"]


class ColpaFailed(Exception):
    r"""
    Raised when Diagonalization via Colpa fails.
    """

    def __init__(self):
        self.message = "Diagonalization via Colpa failed."

    def __str__(self):
        return self.message


class NotationError(ValueError):
    r"""
    Raised when the notation (or individual property) is not defined.

    Gives a summary of the notation (or individual property) and how to set it.

    Parameters
    ----------
    name : str
        Name of the corresponding attribute.

    """

    def __init__(self, name):
        self.message = (
            f"\n\nNotation`s interpretation is not set for the property {name}.\n"
            + f"Set the notation first:\n"
            + f"    SpinHamiltonian.{name} = True  "
            + f"or  SpinHamiltonian.{name} = False\n\n"
            + f"Note: When the attribute is set for the first time it sets the interpretation, "
            + "afterwards it change the notation.\n\n"
            + f"If you want to set the interpretation again, use \n"
            + f"    SpinHamiltonian.set_interpretation({name} = True)"
            + "\nor\n"
            + f"    SpinHamiltonian.set_interpretation({name} = False)\n"
        )

    def __str__(self):
        return self.message

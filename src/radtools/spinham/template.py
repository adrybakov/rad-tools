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

r"""
Exchange template.

Write a tutorial with docstring here.
"""
__all__ = ["ExchangeTemplate"]


class ExchangeTemplate:
    r"""
    Store a template for the :py:class:`.SpinHamiltonian`.

    In addition stores the technical details for plotting,
    orbital decomposition, etc.

    Attributes
    ----------
    names : dict
        Dictionary of neighbors from the template file.

        .. code-block:: python

            {name : [(atom1, atom2, R), ...], ...}

    latex_names : dict
        The dictionary of Latex version of names from `names`.

        .. code-block:: python

            {name : latex_name, ...}
    """

    def __init__(self) -> None:
        self.names = {}
        self.latex_names = {}

    def get_list(self):
        r"""
        Translate bonds present in the template into the list.

        Returns
        -------
        template_list : list
            List with the bond specifications:

            .. code-block:: python

                [(atom_1, atom_2, R), ...]
        """
        template_list = []
        for name in self.names:
            for atom1, atom2, R in self.names[name]:
                template_list.append((atom1, atom2, R))
        return template_list

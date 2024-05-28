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
Script interface to the radtools package.

The behaviour of the command line interface is defined by the
functions in this module.

The functions are called with the same names as the scripts,
but the prefix "rad-" is removed and "-" are substituted by "_".
Function`s arguments directly correspond to the full names of the
arguments of the script (i.e. the argument :ref:`rad-extract-tb2j_input-filename`
of the script :ref:`rad-extract-tb2j` is passed to the function :py:func:`.extract_tb2j` as the
argument ``input_filename``).

Full documentation on the behaviour is available in the
:ref:`scripts-guide`.


.. admonition:: Example

    Identification of Wannier centres from the file "seedname_centres.xyz"
    with increased span of 0.2, saving the result in the file
    "identified_centres" of the current directory:

    .. code-block:: bash

        rad-identify-wannier-centres.py seedname_centres.xyz -s 0.2 -on identified_centres

    The same result could be achieved by calling the function
    :py:func:`.identify_wannier_centres`:

    .. code-block:: python

        from radtools import identify_wannier_centres
        identify_wannier_centres("seedname_centres.xyz",
            span = 0.2,
            output_name="identified_centres")
"""

from .extract_tb2j import manager as extract_tb2j
from .identify_wannier_centres import manager as identify_wannier_centres
from .make_template import manager as make_template
from .plot_dos import manager as plot_dos
from .plot_fatbands import manager as plot_fatbands
from .plot_tb2j import manager as plot_tb2j
from .plot_tb2j_magnons import manager as plot_tb2j_magnons

__all__ = [
    "identify_wannier_centres",
    "plot_tb2j",
    "extract_tb2j",
    "plot_dos",
    "make_template",
    "plot_tb2j_magnons",
    "plot_fatbands",
]

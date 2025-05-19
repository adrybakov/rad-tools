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
from radtools._redirect import *

__all__ = [
    "CUB",
    "FCC",
    "BCC",
    "TET",
    "BCT",
    "ORC",
    "ORCF",
    "ORCI",
    "ORCC",
    "HEX",
    "RHL",
    "MCL",
    "MCLC",
    "TRI",
]


def CUB(a, return_cell=False):
    redirect_to_wulfric()


def FCC(a, return_cell=False):
    redirect_to_wulfric()


def BCC(a, return_cell=False):
    redirect_to_wulfric()


def TET(a, c, return_cell=False):
    redirect_to_wulfric()


def BCT(a, c, return_cell=False):
    redirect_to_wulfric()


def ORC(a, b, c, return_cell=False):
    redirect_to_wulfric()


def ORCF(a, b, c, return_cell=False):
    redirect_to_wulfric()


def ORCI(a, b, c, return_cell=False):
    redirect_to_wulfric()


def ORCC(a, b, c, return_cell=False):
    redirect_to_wulfric()


def HEX(a, c, return_cell=False):
    redirect_to_wulfric()


def RHL(a, alpha, return_cell=False):
    redirect_to_wulfric()


def MCL(a, b, c, alpha, return_cell=False):
    redirect_to_wulfric()


def MCLC(a, b, c, alpha, return_cell=False):
    redirect_to_wulfric()


def TRI(a, b, c, alpha, beta, gamma, reciprocal=False, return_cell=False):
    redirect_to_wulfric()

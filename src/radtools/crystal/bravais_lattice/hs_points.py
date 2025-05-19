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
    "CUB_hs_points",
    "FCC_hs_points",
    "BCC_hs_points",
    "TET_hs_points",
    "BCT_hs_points",
    "ORC_hs_points",
    "ORCF_hs_points",
    "ORCI_hs_points",
    "ORCC_hs_points",
    "HEX_hs_points",
    "RHL_hs_points",
    "MCL_hs_points",
    "MCLC_hs_points",
    "TRI_hs_points",
]


def CUB_hs_points():
    redirect_to_wulfric()


def FCC_hs_points():
    redirect_to_wulfric()


def BCC_hs_points():
    redirect_to_wulfric()


def TET_hs_points():
    redirect_to_wulfric()


def BCT_hs_points(variation, conv_a, conv_c):
    redirect_to_wulfric()


def ORC_hs_points():
    redirect_to_wulfric()


def ORCF_hs_points(variation, conv_a, conv_b, conv_c):
    redirect_to_wulfric()


def ORCI_hs_points(conv_a, conv_b, conv_c):
    redirect_to_wulfric()


def ORCC_hs_points(conv_a, conv_b):
    redirect_to_wulfric()


def HEX_hs_points():
    redirect_to_wulfric()


def RHL_hs_points(variation, conv_alpha):
    redirect_to_wulfric()


def MCL_hs_points(conv_b, conv_c, conv_alpha):
    redirect_to_wulfric()


def MCLC_hs_points(variation, conv_a, conv_b, conv_c, conv_alpha):
    redirect_to_wulfric()


def TRI_hs_points(variation):
    redirect_to_wulfric()

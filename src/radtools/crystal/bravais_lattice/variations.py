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
    "BCT_variation",
    "ORCF_variation",
    "RHL_variation",
    "MCLC_variation",
    "TRI_variation",
]


def BCT_variation(conv_a, conv_c):
    redirect_to_wulfric()


def ORCF_variation(conv_a, conv_b, conv_c, eps):
    redirect_to_wulfric()


def RHL_variation(conv_alpha, eps):
    redirect_to_wulfric()


def MCLC_variation(conv_a, conv_b, conv_c, conv_alpha, k_gamma, eps):
    redirect_to_wulfric()


def TRI_variation(k_alpha, k_beta, k_gamma, eps):
    redirect_to_wulfric()

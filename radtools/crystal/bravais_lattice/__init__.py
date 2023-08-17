from .constructor import *
from .examples import *
from .hs_points import *
from .standardize import *
from .variations import *

__all__ = []
__all__.extend(examples.__all__)
__all__.extend(standardize.__all__)
__all__.extend(constructor.__all__)
__all__.extend(hs_points.__all__)
__all__.extend(variations.__all__)


r"""
Bravais lattices.
"""

# High symmetry k points getters


# # 5
# class BCT(Lattice):
#     def __init__(self) -> None:
#         self._PLOT_NAMES["S"] = "$\\Sigma$"
#         self._PLOT_NAMES["S1"] = "$\\Sigma_1$"


# # 14
# class TRI(Lattice):
#     def __init__(
#         self,
#         a: float,
#         b: float,
#         c: float,
#         alpha: float,
#         beta: float,
#         gamma: float,
#         reciprocal=False,
#     ) -> None:
#         if not reciprocal:
#             tmp = sorted([(a, alpha), (b, beta), (c, gamma)], key=lambda x: x[0])
#             a = tmp[0][0]
#             alpha = tmp[0][1]
#             b = tmp[1][0]
#             beta = tmp[1][1]
#             c = tmp[2][0]
#             gamma = tmp[2][1]
#         super().__init__(a, b, c, alpha, beta, gamma)
#         if reciprocal:
#             self.cell = self.reciprocal_cell
#         self.conv_a = a
#         self.conv_b = b
#         self.conv_c = c
#         self.conv_alpha = alpha
#         self.conv_beta = beta
#         self.conv_gamma = gamma
#         if (
#             self.k_alpha < self.k_gamma < self.k_beta
#             or self.k_beta < self.k_gamma < self.k_alpha
#         ):
#             raise RuntimeError("k_gamma is not minimal nor maximal.")
#         self.conv_cell = self.cell

from math import cos, sin

from radtools.constants import TORADIANS
from radtools.numerical import compare_numerically

__all__ = [
    "BCT_variation",
    "ORCF_variation",
    "RHL_variation",
    "MCLC_variation",
    "TRI_variation",
]


def BCT_variation(conv_a: float, conv_c: float):
    r"""
    Two variations of the BCT lattice.

    Condition :math:`a \ne c` is assumed.

    :math:`\text{BCT}_1: c < a` and :math:`\text{BCT}_2: c > a`

    Parameters
    ----------
    conv_a : float
        Length of the :math:`a_1 == a_2` vector of the conventional cell.
    conv_c : float
        Length of the :math:`a_3` vector of the conventional cell.

    Returns
    -------
    variation : str
        Variation of the lattice. "BCT1" or "BCT2".

    Raises
    ------
    ValueError
        If :math:`a == c`.
    """
    if conv_a > conv_c:
        return "BCT1"
    elif conv_a < conv_c:
        return "BCT2"
    else:
        raise ValueError("a == c")


def ORCF_variation(conv_a: float, conv_b: float, conv_c: float, eps: float):
    r"""
    Three variations of the ORCF lattice.

    Ordering :math:`a < b < c` is assumed.

    :math:`\text{ORCF}_1: \dfrac{1}{a^2} > \dfrac{1}{b^2} + \dfrac{1}{c^2}`,
    :math:`\text{ORCF}_2: \dfrac{1}{a^2} < \dfrac{1}{b^2} + \dfrac{1}{c^2}`,
    :math:`\text{ORCF}_3: \dfrac{1}{a^2} = \dfrac{1}{b^2} + \dfrac{1}{c^2}`,

    Parameters
    ----------
    conv_a : float
        Length of the :math:`a_1` vector of the conventional cell.
    conv_b : float
        Length of the :math:`a_2` vector of the conventional cell.
    conv_c : float
        Length of the :math:`a_3` vector of the conventional cell.
    eps : float
        Tolerance for numerical comparison.

    Returns
    -------
    variation : str
        Variation of the lattice. "ORCF1", "ORCF2" or "ORCF3".

    Raises
    ------
    ValueError
        If :math:`a < b < c` is not satisfied.
    """
    if compare_numerically(conv_a, ">=", conv_b, eps) or compare_numerically(
        conv_b, ">=", conv_c, eps
    ):
        raise ValueError(f"a < b < c is not satisfied with {eps} tolerance.")

    expression = 1 / conv_a**2 - 1 / conv_b**2 - 1 / conv_c**2
    if compare_numerically(expression, "==", 0, eps):
        return "ORCF3"
    elif compare_numerically(expression, ">", 0, eps):
        return "ORCF1"
    elif compare_numerically(expression, "<", 0, eps):
        return "ORCF2"


def RHL_variation(conv_alpha: float, eps: float):
    r"""
    Two variations of the RHL lattice.

    Condition :math:`\alpha \ne 90^{\circ}` is assumed.

    :math:`\text{RHL}_1 \alpha < 90^{\circ}`,
    :math:`\text{RHL}_2 \alpha > 90^{\circ}`

    Parameters
    ----------
    conv_alpha : float
        Angle between vectors :math:`a_1` and :math:`a_2` of the conventional cell. In degrees.
    eps : float
        Tolerance for numerical comparison.

    Returns
    -------
    variation : str
        Variation of the lattice. Either "RHL1" or "RHL2".

    Raises
    ------
    ValueError
        If :math:`\alpha == 90^{\circ}` with given tolerance ``eps``.
    """
    if compare_numerically(conv_alpha, "<", 90, eps):
        return "RHL1"
    elif compare_numerically(conv_alpha, ">", 90, eps):
        return "RHL2"
    else:
        raise ValueError(f"alpha == 90 with {eps} tolerance.")


def MCLC_variation(
    conv_a: float,
    conv_b: float,
    conv_c: float,
    conv_alpha: float,
    k_gamma: float,
    eps: float,
):
    r"""
    Five variation of the MCLC lattice.

    Ordering :math:`a \le c` and :math:`b \le c` and :math:`\alpha < 90^{\circ}` is assumed.

    :math:`\text{MCLC}_1: k_{\gamma} > 90^{\circ}`,
    :math:`\text{MCLC}_2: k_{\gamma} = 90^{\circ}`,
    :math:`\text{MCLC}_3: k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} < 1`
    :math:`\text{MCLC}_4: k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} = 1`
    :math:`\text{MCLC}_5: k_{\gamma} < 90^{\circ}, \dfrac{b\cos(\alpha)}{c} + \dfrac{b^2\sin(\alpha)^2}{a^2} > 1`

    Parameters
    ----------
    conv_a : float
        Length of the :math:`a_1` vector of the conventional cell.
    conv_b : float
        Length of the :math:`a_2` vector of the conventional cell.
    conv_c : float
        Length of the :math:`a_3` vector of the conventional cell.
    conv_alpha : float
        Angle between vectors :math:`a_2` and :math:`a_3` of the conventional cell. In degrees.
    k_gamma : float
        Angle between reciprocal vectors :math:`b_1` and :math:`b_2`. In degrees.
    eps : float
        Tolerance for numerical comparison.

    Returns
    -------
    variation : str
        Variation of the lattice.
        Either "MCLC1", "MCLC2", "MCLC3", "MCLC4" or "MCLC5".

    Raises
    ------
    ValueError
        If :math:`\alpha > 90^{\circ}` or :math:`a > c` or :math:`b > c` with given tolerance ``eps``.
    """

    if compare_numerically(conv_alpha, ">", 90, eps) or compare_numerically(
        conv_b, ">", conv_c, eps
    ):
        raise ValueError(
            f"alpha > 90 or  or b > c with {eps} tolerance:\n"
            + f"  alpha = {conv_alpha}\n"
            + f"  b = {conv_b}\n"
            + f"  c = {conv_c}\n"
        )

    conv_alpha *= TORADIANS

    if compare_numerically(k_gamma, "==", 90, eps):
        return "MCLC2"
    elif compare_numerically(k_gamma, ">", 90, eps):
        return "MCLC1"
    elif compare_numerically(k_gamma, "<", 90, eps):
        expression = (
            conv_b * cos(conv_alpha) / conv_c
            + conv_b**2 * sin(conv_alpha) ** 2 / conv_a**2
        )
        if compare_numerically(expression, "==", 1, eps):
            return "MCLC4"
        elif compare_numerically(expression, "<", 1, eps):
            return "MCLC3"
        elif compare_numerically(expression, ">", 1, eps):
            return "MCLC5"


def TRI_variation(k_alpha: float, k_beta: float, k_gamma: float, eps: float):
    r"""
    Four variations of the TRI lattice.

    Conditions :math:`k_{\alpha} \ne 90^{\circ}` and :math:`k_{\beta} \ne 90^{\circ}` are assumed.

    :math:`\text{TRI}_{1a} k_{\alpha} > 90^{\circ}, k_{\beta} > 90^{\circ}, k_{\gamma} > 90^{\circ}, k_{\gamma} = \min(k_{\alpha}, k_{\beta}, k_{\gamma})`

    :math:`\text{TRI}_{1b} k_{\alpha} < 90^{\circ}, k_{\beta} < 90^{\circ}, k_{\gamma} < 90^{\circ}, k_{\gamma} = \max(k_{\alpha}, k_{\beta}, k_{\gamma})`

    :math:`\text{TRI}_{2a} k_{\alpha} > 90^{\circ}, k_{\beta} > 90^{\circ}, k_{\gamma} = 90^{\circ}`

    :math:`\text{TRI}_{2b} k_{\alpha} < 90^{\circ}, k_{\beta} < 90^{\circ}, k_{\gamma} = 90^{\circ}`

    Parameters
    ----------
    k_alpha : float
        Angle between reciprocal vectors :math:`b_2` and :math:`b_3`. In degrees.
    k_beta : float
        Angle between reciprocal vectors :math:`b_1` and :math:`b_3`. In degrees.
    k_gamma : float
        Angle between reciprocal vectors :math:`b_1` and :math:`b_2`. In degrees.
    eps : float
        Tolerance for numerical comparison.

    Returns
    -------
    variation : str
        Variation of the lattice.
        Either "TRI1a", "TRI1b", "TRI2a" or "TRI2b".

    Raises
    ------
    ValueError
        If :math:`k_{\alpha} == 90^{\circ}` or :math:`k_{\beta} == 90^{\circ}` with given tolerance ``eps``.
    """

    if compare_numerically(k_alpha, "==", 90, eps) or compare_numerically(
        k_beta, "==", 90, eps
    ):
        raise ValueError(f"k_alpha == 90 or k_beta == 90 with {eps} tolerance.")

    if compare_numerically(k_gamma, "==", 90, eps):
        if compare_numerically(k_alpha, ">", 90, eps) and compare_numerically(
            k_beta, ">", 90, eps
        ):
            return "TRI2a"
        elif compare_numerically(k_alpha, "<", 90, eps) and compare_numerically(
            k_beta, "<", 90, eps
        ):
            return "TRI2b"
    elif compare_numerically(min(k_gamma, k_beta, k_alpha), ">", 90, eps):
        return "TRI1a"
    elif compare_numerically(max(k_gamma, k_beta, k_alpha), "<", 90, eps):
        return "TRI1b"
    else:
        return "TRI"

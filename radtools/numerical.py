r"""
Collection of small routines and constants, 
which are used across the whole package.

It's purpose is to serve as an "other" folder.
"""

from radtools.crystal.constants import ABS_TOL, REL_TOL

__all__ = [
    "compare_numerically",
]


def compare_numerically(x, condition, y, eps=None, rtol=REL_TOL, atol=ABS_TOL):
    r"""
    Compare two numbers numerically.

    The approach is taken from [1]_:

    .. math::

        \begin{matrix}
            x < y  & x < y - \varepsilon \\
            x > y  & y < x - \varepsilon \\
            x \le y & \text{ not } (y < x - \varepsilon) \\
            x \ge y & \text{ not } (x < y - \varepsilon) \\
            x = y  & \text{ not } (x < y - \varepsilon \text{ or } y < x - \varepsilon) \\
            x \ne y & x < y - \varepsilon \text{ or } y < x - \varepsilon
        \end{matrix}
        
    Parameters
    ----------
    x : float
        First number.
    condition : str
        Condition to compare with. One of "<", ">", "<=", ">=", "==", "!=".
    y : float
        Second number.
    eps : float, optional
        Tolerance. Used for the comparison if provided. If ``None``, then computed as:

        .. code-block:: python

            eps = atol + rtol * abs(y)

    rtol : float, default 1e-04
        Relative tolerance.
    atol : float, default 1e-08
        Absolute tolerance.

    Returns
    -------
    result: bool
        Whether the condition is satisfied.

    Raises
    ------
    ValueError
        If ``condition`` is not one of "<", ">", "<=", ">=", "==", "!=".

    References
    ----------
    .. [1] Grosse-Kunstleve, R.W., Sauter, N.K. and Adams, P.D., 2004.
        Numerically stable algorithms for the computation of reduced unit cells.
        Acta Crystallographica Section A: Foundations of Crystallography,
        60(1), pp.1-6.
    """

    if eps is None:
        eps = atol + rtol * abs(y)

    if condition == "<":
        return x < y - eps
    if condition == ">":
        return y < x - eps
    if condition == "<=":
        return not y < x - eps
    if condition == ">=":
        return not x < y - eps
    if condition == "==":
        return not (x < y - eps or y < x - eps)
    if condition == "!=":
        return x < y - eps or y < x - eps

    raise ValueError(f'Condition must be one of "<", ">", "<=", ">=", "==", "!=".')

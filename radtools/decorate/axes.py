from typing import Iterable

__all__ = ["plot_hlines", "plot_vlines"]


def plot_hlines(ax, positions, **kwargs):
    r"""
    Plot horizontal lines.

    Parameters
    ----------
    ax : :matplotlib:`axes`
        Axes to plot on.
    positions : list
        List of y positions.
    **kwargs
        Keyword arguments to be passed to :meth:`matplotlib.axes.Axes.hlines`.
    """

    if not isinstance(positions, Iterable):
        positions = [positions]

    if "transfrom" in kwargs:
        kwargs.pop("transform")
    if "color" not in kwargs:
        kwargs["color"] = "grey"
    if "linewidths" not in kwargs:
        kwargs["linewidths"] = 0.5
    if "linestyles" not in kwargs:
        kwargs["linestyles"] = "dashed"

    ax.hlines(positions, 0, 1, transform=ax.get_yaxis_transform(), **kwargs)


def plot_vlines(ax, positions, **kwargs):
    r"""
    Plot vertical lines.

    Parameters
    ----------
    ax : :matplotlib:`axes`
        Axes to plot on.
    positions : list
        List of x positions.
    **kwargs
        Keyword arguments to be passed to :meth:`matplotlib.axes.Axes.vlines`.
    """

    if not isinstance(positions, Iterable):
        positions = [positions]

    if "transfrom" in kwargs:
        kwargs.pop("transform")
    if "color" not in kwargs:
        kwargs["color"] = "grey"
    if "linewidths" not in kwargs:
        kwargs["linewidths"] = 0.5
    if "linestyles" not in kwargs:
        kwargs["linestyles"] = "dashed"

    ax.vlines(positions, 0, 1, transform=ax.get_xaxis_transform(), **kwargs)

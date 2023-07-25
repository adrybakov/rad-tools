from typing import Iterable


__all__ = ["plot_hlines", "plot_vlines"]


def plot_hlines(ax, positions):
    r"""
    Plot horizontal lines.

    Parameters
    ----------
    ax : :matplotlib:`axes`
        Axes to plot on.
    positions : list
        List of y positions.
    """

    if not isinstance(positions, Iterable):
        positions = [positions]

    ax.hlines(
        positions,
        0,
        1,
        transform=ax.get_yaxis_transform(),
        color="grey",
        linewidths=0.5,
        linestyles="dashed",
    )


def plot_vlines(ax, positions):
    r"""
    Plot vertical lines.

    Parameters
    ----------
    ax : :matplotlib:`axes`
        Axes to plot on.
    positions : list
        List of x positions.
    """

    if not isinstance(positions, Iterable):
        positions = [positions]

    ax.vlines(
        positions,
        0,
        1,
        transform=ax.get_xaxis_transform(),
        color="grey",
        linewidths=0.5,
        linestyles="dashed",
    )

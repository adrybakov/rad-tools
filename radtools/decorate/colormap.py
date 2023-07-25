from matplotlib.colors import LinearSegmentedColormap, to_rgb

__all__ = ["custom_cmap"]


def custom_cmap(start_color, finish_color):
    r"""
    Prepare custom colormap. From one color to the other.

    Parameters
    ----------
    start_color : Color
        Start color.
    finish_color : Color
        Finish color.
    """

    r1, g1, b1 = to_rgb(start_color)

    r2, g2, b2 = to_rgb(finish_color)

    cdict = {
        "red": ((0, r1, r1), (1, r2, r2)),
        "green": ((0, g1, g1), (1, g2, g2)),
        "blue": ((0, b1, b1), (1, b2, b2)),
    }

    return LinearSegmentedColormap("custom_cmap", cdict, N=256)

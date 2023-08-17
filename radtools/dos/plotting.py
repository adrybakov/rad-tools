import matplotlib.pyplot as plt
import numpy as np

from radtools.decorate.axes import plot_hlines, plot_vlines
from radtools.dos.pdos import PDOS

__all__ = ["COLOURS", "plot_projected"]

COLOURS = [
    "#00FFFF",
    "#FF9720",
    "#CD00FF",
    "#FFFF2B",
    "#00B9FF",
    "#FF163D",
    "#79FF35",
    "#FF0BEA",
    "#0200FF",
]


def plot_projected(
    pdos: PDOS,
    efermi=0.0,
    output_name="pdos",
    title=None,
    xlim=None,
    ylim=None,
    relative=False,
    normalize=False,
    interactive=False,
    save_pickle=False,
    colours=COLOURS,
    total_label="default",
    axes_labels_fontsize=14,
    legend_fontsize=12,
    title_fontsize=18,
):
    r"""
    Plot PDOS.

    Parameters
    ----------
    pdos : :py:class:`.PDOS`
        PDOS for the plot.
    efermi : float, default 0
        Fermi energy.
    output_name : str, default "pdos"
        output_name for the plot file. Extension ".png" is added at the end.
    title : str, optional
        Title of the plot. Passed to the ``ax.set_title()``.
    xlim : tuple
        limits for the x (Energy) axis
    ylim : tuple
        limits for the y (PDOS) axis
    relative : bool, default False
        Relative plot style.
    normalize : bool, default False
        Whether to norma;ize relative plot style.
    interactive : bool, default False
        Whether to use interactive plotting mode.
    save_pickle : bool, default False
        Whether to save figure as a .pickle file.
        Helps for custom modification of particular figures.
    colours : list
        List of colours to be used. values are passed directly to matplotlib
    total_label : str or ``None``, default "default"
        Label for the total data. If None , then the label is not added
    axes_label_fontsize : int, default 14
        Fontsize of the axes labels.
    legend_fontsize : int, default 12
        Fontsize of the legend.
    title_fontsize : int, default 18
        Fontsize of the title.
    """

    n = len(pdos.projectors)
    pdos = pdos.squeezed()

    if relative:
        fig, ax = plt.subplots(figsize=(9, 4))
    else:
        fig, axs = plt.subplots(n, 1, figsize=(9, n * 2))
        if n == 1:
            axs = [axs]
    fig.subplots_adjust(hspace=0)

    def set_up_axis(ax, i):
        if normalize:
            ax.set_ylabel("PDOS / LDOS", fontsize=axes_labels_fontsize)
        else:
            ax.set_ylabel("DOS, states/eV", fontsize=axes_labels_fontsize)
        if i == n - 1:
            if efermi == 0:
                ax.set_xlabel("E, eV", fontsize=axes_labels_fontsize)
            else:
                ax.set_xlabel("E - E$_{Fermi}$, eV", fontsize=axes_labels_fontsize)
        else:
            ax.axes.get_xaxis().set_visible(False)
        if ylim is not None:
            ax.set_ylim(*tuple(ylim))
        if xlim is not None:
            ax.set_xlim(*tuple(xlim))
        else:
            ax.set_xlim(np.amin(pdos.energy), np.amax(pdos.energy))
        if efermi != 0:
            plot_vlines(ax, 0)

        if title is not None and (i == 0 or relative):
            ax.set_title(title, fontsize=title_fontsize)

    if normalize:
        pdos = pdos.normalized(zeros_to_none=False)

    if relative:
        set_up_axis(ax, n - 1)
        plot_hlines(ax, 0)

    for i, projector in enumerate(pdos):
        if not relative:
            ax = axs[i]
            set_up_axis(ax, i)
        if pdos.spin_pol:
            if total_label == "default":
                label_up = f"{pdos.projectors_group} (up)"
                label_down = f"{pdos.projectors_group} (down)"
            elif total_label is None:
                label_up = None
                label_down = None
            else:
                label_up = f"{total_label} (up)"
                label_down = f"{total_label} (down)"
            if relative:
                if i == 0 and total_label is not None:
                    ax.plot(
                        pdos.energy,
                        pdos.ldos[0],
                        "-",
                        lw=1,
                        color="blue",
                        alpha=0.8,
                        label=label_up,
                    )
                    ax.plot(
                        pdos.energy,
                        -pdos.ldos[1],
                        "-",
                        lw=1,
                        color="red",
                        alpha=0.8,
                        label=label_down,
                    )
                ax.fill_between(
                    pdos.energy,
                    np.sum(pdos[:i], axis=0)[0],
                    np.sum(pdos[: i + 1], axis=0)[0],
                    lw=0,
                    color=colours[i % len(colours)],
                    # alpha=0.5,
                    label=f"{projector}",
                )
                ax.fill_between(
                    pdos.energy,
                    -np.sum(pdos[:i], axis=0)[1],
                    -np.sum(pdos[: i + 1], axis=0)[1],
                    lw=0,
                    color=colours[i % len(colours)],
                    # alpha=0.5,
                )
            else:
                if total_label is not None:
                    ax.fill_between(
                        pdos.energy,
                        0,
                        pdos.ldos[0],
                        lw=0,
                        color="blue",
                        alpha=0.2,
                        label=label_up,
                    )
                    ax.fill_between(
                        pdos.energy,
                        0,
                        -pdos.ldos[1],
                        lw=0,
                        color="red",
                        alpha=0.2,
                        label=label_down,
                    )

                ax.plot(
                    pdos.energy,
                    pdos[projector][0],
                    "-",
                    lw=0.8,
                    color="blue",
                    alpha=0.8,
                    label=f"{projector} (up)",
                )
                ax.plot(
                    pdos.energy,
                    -pdos[projector][1],
                    "-",
                    lw=0.8,
                    color="red",
                    alpha=0.8,
                    label=f"{projector} (down)",
                )
        else:
            if total_label == "default":
                total_label = pdos.projectors_group
            if relative:
                if i == 0 and total_label is not None:
                    ax.plot(
                        pdos.energy,
                        pdos.ldos,
                        "-",
                        lw=1,
                        color="black",
                        alpha=0.8,
                        label=total_label,
                    )
                ax.fill_between(
                    pdos.energy,
                    np.sum(pdos[:i], axis=0),
                    np.sum(pdos[: i + 1], axis=0),
                    lw=0,
                    color=colours[i % len(colours)],
                    # alpha=0.5,
                    label=projector,
                )

            else:
                if total_label is not None:
                    ax.fill_between(
                        pdos.energy,
                        0,
                        pdos.ldos,
                        lw=0,
                        color="black",
                        alpha=0.3,
                        label=total_label,
                    )
                ax.plot(
                    pdos.energy,
                    pdos[projector],
                    "-",
                    lw=0.8,
                    color="black",
                    alpha=0.8,
                    label=projector,
                )
        if interactive:
            ax.legend(
                loc=(1.025, 0.2),
                bbox_transform=ax.transAxes,
                draggable=True,
                fontsize=legend_fontsize,
            )
        else:
            ax.legend(
                loc=(1.025, 0.2), bbox_transform=ax.transAxes, fontsize=legend_fontsize
            )

    if interactive:
        plt.show()
    else:
        plt.savefig(f"{output_name}.png", dpi=600, bbox_inches="tight")
        if save_pickle:
            import pickle

            with open(f"{output_name}.pickle", "wb") as file:
                pickle.dump(fig, file)
    plt.close()

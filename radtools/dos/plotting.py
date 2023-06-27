import matplotlib.pyplot as plt
import numpy as np

from radtools.dos.pdos import PDOS

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
        Fontsize of the title
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
        ax.vlines(
            0,
            0,
            1,
            transform=ax.get_xaxis_transform(),
            color="grey",
            linewidths=0.5,
            linestyles="dashed",
        )
        if title is not None and (i == 0 or relative):
            ax.set_title(title, fontsize=title_fontsize)

    if normalize:
        pdos = pdos.normalized()

    if relative:
        set_up_axis(ax, n - 1)
        ax.hlines(
            0,
            0,
            1,
            transform=ax.get_yaxis_transform(),
            color="black",
            linewidths=1,
        )

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
                if i == 0:
                    ax.plot(
                        pdos.energy,
                        np.where(pdos.ldos[0] > 1e-5, pdos.ldos[0], None),
                        "-",
                        lw=1,
                        color="blue",
                        alpha=0.8,
                        label=label_up,
                    )
                    ax.plot(
                        pdos.energy,
                        np.where(pdos.ldos[1] > 1e-5, -pdos.ldos[1], None),
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
                if i == 0:
                    ax.plot(
                        pdos.energy,
                        np.where(pdos.ldos > 1e-5, pdos.ldos, None),
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

            with open(f"{output_name}.png.pickle", "wb") as file:
                pickle.dump(fig, file)
    plt.close()


def plot_fatbands(
    pdos: PDOS,
    efermi=0.0,
    output_name="fatbands",
    ylim=None,
    interactive=False,
    save_pickle=False,
    colours=COLOURS,
    total_label="default",
    legend_fontsize=14,
    axes_labels_fontsize=12,
    title_fontsize=18,
):
    r"""
    Plot fatbands.

    Parameters
    ----------
    pdos : :py:class:`.PDOS`
        PDOS for the plot.
    efermi : float, default 0
        Fermi energy.
    output_name : str, default "pdos"
        output_name for the plot file. Extension ".png" is added at the end.
    ylim : tuple
        limits for the y (PDOS) axis
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
        Fontsize of the title
    """

    POINT_MAX_SIZE = 16
    norm = np.amax(pdos.pdos)
    if abs(norm) > 1e-8:
        pdos.pdos = pdos.pdos / norm

    kpoints = np.tile(pdos.kpoints, pdos.n_e).reshape((pdos.n_e, pdos.n_k)).T
    energy = np.tile(pdos.energy, pdos.n_k).reshape((pdos.n_k, pdos.n_e))

    if pdos.spin_pol:
        fig, axs = plt.subplots(2, 1, figsize=(10, 9))
        axs[1].get_yaxis().set_visible(False)
        axs[1].set_xlabel("k path", fontsize=axes_labels_fontsize)
        if ylim is not None:
            axs[1].set_ylim(*tuple(ylim))
        axs[0].title("Spin-up", fontsize=title_fontsize)
        axs[1].title("Spin-down", fontsize=title_fontsize)
        axs[1].set_xlim(kpoints[0][0], kpoints[-1][0])
    else:
        fig, axs = plt.subplots(figsize=(5, 9))
        axs = [axs]

    fig.subplots_adjust(wspace=0)

    if efermi == 0:
        axs[0].set_ylabel("E, eV", fontsize=axes_labels_fontsize)
    else:
        axs[0].set_ylabel("E - E$_{Fermi}$, eV", fontsize=axes_labels_fontsize)

    axs[0].set_xlabel("k path", fontsize=axes_labels_fontsize)
    if ylim is not None:
        axs[0].set_ylim(*tuple(ylim))
    axs[0].set_xlim(kpoints[0][0], kpoints[-1][0])

    for i, projector in enumerate(pdos):
        if pdos.spin_pol:
            axs[0].scatter(
                kpoints,
                energy,
                s=pdos[projector][0] * POINT_MAX_SIZE,
                color=colours[i % len(colours)],
                alpha=0.8,
                linewidth=0,
            )
            axs[1].scatter(
                kpoints,
                energy,
                s=pdos[projector][1] * POINT_MAX_SIZE,
                color=colours[i % len(colours)],
                label=f"{projector}",
                alpha=0.8,
                linewidth=0,
            )
        else:
            axs[0].scatter(
                kpoints,
                energy,
                s=pdos[projector] * POINT_MAX_SIZE,
                color=colours[i % len(colours)],
                label=f"{projector}",
                alpha=0.8,
                linewidth=0,
            )

    if interactive:
        if pdos.spin_pol:
            axs[1].legend(
                loc=(1.025, 0.2),
                bbox_transform=axs[1].transAxes,
                draggable=True,
                fontsize=legend_fontsize,
            )
        else:
            axs[0].legend(
                loc=(1.025, 0.2),
                bbox_transform=axs[0].transAxes,
                draggable=True,
                fontsize=legend_fontsize,
            )
    else:
        if pdos.spin_pol:
            axs[1].legend(
                loc=(1.025, 0.2),
                bbox_transform=axs[1].transAxes,
                fontsize=legend_fontsize,
            )
        else:
            axs[0].legend(
                loc=(1.025, 0.2),
                bbox_transform=axs[0].transAxes,
                fontsize=legend_fontsize,
            )

    if interactive:
        plt.show()
    else:
        plt.savefig(f"{output_name}.png", dpi=600, bbox_inches="tight")
        if save_pickle:
            import pickle

            with open(f"{output_name}.png.pickle", "wb") as file:
                pickle.dump(fig, file)
    plt.close()

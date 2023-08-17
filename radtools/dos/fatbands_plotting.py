from os.path import abspath, isfile, join

import matplotlib.pyplot as plt
import numpy as np
from termcolor import cprint

from radtools.decorate.axes import plot_hlines, plot_vlines
from radtools.decorate.colormap import custom_cmap
from radtools.dos.dos import prepare_custom_pdos
from radtools.dos.pdos import PDOS
from radtools.dos.plotting import COLOURS

__all__ = [
    "plot_custom_fatbands",
    "plot_fatbands",
]


def plot_custom_fatbands(
    dos,
    custom,
    output_root=".",
    energy_window=None,
    k_window=None,
    efermi=0.0,
    interactive=False,
    separate=False,
    save_pickle=False,
    save_txt=False,
    colours=COLOURS,
    labels=None,
    legend_fontsize=12,
    axes_labels_fontsize=14,
    title_fontsize=18,
    k_points=None,
):
    cprint("Plotting custom plot", "green")
    print("Input is understood as:")
    projectors = []

    # Check if labels are provided and correct
    if labels is not None and len(labels) != len(custom):
        raise ValueError(
            f"Got {len(labels)} labels, but {len(custom)} PDOS, have to be the same."
        )

    # Set projectors
    if labels is None:
        projectors = custom
    else:
        projectors = labels

    # Get PDOS array
    pdos = prepare_custom_pdos(dos=dos, custom=custom, quiet=False)

    # Create PDOS object
    pdos = PDOS(
        energy=dos.energy,
        pdos=pdos,
        projectors_group="Total (sum)",
        projectors=projectors,
        spin_pol=dos.case in [2, 3],
    )

    # Compute output name
    if isfile(join(output_root, "custom.png")):
        i = 1
        while isfile(join(output_root, f"custom{i}.png")):
            i += 1
        output_name = f"custom{i}"
    else:
        output_name = "custom"

    # Save txt
    if save_txt:
        pdos.dump_txt(join(output_root, f"{output_name}.txt"))

    # Plot
    plot_fatbands(
        pdos=pdos,
        efermi=efermi,
        output_name=join(output_root, output_name),
        ylim=energy_window,
        xlim=k_window,
        interactive=interactive,
        separate=separate,
        save_pickle=save_pickle,
        colours=colours,
        legend_fontsize=legend_fontsize,
        axes_labels_fontsize=axes_labels_fontsize,
        title_fontsize=title_fontsize,
        k_points=k_points,
    )
    cprint(f"Result is in {abspath(join(output_root, f'{output_name}.png'))}", "blue")


def plot_fatbands(
    pdos: PDOS,
    efermi=0.0,
    output_name="fatbands",
    ylim=None,
    xlim=None,
    interactive=False,
    separate=False,
    save_pickle=False,
    colours=COLOURS,
    legend_fontsize=14,
    axes_labels_fontsize=12,
    title_fontsize=18,
    k_points=None,
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
    separate : bool, default False
        Whether to plot each entry in an individual figure.
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

    norm = np.amax(pdos.pdos)

    if xlim is None:
        xlim = (pdos.kpoints[0], pdos.kpoints[-1])
    if ylim is None:
        ylim = (pdos.energy[0], pdos.energy[-1])

    if k_points is not None:
        tmp = k_points
        k_points = [[], []]
        for i in range(len(tmp) // 2):
            k_points[1].append(tmp[2 * i])
            k_points[0].append(float(tmp[2 * i + 1]))

    if not separate:
        pass

    def plot_entry(axs, pdos, projector, xlim, ylim, efermi):
        if pdos.spin_pol:
            axs[1].get_yaxis().set_visible(False)
            if k_points is not None:
                axs[1].set_xticks(
                    k_points[0], k_points[1], fontsize=axes_labels_fontsize
                )
                plot_vlines(axs[1], k_points[0])
            else:
                axs[1].set_xlabel("k path", fontsize=axes_labels_fontsize)
            axs[1].set_ylim(*tuple(ylim))
            axs[1].set_xlim(*tuple(xlim))
        if k_points is not None:
            axs[0].set_xticks(k_points[0], k_points[1], fontsize=axes_labels_fontsize)
            plot_vlines(axs[0], k_points[0])
        else:
            axs[0].set_xlabel("k path", fontsize=axes_labels_fontsize)
        axs[0].set_ylim(*tuple(ylim))
        axs[0].set_xlim(*tuple(xlim))

        if efermi == 0:
            axs[0].set_ylabel("E, eV", fontsize=axes_labels_fontsize)
        else:
            axs[0].set_ylabel("E - E$_{Fermi}$, eV", fontsize=axes_labels_fontsize)

        if pdos.spin_pol:
            axs[0].imshow(
                pdos[projector][0].T,
                cmap=custom_cmap((1, 1, 1), (0, 0, 1)),
                origin="lower",
                extent=(
                    xlim[0],
                    xlim[1],
                    ylim[0],
                    ylim[1],
                ),
                aspect="auto",
                vmax=norm,
                vmin=0,
            )
            axs[0].set_title(f"{projector} (up)", fontsize=title_fontsize)
            axs[1].imshow(
                pdos[projector][1].T,
                cmap=custom_cmap((1, 1, 1), (1, 0, 0)),
                origin="lower",
                extent=(
                    xlim[0],
                    xlim[1],
                    ylim[0],
                    ylim[1],
                ),
                aspect="auto",
                vmax=norm,
                vmin=0,
            )
            axs[1].set_title(f"{projector} (down)", fontsize=title_fontsize)
        else:
            axs[0].imshow(
                pdos[projector].T,
                cmap="inferno",
                origin="lower",
                extent=(
                    xlim[0],
                    xlim[1],
                    ylim[0],
                    ylim[1],
                ),
                aspect="auto",
                vmax=norm,
                vmin=0,
            )
            axs[0].set_title(projector, fontsize=title_fontsize)

        if efermi != 0:
            plot_hlines(axs[0], 0)
            if pdos.spin_pol:
                plot_hlines(axs[1], 0)

    for i, projector in enumerate(pdos):
        if separate:
            if pdos.spin_pol:
                fig, axs = plt.subplots(1, 2, figsize=(8, 7))
            else:
                fig, axs = plt.subplots(figsize=(4, 7))
                axs = [axs]
            fig.subplots_adjust(wspace=0)
            plot_entry(axs, pdos, projector, xlim, ylim, efermi)
            if interactive:
                plt.show()
            else:
                plt.savefig(f"{output_name}_{i}.png", dpi=600, bbox_inches="tight")
                if save_pickle:
                    import pickle

                    with open(f"{output_name}_{i}.png.pickle", "wb") as file:
                        pickle.dump(fig, file)
            plt.close()
        else:
            plot_entry(axs[i], pdos, projector, xlim, ylim, efermi)

    if not separate:
        if interactive:
            plt.show()
        else:
            plt.savefig(f"{output_name}.png", dpi=600, bbox_inches="tight")
            if save_pickle:
                import pickle

                with open(f"{output_name}.png.pickle", "wb") as file:
                    pickle.dump(fig, file)
        plt.close()

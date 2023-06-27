import matplotlib.pyplot as plt
import numpy as np
from termcolor import cprint
from os.path import isfile, join, abspath
from os import makedirs
from tqdm import tqdm

from radtools.dos.pdos import PDOS

__all__ = [
    "COLOURS",
    "plot_projected",
    "plot_custom_pdos",
    "prepare_custom_pdos",
    "plot_custom_fatbands",
    "plot_fatbands",
]

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


def prepare_custom_pdos(
    dos,
    custom,
):
    r"""
    Prepare custom PDOS. Based on the input line.

    Parameters
    ----------
    dos : :py:class:`.DOSQE`
        DOS input files wrapper.
    custom : list of str
        List of strings describing the custom PDOS.

    Returns
    -------
    pdos : list of :numpy:`ndarray`
        PDOS array.
    """
    pdos = []

    for entry in custom:
        cprint(f'"{entry}":', "green")
        entry = entry.replace(" ", "").replace(")", "")
        subentries = entry.split(";")
        pdos_element = None
        for subentry in subentries:
            atom_part = subentry.split("(")[0]
            atom = atom_part.split("#")[0]
            if "#" in subentry:
                atom_numbers = list(map(int, atom_part.split("#")[1:]))
            else:
                atom_numbers = dos.atom_numbers(atom)

            cprint(
                f"  * PDOS is summed among the following {atom} atoms:",
                "green",
                end="\n      ",
            )
            for i, number in enumerate(atom_numbers):
                print(f"{atom}#{number}", end="")
                if i != len(atom_numbers) - 1:
                    print(end=", ")
                else:
                    print()

            if "(" in subentry:
                wfc_parts = subentry.split("(")[1].split(",")
                wfcs = []
                for wfc_part in wfc_parts:
                    wfc = wfc_part.split("#")[0]
                    if "#" in wfc_part:
                        wfc_numbers = list(map(int, wfc_part.split("#")[1:]))
                        for number in wfc_numbers:
                            wfcs.append((wfc, number))
                    else:
                        wfc_list = dos.wfcs(atom)
                        for name, number in wfc_list:
                            if name == wfc:
                                wfcs.append((wfc, number))
            else:
                wfcs = dos.wfcs(atom)

            print(
                f"  * For each {atom} atom PDOS is summed among the following projections:",
                end="\n      ",
            )
            for i, (name, number) in enumerate(wfcs):
                print(f"{name}#{number}", end="")
                if i != len(wfcs) - 1:
                    print(end=", ")
                else:
                    print()

            for name, number in wfcs:
                if pdos_element is None:
                    pdos_element = dos.pdos(
                        atom=atom,
                        wfc=name,
                        wfc_number=number,
                        atom_numbers=atom_numbers,
                    ).ldos
                else:
                    pdos_element += dos.pdos(
                        atom=atom,
                        wfc=name,
                        wfc_number=number,
                        atom_numbers=atom_numbers,
                    ).ldos

        pdos.append(pdos_element)

    return pdos


def plot_custom_pdos(
    dos,
    custom,
    output_root=".",
    energy_window=None,
    dos_window=None,
    efermi=0.0,
    relative=False,
    normalize=False,
    interactive=False,
    save_pickle=False,
    save_txt=False,
    background_total=False,
    colours=COLOURS,
    labels=None,
    legend_fontsize=12,
    axes_labels_fontsize=14,
    title_fontsize=18,
):
    cprint("Plotting custom plot", "green")
    print("Input is understood as:")
    projectors = []

    # Check if labels are provided and correct
    if (
        labels is not None
        and len(labels) != len(custom)
        and len(labels) != len(custom) + 1
    ):
        raise ValueError(
            f"Got {len(labels)} labels, but {len(custom)} PDOS, have to be the same or n custom, n+1 labels."
        )

    # Process some of the predefined labels
    if labels is not None and len(labels) == len(custom) + 1:
        if labels[0].lower() == "none":
            total_label = None
        elif labels[0].lower() == "default":
            total_label = "default"
        else:
            total_label = labels[0]
        labels = labels[1:]
    else:
        total_label = "default"

    # Set projectors
    if labels is None:
        projectors = custom
    else:
        projectors = labels

    # Get PDOS array
    pdos = prepare_custom_pdos(dos=dos, custom=custom)

    # Create PDOS object
    if background_total:
        pdos = PDOS(
            energy=dos.energy,
            pdos=pdos,
            ldos=dos.total_pdos(),
            projectors_group="Total PDOS",
            projectors=projectors,
            spin_pol=dos.spin_pol,
        )
    else:
        pdos = PDOS(
            energy=dos.energy,
            pdos=pdos,
            projectors_group="Total (sum)",
            projectors=projectors,
            spin_pol=dos.spin_pol,
        )

    # Compute output_name
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
    plot_projected(
        pdos=pdos,
        efermi=efermi,
        output_name=join(output_root, output_name),
        xlim=energy_window,
        ylim=dos_window,
        relative=relative,
        normalize=normalize,
        interactive=interactive,
        save_pickle=save_pickle,
        colours=colours,
        total_label=total_label,
        legend_fontsize=legend_fontsize,
        axes_labels_fontsize=axes_labels_fontsize,
        title_fontsize=title_fontsize,
    )
    cprint(f"Result is in {abspath(join(output_root, f'{output_name}.png'))}", "blue")


def plot_orbital_resolved(
    dos,
    output_root=".",
    energy_window=None,
    dos_window=None,
    efermi=0.0,
    separate=False,
    relative=False,
    normalize=False,
    interactive=False,
    save_pickle=False,
    save_txt=False,
    background_total=False,
    colours=COLOURS,
    legend_fontsize=12,
    axes_labels_fontsize=14,
    title_fontsize=18,
):
    cprint("Orbital-resolved PDOS:", "green")
    local_output = join(output_root, "orbital-resolved")
    makedirs(local_output, exist_ok=True)

    data = {}
    for atom, atom_number, wfc, wfc_number in dos:
        if atom not in data:
            data[atom] = []
        data[atom].append((wfc, wfc_number))

    # Avoid repetitions
    for atom in data:
        data[atom] = list(set(data[atom]))

    for atom in data:
        for wfc, wfc_number in data[atom]:
            if separate:
                atom_numbers = dos.atom_numbers(atom)
            else:
                atom_numbers = [None]

            for atom_number in tqdm(atom_numbers, desc=f"  {atom} {wfc} #{wfc_number}"):
                if separate:
                    atom_name = f"{atom}#{atom_number}"
                else:
                    atom_name = atom

                title = f"PDOS for {atom_name} ({wfc} #{wfc_number})"

                pdos = dos.pdos(
                    atom=atom,
                    wfc=wfc,
                    wfc_number=wfc_number,
                    atom_numbers=atom_number,
                    background_total=background_total,
                )
                if background_total:
                    pdos.projectors_group = "Total PDOS"
                if save_txt:
                    pdos.dump_txt(
                        join(local_output, f"{atom_name}_{wfc}#{wfc_number}.txt")
                    )
                plot_projected(
                    pdos=pdos,
                    efermi=efermi,
                    output_name=join(local_output, f"{atom_name}_{wfc}#{wfc_number}"),
                    title=title,
                    xlim=energy_window,
                    ylim=dos_window,
                    relative=relative,
                    normalize=normalize,
                    interactive=interactive,
                    save_pickle=save_pickle,
                    colours=colours,
                    legend_fontsize=legend_fontsize,
                    axes_labels_fontsize=axes_labels_fontsize,
                    title_fontsize=title_fontsize,
                )
    cprint(f"Results are in {abspath(local_output)}", "blue")


def plot_atom_resolved(
    dos,
    output_root=".",
    energy_window=None,
    dos_window=None,
    efermi=0.0,
    separate=False,
    relative=False,
    normalize=False,
    interactive=False,
    save_pickle=False,
    save_txt=False,
    background_total=False,
    colours=COLOURS,
    legend_fontsize=12,
    axes_labels_fontsize=14,
    title_fontsize=18,
):
    cprint("Orbital's contribution for each atom.", "green")
    local_output = join(output_root, "atom-resolved")
    makedirs(local_output, exist_ok=True)

    for atom in dos.atoms:
        if separate:
            atom_numbers = dos.atom_numbers(atom)
        else:
            atom_numbers = [None]
        for atom_number in tqdm(atom_numbers, desc=f"  {atom}"):
            if separate:
                atom_name = f"{atom}#{atom_number}"
            else:
                atom_name = atom
            projectors = []
            pdos = []
            for wfc, wfc_number in dos.wfcs(atom, atom_number):
                projectors.append(f"{wfc} #{wfc_number}")
                pdos.append(dos.pdos(atom, wfc, wfc_number, atom_number).ldos)

            if background_total:
                pdos = PDOS(
                    energy=dos.energy,
                    pdos=pdos,
                    ldos=dos.total_pdos(),
                    projectors_group=atom_name,
                    projectors=projectors,
                    spin_pol=dos.spin_pol,
                )
                pdos.projectors_group = "Total PDOS"
            else:
                pdos = PDOS(
                    energy=dos.energy,
                    pdos=pdos,
                    projectors_group=atom_name,
                    projectors=projectors,
                    spin_pol=dos.spin_pol,
                )
            title = f"PDOS for {atom_name}"
            if save_txt:
                pdos.dump_txt(join(local_output, f"{atom_name}.txt"))
            plot_projected(
                pdos=pdos,
                efermi=efermi,
                output_name=join(local_output, atom_name),
                title=title,
                xlim=energy_window,
                ylim=dos_window,
                relative=relative,
                normalize=normalize,
                interactive=interactive,
                save_pickle=save_pickle,
                colours=colours,
                legend_fontsize=legend_fontsize,
                axes_labels_fontsize=axes_labels_fontsize,
                title_fontsize=title_fontsize,
            )
    cprint(f"Results are in {abspath(local_output)}", "blue")


def plot_atom_to_total(
    dos,
    output_root=".",
    energy_window=None,
    dos_window=None,
    efermi=0.0,
    separate=False,
    relative=False,
    normalize=False,
    interactive=False,
    save_pickle=False,
    save_txt=False,
    background_total=False,
    colours=COLOURS,
    legend_fontsize=12,
    axes_labels_fontsize=14,
    title_fontsize=18,
):
    cprint("Atom's contributions into total PDOS:", "green")
    projectors = []
    pdos = []
    for atom in dos.atoms:
        if separate:
            atom_numbers = dos.atom_numbers(atom)
        else:
            atom_numbers = [None]
        for atom_number in atom_numbers:
            if separate:
                atom_name = f"{atom}#{atom_number}"
            else:
                atom_name = atom
            projectors.append(atom_name)
            for i, (wfc, wfc_number) in enumerate(dos.wfcs(atom, atom_number)):
                if i == 0:
                    ldos = dos.pdos(atom, wfc, wfc_number, atom_number).ldos
                else:
                    ldos += dos.pdos(atom, wfc, wfc_number, atom_number).ldos
            pdos.append(ldos)
    if background_total:
        pdos = PDOS(
            energy=dos.energy,
            pdos=pdos,
            ldos=dos.total_pdos(),
            projectors_group="Total PDOS",
            projectors=projectors,
            spin_pol=dos.spin_pol,
        )
        pdos.projectors_group = "Total PDOS"
    else:
        pdos = PDOS(
            energy=dos.energy,
            pdos=pdos,
            projectors_group="Total PDOS",
            projectors=projectors,
            spin_pol=dos.spin_pol,
        )

    title = f"Atom contribution in PDOS"
    if save_txt:
        pdos.dump_txt(join(output_root, "atomic-contributions.txt"))
    plot_projected(
        pdos=pdos,
        efermi=efermi,
        output_name=join(output_root, "atomic-contributions"),
        title=title,
        xlim=energy_window,
        ylim=dos_window,
        relative=relative,
        normalize=normalize,
        interactive=interactive,
        save_pickle=save_pickle,
        colours=colours,
        legend_fontsize=legend_fontsize,
        axes_labels_fontsize=axes_labels_fontsize,
        title_fontsize=title_fontsize,
    )
    cprint(f"Result is in {abspath(output_root)}", "blue")


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


def plot_custom_fatbands(
    dos,
    custom,
    output_root=".",
    energy_window=None,
    k_window=None,
    efermi=0.0,
    interactive=False,
    save_pickle=False,
    save_txt=False,
    colours=COLOURS,
    labels=None,
    legend_fontsize=12,
    axes_labels_fontsize=14,
    title_fontsize=18,
    scale_points=1,
):
    cprint("Plotting custom plot", "green")
    print("Input is understood as:")
    projectors = []

    # Check if labels are provided and correct
    if (
        labels is not None
        and len(labels) != len(custom)
        and len(labels) != len(custom) + 1
    ):
        raise ValueError(
            f"Got {len(labels)} labels, but {len(custom)} PDOS, have to be the same or n custom, n+1 labels."
        )

    # Process some of the predefined labels
    if labels is not None and len(labels) == len(custom) + 1:
        if labels[0].lower() == "none":
            total_label = None
        elif labels[0].lower() == "default":
            total_label = "default"
        else:
            total_label = labels[0]
        labels = labels[1:]
    else:
        total_label = "default"

    # Set projectors
    if labels is None:
        projectors = custom
    else:
        projectors = labels

    # Get PDOS array
    pdos = prepare_custom_pdos(dos=dos, custom=custom)

    # Create PDOS object
    pdos = PDOS(
        energy=dos.energy,
        pdos=pdos,
        ldos=dos.total_pdos(),
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
        save_pickle=save_pickle,
        colours=colours,
        total_label=total_label,
        legend_fontsize=legend_fontsize,
        axes_labels_fontsize=axes_labels_fontsize,
        title_fontsize=title_fontsize,
        scale_points=scale_points,
    )
    cprint(f"Result is in {abspath(join(output_root, f'{output_name}.png'))}", "blue")


def plot_fatbands(
    pdos: PDOS,
    efermi=0.0,
    output_name="fatbands",
    ylim=None,
    xlim=None,
    interactive=False,
    save_pickle=False,
    colours=COLOURS,
    total_label="default",
    legend_fontsize=14,
    axes_labels_fontsize=12,
    title_fontsize=18,
    scale_points=1,
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
    point_max_size : int, default 16
        Maximum size of the points in the plot.
    """

    point_max_size = 16 * scale_points

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
        if xlim is not None:
            axs[1].set_xlim(*tuple(xlim))
        else:
            axs[1].set_xlim(kpoints[0][0], kpoints[-1][0])
        axs[0].title("Spin-up", fontsize=title_fontsize)
        axs[1].title("Spin-down", fontsize=title_fontsize)
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
    if xlim is not None:
        axs[0].set_xlim(*tuple(xlim))
    else:
        axs[0].set_xlim(kpoints[0][0], kpoints[-1][0])
    axs[0].set_xlim(kpoints[0][0], kpoints[-1][0])

    for i, projector in enumerate(pdos):
        if pdos.spin_pol:
            axs[0].scatter(
                kpoints,
                energy,
                s=pdos[projector][0] * point_max_size,
                color=colours[i % len(colours)],
                alpha=0.8,
                linewidth=0,
            )
            axs[1].scatter(
                kpoints,
                energy,
                s=pdos[projector][1] * point_max_size,
                color=colours[i % len(colours)],
                label=f"{projector}",
                alpha=0.8,
                linewidth=0,
            )
        else:
            axs[0].scatter(
                kpoints,
                energy,
                s=pdos[projector] * point_max_size,
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

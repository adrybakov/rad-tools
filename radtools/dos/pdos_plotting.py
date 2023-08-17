from os import makedirs
from os.path import abspath, isfile, join

from termcolor import cprint
from tqdm import tqdm

from radtools.dos.dos import prepare_custom_pdos
from radtools.dos.pdos import PDOS
from radtools.dos.plotting import COLOURS, plot_projected

__all__ = [
    "plot_custom_pdos",
    "plot_orbital_resolved",
    "plot_atom_resolved",
    "plot_atom_to_total",
]


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
    pdos = prepare_custom_pdos(dos=dos, custom=custom, quiet=False)

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
                    wfc_numbers=wfc_number,
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
            for wfc in dos.wfcs(atom, atom_number):
                for wfc_number in dos.wfc_numbers(atom, wfc, atom_number):
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
            for i, wfc in enumerate(dos.wfcs(atom, atom_number)):
                if i == 0:
                    ldos = dos.pdos(atom, wfc, None, atom_number).ldos
                else:
                    ldos += dos.pdos(atom, wfc, None, atom_number).ldos
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

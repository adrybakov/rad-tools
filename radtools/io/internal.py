r"""
Input-output for the files related with this package.
"""

__all__ = [
    "load_template",
    "read_template",
    "dump_spinham_txt",
    "dump_pickle",
    "load_pickle",
]

import numpy as np
from radtools.spinham.template import ExchangeTemplate
from radtools.decorate.array import print_2d_array
from radtools.decorate.stats import logo
from radtools.spinham.constants import TXT_FLAGS
from radtools.spinham.hamiltonian import SpinHamiltonian

meV_TO_J = 1.602176634e-22


def dump_pickle(object, filename):
    """
    Save any python Object in a binary format.

    Parameters
    ----------
    filename : str
        Name of the file for the Object to be saved in.
        ".pickle" is added automatically.
    """
    import pickle

    with open(f"{filename}.pickle", "wb") as file:
        pickle.dump(object, file)


def load_pickle(filename):
    r"""
    Load any python Object from a binary format.

    Parameters
    ----------
    filename : str
        Name of the .pickle file for the Object to be loaded from.
    """

    import pickle

    with open(filename, "rb") as file:
        object = pickle.load(file)
    return object


def load_template(filename):
    r"""
    Read template from the template file.

    .. versionchanged:: 0.8.1 Renamed from ``read_template``

    See :ref:`template-draft` for the description of the template file format.

    Parameters
    ----------
    filename : str
        Path to the template file.

    Returns
    -------
    template : :py:class:`.ExchangeTemplate`
        Exchange template.

    Notes
    -----
    See also :ref:`Template specification <api>`.
    """

    # Constants
    major_sep = "=" * 20
    minor_sep = "-" * 20
    neighbors_flag = "Neighbors template:"

    template = ExchangeTemplate()

    with open(filename) as file:
        line = True
        while line:
            line = file.readline()

            # Read the neighbors names
            if line and neighbors_flag in line:
                while line and major_sep not in line:
                    if line and minor_sep in line:
                        line = file.readline()
                        name = line.split()[0]
                        try:
                            latex_name = line.split()[1]
                        except IndexError:
                            latex_name = name
                        template.names[name] = []
                        template.latex_names[name] = latex_name
                        line = file.readline()
                        while line and minor_sep not in line and major_sep not in line:
                            atom1 = line.split()[0]
                            atom2 = line.split()[1]

                            R = tuple(map(int, line.split()[2:]))
                            template.names[name].append((atom1, atom2, R))
                            line = file.readline()
                    if line and major_sep in line:
                        break
                    if line and minor_sep not in line:
                        line = file.readline()
    return template


# For backward compatibility
read_template = load_template


def dump_spinham_txt(
    spinham: SpinHamiltonian,
    filename=None,
    anisotropic=True,
    matrix=True,
    dmi=True,
    template=None,
    decimals=4,
    additional_stats=None,
):
    """
    Save the :py:class:`.SpinHamiltonian` in a human-readable format.

    Parameters
    ----------
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian to be saved.
    filename : str, optional
        Name of the file for the Hamiltonian to be saved in.
        If not given, the Hamiltonian will be printed in the console.
    anisotropic : bool, default True
        Whether to output anisotropic exchange.
    matrix : bool, default True
        Whether to output whole matrix exchange.
    dmi : bool, default True
        Whether to output DMI exchange.
    template : :py:class:`.ExchangeTemplate`, optional
        If provided, then not the SpinHamiltonian will be written, but the model
        based on the template.
    decimals : int, default 4
        Number of decimals to be printed (only for the exchange values).
    additional_stats : str, optional
        Additional info, which will be printed right after the logo, before the
        main separator.
    """

    main_separator = "=" * 80 + "\n"
    separator = "-" * 80 + "\n"
    spinham_txt = []
    fmt = f"{decimals+4}.{decimals}f"

    spinham_txt.append(main_separator)
    spinham_txt.append(logo(date_time=True, line_length=80) + "\n")
    if additional_stats is not None:
        spinham_txt.append(additional_stats)
    spinham_txt.append(main_separator)
    spinham_txt.append(TXT_FLAGS["cell"] + "\n")
    spinham_txt.append(
        print_2d_array(
            spinham.cell,
            borders=False,
            fmt="^.8f",
            print_result=False,
            header_row=["x", "y", "z"],
        )
        + "\n"
    )
    spinham_txt.append(main_separator)
    spinham_txt.append(TXT_FLAGS["atoms"] + "\n")
    header_column = []
    header_row = [
        f"Index Name",
        "a1 (rel)",
        "a2 (rel)",
        "a3 (rel)",
        "x",
        "y",
        "z",
    ]
    atom_data = np.zeros((len(spinham.atoms), 6))
    for a_i, atom in enumerate(spinham.atoms):
        header_column.append(f"{atom.index:<5} {atom.name:<4}")
        atom_data[a_i, :3] = spinham.get_atom_coordinates(atom, relative=True)
        atom_data[a_i, 3:] = spinham.get_atom_coordinates(atom, relative=False)
    spinham_txt.append(
        print_2d_array(
            atom_data,
            borders=False,
            fmt="^.8f",
            print_result=False,
            header_row=header_row,
            header_column=header_column,
        )
        + "\n"
    )
    spinham_txt.append(main_separator)

    write_spins = True
    write_magmoms = True
    for atom in spinham.magnetic_atoms:
        try:
            atom.magmom
        except ValueError:
            write_magmoms = False
        try:
            atom.spin
            atom.spin_vector
        except ValueError:
            write_spins = False
    if write_magmoms:
        spinham_txt.append(TXT_FLAGS["magmoms"] + "\n")
        header_column = []
        header_row = [f"Index Name", "m_x", "m_y", "m_z"]
        magmom_data = np.zeros((len(spinham.magnetic_atoms), 3))
        for a_i, atom in enumerate(spinham.magnetic_atoms):
            header_column.append(f"{atom.index:<5} {atom.name:<4}")
            magmom_data[a_i] = atom.magmom
        spinham_txt.append(
            print_2d_array(
                magmom_data,
                borders=False,
                fmt="^.8f",
                print_result=False,
                header_row=header_row,
                header_column=header_column,
            )
            + "\n"
        )
        spinham_txt.append(main_separator)
    if write_spins:
        spinham_txt.append(TXT_FLAGS["spins"] + "\n")
        header_column = []
        header_row = [f"Index Name", "S", "S_x", "S_y", "S_z"]
        spin_data = np.zeros((len(spinham.magnetic_atoms), 4))
        for a_i, atom in enumerate(spinham.magnetic_atoms):
            header_column.append(f"{atom.index:<5} {atom.name:<4}")
            spin_data[a_i, 0] = atom.spin
            spin_data[a_i, 1:] = atom.spin_vector
        spinham_txt.append(
            print_2d_array(
                spin_data,
                borders=False,
                fmt="^.8f",
                print_result=False,
                header_row=header_row,
                header_column=header_column,
            )
            + "\n"
        )
        spinham_txt.append(main_separator)

    spinham_txt.append(TXT_FLAGS["notation"] + "\n")
    spinham_txt.append(f"{spinham.notation}\n")
    spinham_txt.append(spinham.notation_string + "\n")
    spinham_txt.append(main_separator)

    if template is None:
        spinham_txt.append(TXT_FLAGS["spinham"] + "\n")
        spinham_txt.append(
            f"{'Atom1':6} {'Atom2':6} (  i,   j,   k) {'J_iso':^{decimals+4}} {'Distance':^8}\n"
        )
        bonds_data = []
        for atom1, atom2, (i, j, k), J in spinham:
            data_entry = []
            distance = spinham.get_distance(atom1, atom2, (i, j, k))
            atom1 = f"{atom1.name}({atom1.index})"
            atom2 = f"{atom2.name}({atom2.index})"
            data_entry.append(separator)
            data_entry.append(
                f"{atom1:6} {atom2:6} "
                + f"({i:>3}, {j:>3}, {k:>3}) "
                + f"{J.iso:^{decimals+4}.{decimals}f} "
                + f"{distance:^8.4f}\n"
            )
            if matrix:
                data_entry.append(TXT_FLAGS["matrix"] + "\n")
                data_entry.append(
                    print_2d_array(
                        J.matrix, fmt=fmt, borders=False, print_result=False, shift=2
                    )
                    + "\n"
                )
            if anisotropic:
                data_entry.append(TXT_FLAGS["aniso"] + "\n")
                data_entry.append(
                    print_2d_array(
                        J.aniso, fmt=fmt, borders=False, print_result=False, shift=2
                    )
                    + "\n"
                )
            if dmi:
                data_entry.append(TXT_FLAGS["dmi_module"] + "\n")
                data_entry.append(f"  {J.dmi_module:{decimals+4}.{decimals}f}\n")
                data_entry.append(TXT_FLAGS["dmi_relative"] + "\n")
                data_entry.append(f"  {J.rel_dmi:{decimals+4}.{decimals}f}\n")
                data_entry.append(TXT_FLAGS["dmi"] + "\n")
                data_entry.append(
                    print_2d_array(
                        J.dmi, fmt=fmt, borders=False, print_result=False, shift=2
                    )
                    + "\n"
                )
            bonds_data.append(["".join(data_entry), distance])
        bonds_data = sorted(bonds_data, key=lambda x: x[1])
        bonds_data = [x[0] for x in bonds_data]
        spinham_txt.append("".join(bonds_data))
    else:
        spinham = spinham.formed_model(template)
        spinham_txt.append(TXT_FLAGS["spinmodel"] + "\n")
        for parameter in template.names:
            bonds = template.names[parameter]
            J = spinham[bonds[0]]
            spinham_txt.append(separator)
            spinham_txt.append(f"{parameter} contains {len(bonds)} bonds\n")
            spinham_txt.append(TXT_FLAGS["iso"])
            spinham_txt.append(f"  {J.iso:{decimals+4}.{decimals}f}\n")
            if matrix:
                spinham_txt.append(TXT_FLAGS["matrix"] + "\n")
                spinham_txt.append(
                    print_2d_array(
                        J.matrix, fmt=fmt, borders=False, print_result=False, shift=2
                    )
                    + "\n"
                )
            if anisotropic:
                spinham_txt.append(TXT_FLAGS["aniso"] + "\n")
                spinham_txt.append(
                    print_2d_array(
                        J.aniso, fmt=fmt, borders=False, print_result=False, shift=2
                    )
                    + "\n"
                )
            if dmi:
                spinham_txt.append(TXT_FLAGS["dmi_module"] + "\n")
                spinham_txt.append(f"  {J.dmi_module:{decimals+4}.{decimals}f}\n")
                spinham_txt.append(TXT_FLAGS["dmi_relative"] + "\n")
                spinham_txt.append(f"  {J.rel_dmi:{decimals+4}.{decimals}f}\n")
                spinham_txt.append(TXT_FLAGS["dmis"] + "\n")
                for atom1, atom2, (i, j, k) in bonds:
                    J = spinham[atom1, atom2, (i, j, k)]
                    atom1 = spinham.get_atom(atom1)
                    atom2 = spinham.get_atom(atom2)
                    atom1 = f"{atom1.name}({atom1.index})"
                    atom2 = f"{atom2.name}({atom2.index})"
                    spinham_txt.append(
                        print_2d_array(
                            J.dmi, fmt=fmt, borders=False, print_result=False, shift=2
                        )
                    )
                    spinham_txt.append(
                        f"   ({atom1:6} {atom2:6} {i:>3}, {j:>3}, {k:>3})\n"
                    )
            else:
                spinham_txt.append(TXT_FLAGS["bonds"] + "\n")
                spinham_txt.append(f"   {'Atom1':6} {'Atom2':6}   i,   j,   k\n")
                for atom1, atom2, (i, j, k) in bonds:
                    atom1 = spinham.get_atom(atom1)
                    atom2 = spinham.get_atom(atom2)
                    atom1 = f"{atom1.name}({atom1.index})"
                    atom2 = f"{atom2.name}({atom2.index})"
                    spinham_txt.append(
                        f"   {atom1:6} {atom2:6} {i:>3}, {j:>3}, {k:>3}\n"
                    )
    spinham_txt.append(main_separator)
    spinham_txt = "".join(spinham_txt)
    if filename is not None:
        with open(filename, "w") as file:
            file.write(spinham_txt)
    else:
        print(spinham_txt)

r"""
Input-output for |Vampire|_.
"""

__all__ = ["dump_vampire", "dump_mat", "dump_ucf"]

from radtools.decorate.array import print_2d_array
from radtools.spinham.parameter import ExchangeParameter
from radtools.decorate.stats import logo

meV_TO_J = 1.602176634e-22


def dump_vampire(
    spinham,
    seedname="vampire",
    anisotropic=True,
    dmi=True,
    custom_mask=None,
    decimals=5,
    materials=None,
    nologo=False,
):
    """
    Save the Hamiltonian in a Human-readable format.

    Parameters
    ----------
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian to be saved.
    seedname : str, default "vampire"
        Seedname for the .UCF and .mat files. Extensions are added automatically.
        Input file always have the name "input".
    anisotropic : bool, default True
        Whether to output anisotropic exchange.
    dmi : bool, default True
        Whether to output DMI exchange.
    custom_mask : func, optional
        Custom mask for the exchange parameter. Function which take (3,3) numpy:`ndarray`
        as an input and returns (3,3) numpy:`ndarray` as an output. If given, then
        ``anisotropic`` and ``dmi`` parameters are ignored.
    decimals : int, default 4
        Number of decimals to be printed (only for the exchange values).
    materials : list of int, optional
        List of materials for the atoms. Has to have the same length as the number of
        magnetic atoms in the ``spinham``. Order is the same as in :py:attr:`.SpinHamiltonian.magnetic_atoms`.
        If not given, each atom will be considered as a separate material. Starting from 0.
    nologo : bool, default False
        Whether to print the logo in the output files.

    Returns
    -------
    content : str
        Content of the .UCF file if ``filename`` is not given.
    """

    dump_ucf(
        spinham,
        filename=seedname + ".UCF",
        anisotropic=anisotropic,
        dmi=dmi,
        custom_mask=custom_mask,
        decimals=decimals,
        materials=materials,
        nologo=nologo,
    )
    dump_mat(
        spinham,
        filename=seedname + ".mat",
        materials=materials,
        nologo=nologo,
    )
    with open("input", "w") as file:
        file.write(
            "#------------------------------------------\n"
            f"material:file={seedname}.mat\n"
            f"material:unit-cell-file = {seedname}.UCF\n"
            "#------------------------------------------\n"
        )


def dump_mat(spinham, filename=None, materials=None, nologo=False):
    """
    Write .mat file for |Vampire|_.

    Parameters
    ----------
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian to be saved.
    filename : str, optional
        Name for the .mat file. No extensions is added automatically.
        If not given, the output is returned as a string.
    materials : list of int, optional
        List of materials for the atoms. Has to have the same length as the number of
        magnetic atoms in the ``spinham``. Order is the same as in :py:attr:`.SpinHamiltonian.magnetic_atoms`.
        If not given, each atom will be considered as a separate material. Starting from 0.
    nologo : bool, default False
        Whether to print the logo in the output files.
    """
    if nologo:
        result = []
    else:
        result = [logo(comment=True, date_time=True) + "\n"]
    if materials is not None:
        result.append(f"material:num-materials = {max(materials)+1}\n")
    else:
        result.append(f"material:num-materials = {len(spinham.magnetic_atoms)}\n")
    for i, atom in enumerate(spinham.magnetic_atoms):
        if materials is None or materials[i] not in materials[:i]:
            if materials is not None:
                m_i = materials[i] + 1
            else:
                m_i = i + 1
            result.append("#---------------------------------------------------\n")
            result.append(f"# Material {m_i} \n")
            result.append("#---------------------------------------------------\n")
            result.append(f"material[{m_i}]:material-name = {atom.name}\n")
            result.append(f"material[{m_i}]:material-element = {atom.type}\n")
            result.append(f"material[{m_i}]:atomic-spin-moment={atom.spin}\n")
            result.append(
                f"material[{m_i}]:initial-spin-direction = {atom.spin_direction[0]:.8f},{atom.spin_direction[1]:.8f},{atom.spin_direction[2]:.8f}\n"
            )
    result.append("#---------------------------------------------------\n")

    result = "".join(result)

    if filename is None:
        return "".join(result)

    with open(filename, "w") as file:
        file.write("".join(result))


def dump_ucf(
    spinham,
    filename=None,
    anisotropic=True,
    dmi=True,
    custom_mask=None,
    decimals=5,
    materials=None,
    nologo=False,
):
    """
    Save the Hamiltonian in a Human-readable format.

    Parameters
    ----------
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian to be saved.
    filename : str, optional
        Name for the .UCF file. No extension is added automatically.
        If not given, the output is returned as a string.
    anisotropic : bool, default True
        Whether to output anisotropic exchange.
    dmi : bool, default True
        Whether to output DMI exchange.
    custom_mask : func, optional
        Custom mask for the exchange parameter. Function which take (3,3) numpy:`ndarray`
        as an input and returns (3,3) numpy:`ndarray` as an output. If given, then
        ``anisotropic`` and ``dmi`` parameters are ignored.
    decimals : int, default 4
        Number of decimals to be printed (only for the exchange values).
    materials : list of int, optional
        List of materials for the atoms. Has to have the same length as the number of
        magnetic atoms in the ``spinham``. Order is the same as in :py:attr:`.SpinHamiltonian.magnetic_atoms`.
        If not given, each atom will be considered as a separate material.
    nologo : bool, default False
        Whether to print the logo in the output files.

    Returns
    -------
    content : str
        Content of the .UCF file if ``filename`` is not given.
    """
    original_notation = spinham.notation
    spinham.notation = "vampire"
    if nologo:
        result = []
    else:
        result = [logo(comment=True, date_time=True) + "\n"]
    result.append("# Unit cell size:\n")
    result.append(f"{spinham.a:.8f} {spinham.b:.8f} {spinham.c:.8f}\n")
    result.append("# Unit cell lattice vectors:\n")
    result.append(
        print_2d_array(
            spinham.a1, borders=False, print_result=False, fmt=".8f"
        ).lstrip()
        + "\n"
    )
    result.append(
        print_2d_array(
            spinham.a2, borders=False, print_result=False, fmt=".8f"
        ).lstrip()
        + "\n"
    )
    result.append(
        print_2d_array(
            spinham.a3, borders=False, print_result=False, fmt=".8f"
        ).lstrip()
        + "\n"
    )
    result.append("# Atoms\n")
    result.append(f"{len(spinham.magnetic_atoms)} 1\n")
    atom_indices = {}
    for a_i, atom in enumerate(spinham.magnetic_atoms):
        result.append(
            f"{a_i} {atom.position[0]:.8f} {atom.position[1]:.8f} {atom.position[2]:.8f}"
        )
        if materials is not None:
            result.append(f" {materials[a_i]}")
        else:
            result.append(f" {a_i}")
        result.append("\n")
        atom_indices[atom] = a_i
    result.append("# Interactions\n")
    result.append(f"{len(spinham)} tensorial\n")

    IID = 0
    for atom1, atom2, (i, j, k), J in spinham:
        if custom_mask is not None:
            J = ExchangeParameter(custom_mask(J))
        else:
            if not dmi:
                J = ExchangeParameter(matrix=J.matrix - J.dmi_matrix)
            if not anisotropic:
                J = ExchangeParameter(matrix=J.matrix - J.aniso)
        J = J * meV_TO_J
        fmt = f"{7+decimals}.{decimals}e"
        result.append(
            f"{IID:<2} {atom_indices[atom1]:>2} {atom_indices[atom2]:>2} {i:>2} {j:>2} {k:>2} "
        )
        result.append(f"{J.xx:{fmt}} {J.xy:{fmt}} {J.xz:{fmt}} ")
        result.append(f"{J.yx:{fmt}} {J.yy:{fmt}} {J.yz:{fmt}} ")
        result.append(f"{J.zx:{fmt}} {J.zy:{fmt}} {J.zz:{fmt}}\n")
        IID += 1

    spinham.notation = original_notation

    result = "".join(result)

    if filename is None:
        return result

    with open(filename, "w") as file:
        file.write(result)

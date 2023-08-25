r"""
Input-output for |Vampire|_.
"""

__all__ = ["dump_vampire"]

from radtools.decorate.array import print_2d_array
from radtools.spinham.parameter import ExchangeParameter

meV_TO_J = 1.602176634e-22


def dump_vampire(
    spinham,
    filename=None,
    anisotropic=True,
    dmi=True,
    custom_mask=None,
    decimals=5,
    materials=None,
):
    """
    Save the Hamiltonian in a Human-readable format.

    Parameters
    ----------
    spinham : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian to be saved.
    filename : str, optional
        Name of the file for the Hamiltonian to be saved in.
        If not given, the result will be printed in the console.
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
        magnetic atoms in the ``spinham``. If not given, all atoms will be considered
        as the same material. Order is the same as in :py:attr:`.SpinHamiltonian.magnetic_atoms`.
    """
    original_notation = spinham.notation
    spinham.notation = "vampire"
    result = []
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
            f"{IID:<2} {atom_indices[atom1]:<2} {atom_indices[atom2]:<2} {i:<2} {j:<2} {k:<2} "
        )
        result.append(f"{J.xx:{fmt}} {J.xy:{fmt}} {J.xz:{fmt}} ")
        result.append(f"{J.yx:{fmt}} {J.yy:{fmt}} {J.yz:{fmt}} ")
        result.append(f"{J.zx:{fmt}} {J.zy:{fmt}} {J.zz:{fmt}}\n")
        IID += 1

    spinham.notation = original_notation

    result = "".join(result)
    if filename is not None:
        with open(filename, "w") as file:
            file.write(result)
    else:
        print(result)

r"""
Input-output from TB2J.
"""

import numpy as np

from rad_tools.exchange.model import Bond, ExchangeModel
from rad_tools.routines import absolute_to_relative


def read_exchange_model(filename, quiet=False) -> ExchangeModel:
    r"""
    Read exchange model from TB2J output file.

    In |TB2J|_ exchange Hamiltonian is define in a following notation:

    .. math::
        H = -\sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}){i,j} \hat{\boldsymbol{S}}_j

    where spin vectors :math:`\boldsymbol{S}_i` are normalized to 1
    and double counting is present (both :math:`ij` and :math:`ji` are in the sum).
    :py:class:`.ExchangeModel` can store exchange values in any notation.
    One could check the notation by calling the attribute :py:func:`.ExchangeModel.notation`.
    This function reads and stores exchange parameters in the notation of Hamiltonian
    mentioned above, since for the conversion to the non-normalized spins one need to
    know the values of spin for each atom and those data are not present in the |TB2J|_
    file. One could have the data of spins for each atom, but this silent conversion may
    be confusing for the user, thus it is avoided here.

    Parameters
    ----------
    filename : str
        Path to the TB2J output file.
    quiet : bool, default True
        Whether to suppress output.

    Returns
    -------
    model : :py:class:`.ExchangeModel`
        Exchange model build from TB2J file.
    """

    major_sep = "=" * 90
    minor_sep = "-" * 88
    garbage = str.maketrans(
        {"(": None, ")": None, "[": None, "]": None, ",": None, "'": None}
    )
    # Do not correct spelling, it is taken from TB2J.
    cell_flag = "Cell (Angstrom):"
    atoms_flag = "Atoms:"
    atom_end_flag = "Total"
    exchange_flag = "Exchange:"
    iso_flag = "J_iso:"
    aniso_flag = "J_ani:"
    dmi_flag = "DMI:"

    file = open(filename, "r")
    model = ExchangeModel()
    line = True
    atoms = {}

    # Read everything before exchange
    while line:
        line = file.readline()

        # Read cell
        if line and cell_flag in line:
            model.cell = np.array(
                [
                    list(map(float, file.readline().split())),
                    list(map(float, file.readline().split())),
                    list(map(float, file.readline().split())),
                ]
            )

        # Read atoms
        if line and atoms_flag in line:
            line = file.readline()
            line = file.readline()
            line = file.readline().split()
            while line and atom_end_flag not in line:
                atoms[line[0]] = absolute_to_relative(
                    model.cell, *tuple(map(float, line[1:4]))
                )
                line = file.readline().split()

        # Check if the exchange section is reached
        if line and exchange_flag in line:
            break

    # Read exchange
    while line:
        while line and minor_sep not in line:
            line = file.readline()
        line = file.readline().translate(garbage).split()
        atom1 = line[0]
        atom2 = line[1]
        R = tuple(map(int, line[2:5]))
        distance = float(line[-1])
        iso = None
        aniso = None
        dmi = None
        while line and minor_sep not in line:
            line = file.readline()

            # Read isotropic exchange
            if line and iso_flag in line:
                iso = float(line.split()[-1])

            # Read anisotropic exchange
            if line and aniso_flag in line:
                aniso = np.array(
                    [
                        list(map(float, file.readline().translate(garbage).split())),
                        list(map(float, file.readline().translate(garbage).split())),
                        list(map(float, file.readline().translate(garbage).split())),
                    ]
                )

            # Read DMI
            if line and dmi_flag in line:
                dmi = tuple(map(float, line.translate(garbage).split()[-3:]))

        # Adding info from the exchange block to the ExchangeModel structure
        if atom1 not in model.magnetic_atoms:
            model.add_atom(atom1, *atoms[atom1])
        if atom2 not in model.magnetic_atoms:
            model.add_atom(atom2, *atoms[atom2])
        bond = Bond(iso=iso, aniso=aniso, dmi=dmi)
        model.add_bond(bond, atom1, atom2, R)
        computed_distance = model.get_distance(atom1, atom2, R)
        if abs(computed_distance - distance) > 0.001 and not quiet:
            print(
                f"\nComputed distance is a different from the read one:\n"
                + f"  Computed: {computed_distance:.4f}\n  "
                + f"Read: {distance:.4f}\n"
            )

    # Fill non-magnetic atoms
    for atom in atoms:
        if atom not in model.magnetic_atoms:
            model.nonmagnetic_atoms[atom] = np.array(atoms[atom])

    return model

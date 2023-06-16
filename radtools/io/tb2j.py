r"""
Input-output from |TB2J|_.
"""

import numpy as np

from radtools.crystal.atom import Atom
from radtools.exchange.hamiltonian import ExchangeHamiltonian
from radtools.exchange.parameter import ExchangeParameter


def read_tb2j_model(filename, quiet=True) -> ExchangeHamiltonian:
    r"""
    Read exchange Hamiltonian from |TB2J|_ output file.

    .. versionchanged:: 0.7

    In |TB2J|_ exchange Hamiltonian is define in a following notation:

    .. math::
        H = -\sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \hat{\boldsymbol{S}}_j

    where spin vectors :math:`\boldsymbol{S}_i` are normalized to 1
    and double counting is present (both :math:`ij` and :math:`ji` are in the sum).
    :py:class:`.ExchangeHamiltonian` can store exchange values in any notation.
    One could check the notation by calling the attribute :py:attr:`.ExchangeHamiltonian.notation`.

    This function reads and stores exchange parameters in the notation of Hamiltonian
    mentioned above.

    Distance between the atoms are not read from the atoms,
    but rather computed from the unit cell and atom positions
    (see :py:meth:`.Crystal.get_distance`). If the distance read from the file
    is different from the computed one and ``quiet=False``, the warning is printed.
    It is usual for the computed and read distances to differ in the last digits, since
    the unit cell is provided with the same precision as the distance in
    |TB2J|_ output file.


    Parameters
    ----------
    filename : str
        Path to the |TB2J|_ output file.
    quiet : bool, default True
        Whether to suppress output.

    Returns
    -------
    model : :py:class:`.ExchangeHamiltonian`
        Exchange Hamiltonian build from |TB2J|_ file. With notation set to "TB2J".
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
    model = ExchangeHamiltonian(notation="TB2J")
    line = True

    # Read everything before exchange
    while line:
        line = file.readline()

        # Read cell
        if line and cell_flag in line:
            model.crystal.lattice.cell = np.array(
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
            i = 1
            while line and atom_end_flag not in line:
                try:
                    # Slicing is not used intentionally.
                    magmom = tuple(map(float, [line[5], line[6], line[7]]))
                except IndexError:
                    magmom = None
                try:
                    charge = float(line[4])
                except IndexError:
                    charge = None
                model.add_atom(
                    Atom(
                        name=line[0],
                        position=tuple(map(float, line[1:4])),
                        index=i,
                        magmom=magmom,
                        charge=charge,
                    )
                )
                line = file.readline().split()
                i += 1

        # Check if the exchange section is reached
        if line and exchange_flag in line:
            break

    # Identify lattice type
    model.crystal.identify()

    # Read exchange
    while line:
        while line and minor_sep not in line:
            line = file.readline()
        line = file.readline().translate(garbage).split()
        # atom1 = model.crystal.get_atom(line[0])
        # atom2 = model.crystal.get_atom(line[1])
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

        # Adding info from the exchange block to the ExchangeHamiltonian structure
        J = ExchangeParameter(iso=iso, aniso=aniso, dmi=dmi)
        model.add_bond(J, atom1, atom2, R)
        computed_distance = model.get_distance(atom1, atom2, R)
        if abs(computed_distance - distance) > 0.001 and not quiet:
            print(
                f"\nComputed distance is a different from the read one:\n"
                + f"  Computed: {computed_distance:.4f}\n  "
                + f"Read: {distance:.4f}\n"
            )

    return model

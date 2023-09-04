r"""
Input-output from |TB2J|_.
"""
__all__ = ["load_tb2j_model", "read_tb2j_model"]

import numpy as np

from radtools.crystal.constants import REL_TOL
from radtools.spinham.hamiltonian import SpinHamiltonian


def load_tb2j_model(filename, quiet=True, bravais_type=None) -> SpinHamiltonian:
    r"""
    Read spin Hamiltonian from |TB2J|_ output file.

    .. versionchanged:: 0.8.1 Renamed from ``read_tb2j_model``


    In |TB2J|_ spin Hamiltonian is define in a following notation:

    .. math::
        H = -\sum_{i,j} \hat{\boldsymbol{S}}_i \cdot \boldsymbol{J}_{i,j} \hat{\boldsymbol{S}}_j

    where spin vectors :math:`\boldsymbol{S}_i` are normalized to 1
    and double counting is present (both :math:`ij` and :math:`ji` are in the sum).
    :py:class:`.SpinHamiltonian` can store exchange values in any notation.
    One could check the notation by calling the attribute :py:attr:`.SpinHamiltonian.notation`.

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
    bravais_type : str, default None
        Expected bravais lattice type. If ``None``, the lattice type is identified
        automatically. See :py:meth:`.Crystal.identify` for more details.
        The bravais lattice type is reached by reducing the accuracy for the :ref:`library_lepage`.
        If the desired lattice type is not reached, the error is raised.

    Returns
    -------
    model : :py:class:`.SpinHamiltonian`
        Spin Hamiltonian build from |TB2J|_ file. With notation set to "TB2J".

    Raises
    ------
    ValueError
        If ``bravais_type`` is provided and can not be reached.
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
    model = SpinHamiltonian(notation="TB2J")
    line = True

    # Read everything before exchange
    while line:
        line = file.readline()

        # Read cell and define eps_rel
        if line and cell_flag in line:
            a = file.readline().split()
            b = file.readline().split()
            c = file.readline().split()

            try:
                n_a = min(*tuple(map(lambda x: len(x.split(".")[1]), a)))
                n_b = min(*tuple(map(lambda x: len(x.split(".")[1]), b)))
                n_c = min(*tuple(map(lambda x: len(x.split(".")[1]), c)))
                model.eps_rel = 10 ** (-min(n_a, n_b, n_c))
            except:
                model.eps_rel = REL_TOL

            model.cell = np.array(
                [
                    list(map(float, a)),
                    list(map(float, b)),
                    list(map(float, c)),
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
                    name=line[0],
                    position=tuple(map(float, line[1:4])),
                    index=i,
                    magmom=magmom,
                    charge=charge,
                    relative=False,
                )
                line = file.readline().split()
                i += 1

        # Check if the exchange section is reached
        if line and exchange_flag in line:
            break

    # Try to reach desired bravais type, if any
    if bravais_type is not None:
        while model.eps_rel < 1 and model.type() != bravais_type:
            model.eps_rel *= 10
        if model.type() != bravais_type:
            raise ValueError(f"Bravais type {bravais_type} could not be reached.")

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

        # Adding info from the exchange block to the SpinHamiltonian structure
        model.add_bond(atom1, atom2, R, iso=iso, aniso=aniso, dmi=dmi)
        computed_distance = model.get_distance(atom1, atom2, R)
        if abs(computed_distance - distance) > 0.001 and not quiet:
            print(
                f"\nComputed distance is a different from the read one:\n"
                + f"  Computed: {computed_distance:.4f}\n  "
                + f"Read: {distance:.4f}\n"
            )

    return model


# For backward compatibility
read_tb2j_model = load_tb2j_model

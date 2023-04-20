r"""
DOS
"""
import re
from os import walk
from os.path import abspath, join

import matplotlib.pyplot as plt
import numpy as np

from rad_tools.dos.pdos import PDOS


class DOS:
    r"""
    Analyse the folder with the Quantum Espresso results and store all information.

    Parameters
    ----------
    seedname : str
        Seedname for the Quantum Espresso output files.
    input_folder : str
        Path to the folder with Quantum Espresso output files.

    Attributes
    ----------
    seedname : str
        Seedname for the Quantum Espresso output files.
    energy : array
        Energy of the DOS/PDOS
    window : tuple, default (0, -1)
    efermi : float, defaul 0
        Fermi energy.
    k_points : array
    total_dos : array
    total_pdos : array
    case : int
        There are 4 cases:

        #. Collinear
        #. Collinear, spin-polarized
        #. Non-collinear, non-spin-orbit
        #. Non-collinear, spin-orbit

        Headers of "seedname".pdos_tot with according to the
        :projwfc:`QE <>`:

        * Collinear:
            E DOS(E) PDOS(E)

        * Collinear, spin-polarized:
            E DOSup(E) DOSdw(E) PDOSup(E) PDOSdw(E)

        * Non-collinear, non-spin-orbit:
            E DOS(E) PDOSup(E) PDOSdw(E)

        * Non-collinear, spin-orbit:
            E DOS(E) PDOS(E)
    filenames : list
    """

    def __init__(self, seedname: str, input_path: str, efermi=0, window=None) -> None:
        self.efermi = 0
        self.seedname = seedname
        self._input_path = input_path
        self.k_resolved = False
        self.case = None
        self._filenames = None
        self._energy = None
        self.nkpoints = None
        self.nepoints = None
        self.window = window
        self._detect_case()
        self._exctract_energy()
        self._atoms = {}
        self._wfcs = {}
        self._decompose_filenames()

    def __iter__(self):
        return DOSIterator(self)

    def __contains__(self, item):
        atom, atom_number, wfc, wfc_number = item
        return (
            atom in self.atoms
            and atom_number in self.atom_numbers(atom)
            and (wfc, wfc_number) in self.wfcs(atom, atom_number)
        )

    @property
    def case_name(self):
        r"""
        Human-readable name of the case.

        Returns
        -------
        case_name : str
            Name of the case.
        """
        case_names = [
            "collinear, spin-unpolarized",
            "collinear, spin-polarized",
            "non-collinear, non spin-orbit",
            "non-collinear, spin-orbit",
        ]
        return case_names[self.case - 1]

    @property
    def filenames(self):
        r"""
        List of filenames (without "seedname".pdos_tot).
        """
        if self._filenames is not None:
            return self._filenames

        pattern = ".pdos_atm#[0-9]*\\([a-zA-Z]*\\)_wfc#[0-9]*\\([spdf_0-9j.]*\\)"

        # Get list of files in the folder
        all_filenames = []
        for dirpath, dirnames, filenames in walk(self._input_path):
            all_filenames.extend(filenames)
            break

        self._filenames = []
        # Search
        for filename in all_filenames:
            if re.fullmatch(
                f"{re.escape(self.seedname)}{pattern}",
                filename,
            ):
                self._filenames.append(filename)

        return self._filenames

    @property
    def energy(self):
        return self._energy[self.window[0] : self.window[1]]

    def total_dos(self, squeeze=False):
        r"""
        Total density of states.

        Density of states, computed directly from plane-wave basis set.
        """
        dos = np.loadtxt(
            join(self._input_path, f"{self.seedname}.pdos_tot"), skiprows=1
        ).T

        if self.case == 2:
            if self.k_resolved:
                if squeeze:
                    return (
                        np.sum(
                            dos[2:4].reshape((2, self.nkpoints, self.nepoints))[
                                :, :, self.window[0] : self.window[1]
                            ],
                            axis=1,
                        )
                        / self.nkpoints
                    )
                return dos[2:4].reshape((2, self.nkpoints, self.nepoints))[
                    :, :, self.window[0] : self.window[1]
                ]
            return dos[1:3][:, self.window[0] : self.window[1]]

        if self.k_resolved:
            if squeeze:
                return (
                    np.sum(
                        dos[2].reshape((self.nkpoints, self.nepoints))[
                            :, self.window[0] : self.window[1]
                        ],
                        axis=0,
                    )
                    / self.nkpoints
                )
            return dos[2].reshape((self.nkpoints, self.nepoints))[
                :, self.window[0] : self.window[1]
            ]
        return dos[1][self.window[0] : self.window[1]]

    def total_pdos(self, squeeze=False):
        r"""
        Total density of states.

        Density of states, computed directly from plane-wave basis set.
        """
        dos = np.loadtxt(f"{self._input_path}/{self.seedname}.pdos_tot", skiprows=1).T

        if self.case == 2:
            if self.k_resolved:
                if squeeze:
                    return (
                        np.sum(
                            dos[4:6].reshape((2, self.nkpoints, self.nepoints))[
                                :, :, self.window[0] : self.window[1]
                            ],
                            axis=1,
                        )
                        / self.nkpoints
                    )
                return dos[4:6].reshape((2, self.nkpoints, self.nepoints))[
                    :, :, self.window[0] : self.window[1]
                ]
            return dos[3:5][:, self.window[0] : self.window[1]]
        elif self.case == 3:
            if self.k_resolved:
                if squeeze:
                    return (
                        np.sum(
                            dos[3:5].reshape((2, self.nkpoints, self.nepoints))[
                                :, :, self.window[0] : self.window[1]
                            ],
                            axis=1,
                        )
                        / self.nkpoints
                    )
                return dos[3:5].reshape((2, self.nkpoints, self.nepoints))[
                    :, :, self.window[0] : self.window[1]
                ]
            return dos[2:4][:, self.window[0] : self.window[1]]
        if self.k_resolved:
            if squeeze:
                return (
                    np.sum(
                        dos[3].reshape((self.nkpoints, self.nepoints))[
                            :, self.window[0] : self.window[1]
                        ],
                        axis=0,
                    )
                    / self.nkpoints
                )

            return dos[3].reshape((self.nkpoints, self.nepoints))[
                :, self.window[0] : self.window[1]
            ]
        return dos[2][self.window[0] : self.window[1]]

    @property
    def atoms(self):
        r"""
        List of atom's types.

        Returns
        -------
        atoms : list
            List of atom's types labels

            .. code-block:: python

                atoms = [atom1, atom2, ...]
        """

        return [atom for atom in self._atoms]

    def atom_numbers(self, atom: str):
        r"""
        List of atom's number for particular atom type.

        Parameters
        ----------
        atom : str
            Label of atom.

        Returns
        -------
        numbers : list

            .. code-block:: python

                numbers = [n_1, n_2, n_3, ...]
        """

        return self._atoms[atom]

    def wfcs(self, atom: str, atom_number: int):
        r"""
        Return list of wave function labels for particular atom.

        Parameters
        ----------
        atom : str
            Label of an atom.
        atom_number : int
            Number of an atom.

        Returns
        -------
        wfcs : list
            List of wave function labels and numbers

            .. code-block:: python

                wfcs = [(wfc_label, wfc_number), ...]

        """

        return self._wfcs[(atom, atom_number)]

    def ldos(
        self, atoms=None, atom_numbers=None, wfcs=None, wfc_numbers=None, squeeze=False
    ):
        r"""
        Local projected density of states for a particular atom.
        """
        if (
            atoms is None
            and atom_numbers is None
            and wfcs is None
            and wfc_numbers is None
        ):
            return self.total_pdos(squeeze=squeeze)

        def load_ldos(atom, atom_number, wfc, wfc_number):
            path = join(
                self._input_path,
                f"{self.seedname}.pdos_atm#{atom_number}({atom})_wfc#{wfc_number}({wfc})",
            )
            ldos = np.loadtxt(path, skiprows=1).T
            if self.case in [1, 4]:
                if self.k_resolved:
                    ldos = ldos[2]
                    ldos = ldos.reshape((self.nkpoints, self.nepoints))[
                        :, self.window[0] : self.window[1]
                    ]
                else:
                    ldos = ldos[1][self.window[0] : self.window[1]]
            else:
                if self.k_resolved:
                    ldos = ldos[2:4]
                    ldos = ldos.reshape((2, self.nkpoints, self.nepoints))[
                        :, :, self.window[0] : self.window[1]
                    ]
                else:
                    ldos = ldos[1:3][:, self.window[0] : self.window[1]]
            return ldos

        for atom, atom_number, wfc, wfc_number in self:
            ldos = np.zeros(load_ldos(atom, atom_number, wfc, wfc_number), dtype=float)
            break

        counter = 0
        for atom, atom_number, wfc, wfc_number in self:
            if (
                (atoms is None or atom == atoms or atom in atoms)
                and (
                    atom_numbers is None
                    or atom_number == atom_numbers
                    or atom_number in atom_numbers
                )
                and (wfcs is None or wfc == wfcs or wfc in wfcs)
                and (
                    wfc_numbers is None
                    or wfc_number == wfc_numbers
                    or wfc in wfc_numbers
                )
            ):
                ldos += load_ldos(atom, atom_number, wfc, wfc_number)
                counter += 1
        if counter:
            ldos /= counter

        if self.k_resolved and squeeze:
            if self.case in [1, 4]:
                return np.sum(ldos, axis=0)
            else:
                return np.sum(ldos, axis=1)

        return ldos

    def pdos(
        self, atom, wfc, atom_numbers=None, wfc_numbers=None, squeeze=False
    ) -> PDOS:
        r"""
        Projected density of states for a particular atom.

        If ``atom_number`` = 0 summ all pdos over the atoms of the same type.
        """

        if atom_numbers is None:
            atom_numbers = self.atom_numbers(atom)
        elif isinstance(atom_numbers, int):
            if atom_numbers in self.atom_numbers(atom):
                atom_numbers = [atom_numbers]
            else:
                atom_numbers = []

        if len(atom_numbers) > 0:
            avail_wfc_numbers = [i[1] for i in self.wfcs(atom, atom_numbers[0])]

            if wfc_numbers is None:
                wfc_numbers = avail_wfc_numbers
            elif isinstance(wfc_numbers, int):
                if wfc_numbers in avail_wfc_numbers:
                    wfc_numbers = [wfc_numbers]
                else:
                    wfc_numbers = []

        def load_pdos(atom, atom_number, wfc, wfc_number):
            if (wfc, wfc_number) not in self.wfcs(atom, atom_number):
                raise KeyError(f"No {wfc}_{wfc_number} in {atom}_{atom_number}")
            path = join(
                self._input_path,
                f"{self.seedname}.pdos_atm#{atom_number}({atom})_wfc#{wfc_number}({wfc})",
            )
            pdos = np.loadtxt(path, skiprows=1).T
            if self.k_resolved:
                pdos = pdos[2:]
                pdos = pdos.reshape((pdos.shape[0], self.nkpoints, self.nepoints))[
                    :, :, self.window[0] : self.window[1]
                ]
            else:
                pdos = pdos[1:][:, self.window[0] : self.window[1]]
            return pdos

        counter = 0
        for atom_number in atom_numbers:
            for wfc_number in wfc_numbers:
                pdos = np.zeros(
                    load_pdos(atom, atom_number, wfc, wfc_number).shape, dtype=float
                )
                break
            break

        for atom_number in atom_numbers:
            for wfc_number in wfc_numbers:
                pdos += load_pdos(atom, atom_number, wfc, wfc_number)
                counter += 1
        if counter:
            pdos /= counter

        if self.k_resolved and squeeze:
            pdos = np.sum(pdos, axis=1)
        return PDOS(
            energy=self.energy,
            pdos=pdos,
            case=self.case,
            projector_family=wfc,
            k_resolved=self.k_resolved and not squeeze,
        )

    def _detect_case(self):
        r"""
        Detects case of the DOS calculations.
        """

        with open(f"{self._input_path}/{self.seedname}.pdos_tot") as file:
            header = file.readline().lower().split()

        # Check for k-resolved
        if "ik" in header:
            self.k_resolved = True

        # Detect the case
        if "dos(e)" in header and "pdos(e)" in header:
            self.case = 1
        if "dos(e)" in header and "pdosup(e)" in header and "pdosdw(e)" in header:
            self.case = 3
        if (
            "dosup(e)" in header
            and "dosdw(e)" in header
            and "pdosup(e)" in header
            and "pdosdw(e)" in header
        ):
            self.case = 2
        if self.case == 1:
            pattern = ".pdos_atm#[0-9]*\\([a-zA-Z]*\\)_wfc#[0-9]*\\([spdf.]*\\)"
            filename = self.filenames[0]
            if not re.fullmatch(
                f"{re.escape(self.seedname)}{pattern}",
                filename,
            ):
                self.case = 4
        if self.case is None:
            raise RuntimeError(
                "Unable to detect case, analysed header:\n" + f"{header}"
            )

    def _exctract_energy(self):
        dos = np.loadtxt(f"{self._input_path}/{self.seedname}.pdos_tot", skiprows=1).T
        if self.k_resolved:
            self.nkpoints = int(dos[0][-1])
            self.nepoints = len(dos[0]) // self.nkpoints
            self._energy = dos[1][0 : self.nepoints] - self.efermi
        else:
            self._energy = dos[0] - self.efermi

        if self.window is None:
            self.window = (0, -1)
        else:
            i, j = 0, 0
            for e_i, e in enumerate(self._energy):
                if e <= self.window[0]:
                    i = e_i
                if e <= self.window[1]:
                    j = e_i
            self.window = (i, j + 1)

    def _decompose_filenames(self):
        r"""
        Decompose filenames and extract information about atoms and projections.

        Parameters
        ----------
        filenames : str
            List of filenames with projected DOS (not pdos_total).

        Note
        ----
        atoms : dist
            Dictionary with the atom labels and their numbers.

            .. code-block:: python

                atoms = {atom1: [n_1, n_2, n_3, ...]}
        wfcs : dict
            Dictionary of projection functions and their numbers.

            .. code-block:: python

                wfcs = {(atom1, atom1_number): (wfc_label, wfc_number), ...}
        """

        for filename in self.filenames:
            # Detect names and numbers
            atom_number = int(filename.split(".pdos_atm#")[1].split("(")[0])
            atom_symbol = filename.split("(")[1].split(")")[0]
            wfc_number = int(filename.split(")_wfc#")[1].split("(")[0])
            wfc_symbol = filename.split(")_wfc#")[1].split("(")[1].split(")")[0]

            if atom_symbol not in self._atoms:
                self._atoms[atom_symbol] = []
            self._atoms[atom_symbol].append(atom_number)

            if (atom_symbol, atom_number) not in self._wfcs:
                self._wfcs[(atom_symbol, atom_number)] = []
            self._wfcs[(atom_symbol, atom_number)].append((wfc_symbol, wfc_number))

        # Sort entries
        for atom in self._atoms:
            self._atoms[atom] = list(set(self._atoms[atom]))
            self._atoms[atom].sort()
        for key in self._wfcs:
            self._wfcs[key].sort(key=lambda x: x[1])

    def plot_pdos_tot(self, output_name, interactive=False, ylim=None):
        ax = plt.subplots(figsize=(8, 4))[1]

        ax.set_xlabel("Energy, eV", fontsize=18)
        ax.set_ylabel("DOS, states/eV", fontsize=18)
        if self.efermi != 0:
            ax.set_title(f"DOS vs PDOS (0 is Fermi energy)", fontsize=18)
        else:
            ax.set_title(f"DOS vs PDOS (0 is 0)", fontsize=18)
        ax.vlines(
            0,
            0,
            1,
            transform=ax.get_xaxis_transform(),
            color="grey",
            linewidths=0.5,
            linestyles="dashed",
        )

        if self.case in [1, 4]:
            total_dos = self.total_dos(squeeze=True)
            total_pdos = self.total_pdos(squeeze=True)

            ax.fill_between(
                self.energy,
                0,
                total_dos,
                lw=0,
                color="grey",
                alpha=0.3,
                label="DOS",
            )
            ax.plot(
                self.energy,
                total_pdos,
                "-",
                lw=1,
                color="black",
                alpha=0.7,
                label="PDOS",
            )
            ncol = 1
        if self.case == 2:
            total_dos = self.total_dos(squeeze=True)
            total_pdos = self.total_pdos(squeeze=True)

            ax.fill_between(
                self.energy,
                0,
                total_dos[0],
                lw=0,
                color="blue",
                alpha=0.3,
                label="DOS (up)",
            )
            ax.fill_between(
                self.energy,
                0,
                -total_dos[1],
                lw=0,
                color="red",
                alpha=0.3,
                label="DOS (down)",
            )

            ax.plot(
                self.energy,
                total_pdos[0],
                "-",
                lw=1,
                color="blue",
                alpha=0.7,
                label="PDOS (up)",
            )
            ax.plot(
                self.energy,
                -total_pdos[1],
                "-",
                lw=1,
                color="red",
                alpha=0.7,
                label="PDOS (down)",
            )
            ncol = 2
        if self.case == 3:
            total_dos = self.total_dos
            total_pdos = self.total_pdos
            if self.k_resolved:
                total_dos = np.sum(total_dos, axis=0) / self.nkpoints
                total_pdos = np.sum(total_pdos, axis=1) / self.nkpoints

            ax.fill_between(
                self.energy, 0, total_dos, lw=0, color="grey", alpha=0.3, label="DOS"
            )
            ax.fill_between(
                self.energy, -total_dos, 0, lw=0, color="grey", alpha=0.3, label="-DOS"
            )

            ax.plot(
                self.energy,
                total_pdos[0],
                "-",
                lw=1,
                color="blue",
                alpha=0.7,
                label="PDOS (up)",
            )
            ax.plot(
                self.energy,
                -total_pdos[1],
                "-",
                lw=1,
                color="red",
                alpha=0.7,
                label="PDOS (down)",
            )
            ncol = 1

        ax.set_xlim(np.amin(self.energy), np.amax(self.energy))
        if ylim is not None:
            ax.set_ylim(*ylim)
        ax.legend(loc="best", ncol=ncol, draggable=True)

        if interactive:
            plt.show()
        else:
            png_path = f"{output_name}.png"
            plt.savefig(png_path, dpi=600, bbox_inches="tight")
            print(f"Total DOS plot is in {abspath(png_path)}")
        plt.close()


class DOSIterator:
    def __init__(self, dos: DOS) -> None:
        self.list = []
        for atom in dos.atoms:
            for atom_number in dos.atom_numbers(atom):
                for wfc, wfc_number in dos.wfcs(atom, atom_number):
                    self.list.append([atom, atom_number, wfc, wfc_number])
        self._index = 0

    def __next__(self) -> str:
        if self._index < len(self.list):
            result = self.list[self._index]
            self._index += 1
            return result
        raise StopIteration

    def __iter__(self):
        return self

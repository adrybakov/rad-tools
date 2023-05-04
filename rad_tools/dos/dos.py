r"""
DOS
"""
import re
from os import walk
from os.path import abspath, join

import matplotlib.pyplot as plt
import numpy as np

from rad_tools.dos.pdos import PDOSQE


class DOSQE:
    r"""
    Analyse the folder with the Quantum Espresso results and store all information.

    Parameters
    ----------
    seedname : str
        Seedname for the Quantum Espresso output files.
    input_folder : str
        Path to the folder with Quantum Espresso output files.
    energy_window : tuple, default None
        Energy limits, necessary for correct plotting.

    Attributes
    ----------
    seedname : str
        Seedname for the Quantum Espresso output files.
    energy : array
        Energy values of the DOS/PDOS.
    k_resolved : bool
    case : int
    casename : str
        Human readable name of the case base on ``case``.
    nkpoints : int
        Number of k points for k resolved DOS. 1 if DOS is not k resolved.
    nepoints : int
        Number of energy points.
    filenames : list
    atoms : list
    """

    def __init__(self, seedname: str, input_path: str, energy_window=None) -> None:
        self.seedname = seedname
        self._input_path = input_path

        # Case and whether DOS is k resolved
        self._k_resolved = None
        self._case = None

        # Detect number of energy and k points and read energy values
        self.nkpoints = None
        self.nepoints = None
        self.energy = None
        self._extract_energy()

        self._filenames = None

        self._atoms = {}
        self._wfcs = {}
        self._decompose_filenames()

        if energy_window is None:
            self.energy_window = (0, -1)
        else:
            i, j = 0, 0
            for i_e, e in enumerate(self.energy):
                if e <= energy_window[0]:
                    i += 1
                if e <= energy_window[1]:
                    j += 1
            self.energy_window = (i, min(j + 1, i_e))
        self.energy = self.energy[self.energy_window[0] : self.energy_window[1]]

    @property
    def case(self):
        r"""
        Detects case of the DOS calculations.

        Detects the case when called for the first time.

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
        """

        if self._case is None:
            # Detect the case
            with open(join(f"{self._input_path}", f"{self.seedname}.pdos_tot")) as file:
                header = file.readline().lower().split()

            if "dos(e)" in header and "pdos(e)" in header:
                self._case = 1
            if "dos(e)" in header and "pdosup(e)" in header and "pdosdw(e)" in header:
                self._case = 3
            if (
                "dosup(e)" in header
                and "dosdw(e)" in header
                and "pdosup(e)" in header
                and "pdosdw(e)" in header
            ):
                self._case = 2
            if self._case == 1:
                pattern = ".pdos_atm#[0-9]*\\([a-zA-Z]*\\)_wfc#[0-9]*\\([spdf.]*\\)"
                filename = self.filenames[0]
                if not re.fullmatch(
                    f"{re.escape(self.seedname)}{pattern}",
                    filename,
                ):
                    self._case = 4
            if self._case is None:
                raise RuntimeError(
                    "Unable to detect case, analysed header:\n" + f"{header}"
                )

        return self._case

    @property
    def casename(self):
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
    def k_resolved(self):
        r"""
        Whether DOS is k-resolved.
        """

        if self._k_resolved is None:
            # Check for k-resolved
            with open(join(f"{self._input_path}", f"{self.seedname}.pdos_tot")) as file:
                self._k_resolved = "ik" in file.readline().lower().split()

        return self._k_resolved

    def _extract_energy(self):
        dos = np.loadtxt(
            join(f"{self._input_path}", f"{self.seedname}.pdos_tot"), skiprows=1
        ).T
        if self.k_resolved:
            self.nkpoints = int(dos[0][-1])
            self.nepoints = len(dos[0]) // self.nkpoints
            self.energy = dos[1][0 : self.nepoints]
        else:
            self.nkpoints = 1
            self.nepoints = len(dos[0])
            self.energy = dos[0]

    @property
    def filenames(self):
        r"""
        List of filenames (without "seedname".pdos_tot).
        """
        if self._filenames is None:
            pattern = ".pdos_atm#[0-9]*\\([a-zA-Z0-9]*\\)_wfc#[0-9]*\\([spdf_0-9j.]*\\)"

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

    def _decompose_filenames(self):
        r"""
        Decompose filenames and extract information about atoms and projectors.

        Note
        ----
        atoms : dist
            Dictionary with the atom labels and their numbers.

            .. code-block:: python

                atoms = {atom1: [n_1, n_2, n_3, ...]}
        wfcs : dict
            Dictionary of projectors functions and their numbers.

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

    def __iter__(self):
        return DOSIterator(self)

    def __contains__(self, item):
        atom, atom_number, wfc, wfc_number = item
        return (
            atom in self.atoms
            and atom_number in self.atom_numbers(atom)
            and (wfc, wfc_number) in self.wfcs(atom, atom_number)
        )

    def total_dos(self, squeeze=False):
        r"""
        Total density of states.

        Density of states, computed directly from plane-wave basis set.

        If DOS is k resolved, then sum over k points.

        Parameters
        ----------
        squeeze : bool, default False
            Whenever to sum over k points. Ignored if DOS is not k resolved.

        Returns
        -------
        total_dos : array
            Shape :math:`(2, n_e)` if DOS is collinear, spin-polarized;
            :math:`(n_e)` otherwise.
        """

        dos = np.loadtxt(
            join(self._input_path, f"{self.seedname}.pdos_tot"), skiprows=1
        ).T

        if self.case == 2:
            if self.k_resolved:
                if squeeze:
                    return (
                        np.sum(
                            dos[2:4].reshape((2, self.nkpoints, self.nepoints)), axis=1
                        )
                        / self.nkpoints
                    )[:, self.energy_window[0] : self.energy_window[1]]
                return dos[2:4].reshape((2, self.nkpoints, self.nepoints))[
                    :, :, self.energy_window[0] : self.energy_window[1]
                ]
            return dos[1:3][:, self.energy_window[0] : self.energy_window[1]]

        if self.k_resolved:
            if squeeze:
                return (
                    np.sum(dos[2].reshape((self.nkpoints, self.nepoints)), axis=0)
                    / self.nkpoints
                )[self.energy_window[0] : self.energy_window[1]]
            return dos[2].reshape((self.nkpoints, self.nepoints))[
                :, self.energy_window[0] : self.energy_window[1]
            ]
        return dos[1][self.energy_window[0] : self.energy_window[1]]

    def total_pdos(self, squeeze=False):
        r"""
        Total partial density of states.

        If DOS is k resolved, then sum over k points.

        Parameters
        ----------
        squeeze : bool, default False
            Whenever to sum over k points. Ignored if DOS is not k resolved.

        Returns
        -------
        total_pdos : array
            Shape :math:`(2, n_e)` if PDOS is collinear, spin-polarized;
            :math:`(n_e)` otherwise.
        """

        dos = np.loadtxt(
            join(f"{self._input_path}", f"{self.seedname}.pdos_tot"), skiprows=1
        ).T

        if self.case == 2:
            if self.k_resolved:
                if squeeze:
                    return (
                        np.sum(
                            dos[4:6].reshape((2, self.nkpoints, self.nepoints)), axis=1
                        )
                        / self.nkpoints
                    )[:, self.energy_window[0] : self.energy_window[1]]
                return dos[4:6].reshape((2, self.nkpoints, self.nepoints))[
                    :, :, self.energy_window[0] : self.energy_window[1]
                ]
            return dos[3:5][:, self.energy_window[0] : self.energy_window[1]]

        elif self.case == 3:
            if self.k_resolved:
                if squeeze:
                    return (
                        np.sum(
                            dos[3:5].reshape((2, self.nkpoints, self.nepoints)), axis=1
                        )
                        / self.nkpoints
                    )[:, self.energy_window[0] : self.energy_window[1]]
                return dos[3:5].reshape((2, self.nkpoints, self.nepoints))[
                    :, :, self.energy_window[0] : self.energy_window[1]
                ]
            return dos[2:4][:, self.energy_window[0] : self.energy_window[1]]

        if self.k_resolved:
            if squeeze:
                return (
                    np.sum(dos[3].reshape((self.nkpoints, self.nepoints)), axis=0)
                    / self.nkpoints
                )[self.energy_window[0] : self.energy_window[1]]
            return dos[3].reshape((self.nkpoints, self.nepoints))[
                :, self.energy_window[0] : self.energy_window[1]
            ]
        return dos[2][self.energy_window[0] : self.energy_window[1]]

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

    def wfcs(self, atom: str, atom_number: int = None):
        r"""
        Return list of wave function labels for particular atom.

        Parameters
        ----------
        atom : str
            Label of an atom.
        atom_number : int, default None
            Number of an atom. If ``None`` then return wfc and wfc numbers of first atom.

        Returns
        -------
        wfcs : list
            List of wave function labels and numbers

            .. code-block:: python

                wfcs = [(wfc_label, wfc_number), ...]

        """
        if atom_number is None:
            atom_number = self.atom_numbers(atom)[0]
        return self._wfcs[(atom, atom_number)]

    def pdos(self, atom, wfc, wfc_number, atom_numbers=None) -> PDOSQE:
        r"""
        Projected density of states for a particular atom.

        Parameters
        ----------
        atom : str
            Name of the atom type.
        wfc : str
            Name of the projector wave function.
        wfc_numbers : int
            Number of wave function projector.
        atom_numbers : list or int, default None
            If ``None``, then PDOS summed over all atom numbers for ``atom``.

        Returns
        -------
        pdos : :py:class:`PDOSQE`
            Partial density of states, orbital-resolved.
        """

        if atom_numbers is None:
            atom_numbers = self.atom_numbers(atom)
        elif isinstance(atom_numbers, int):
            if atom_numbers in self.atom_numbers(atom):
                atom_numbers = [atom_numbers]
            else:
                raise ValueError(f"There is no {atom}{atom_numbers} in PDOS.")

        for i, atom_number in enumerate(atom_numbers):
            path = join(
                self._input_path,
                f"{self.seedname}.pdos_atm#{atom_number}({atom})_wfc#{wfc_number}({wfc})",
            )
            if i == 0:
                pdos = np.loadtxt(path, skiprows=1).T
            else:
                pdos += np.loadtxt(path, skiprows=1).T
        if self.case in [2, 3]:
            if self.k_resolved:
                ldos = pdos[2:4].reshape(2, self.nkpoints, self.nepoints)[
                    :, :, self.energy_window[0] : self.energy_window[1]
                ]
                pdos = pdos[4:].reshape(
                    (pdos.shape[0] - 4) // 2, 2, self.nkpoints, self.nepoints
                )[:, :, :, self.energy_window[0] : self.energy_window[1]]
            else:
                ldos = pdos[1:3][:, self.energy_window[0] : self.energy_window[1]]
                pdos = pdos[3:].reshape((pdos.shape[0] - 3) // 2, 2, self.nepoints)[
                    :, :, self.energy_window[0] : self.energy_window[1]
                ]
        else:
            if self.k_resolved:
                ldos = pdos[2].reshape(self.nkpoints, self.nepoints)[
                    :, self.energy_window[0] : self.energy_window[1]
                ]
                pdos = pdos[3:].reshape(
                    pdos.shape[0] - 3, self.nkpoints, self.nepoints
                )[:, :, self.energy_window[0] : self.energy_window[1]]
            else:
                ldos = pdos[1][self.energy_window[0] : self.energy_window[1]]
                pdos = pdos[2:][:, self.energy_window[0] : self.energy_window[1]]

        return PDOSQE(
            energy=self.energy,
            pdos=pdos,
            projectors_group=wfc,
            ldos=ldos,
            spin_pol=self.case in [2, 3],
        )

    def plot_pdos_tot(
        self,
        output_name,
        interactive=False,
        efermi=0,
        xlim=None,
        ylim=None,
        save_pickle=False,
    ):
        r"""
        Plot total DOS vs total PDOS.

        Parameters
        ----------
        output_name : str
            Name of the output file.
        interactive : bool, default False
            Whether to plot in interactive matplotlib window.
        efermi : float, default 0
            Fermi energy. Zero of energy scale is shifted to it.
        xlim : tuple, default None
            Limits for energy scale.
        ylim : tuple, default None
            Limits for dos scale.
        save_pickle : bool, default False
            Whether to save figure as a .pickle file.
            Helps for custom modification of particular figures.
        """
        fig, ax = plt.subplots(figsize=(8, 4))

        ax.set_xlabel("Energy, eV", fontsize=18)
        ax.set_ylabel("DOS, states/eV", fontsize=18)
        if xlim is None:
            xlim = (np.amin(self.energy), np.amax(self.energy))
        ax.set_xlim(*tuple(xlim))
        if ylim is not None:
            ax.set_ylim(*tuple(ylim))
        if efermi != 0:
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
            ax.fill_between(
                self.energy,
                0,
                self.total_dos(squeeze=True),
                lw=0,
                color="grey",
                alpha=0.3,
                label="DOS",
            )
            ax.plot(
                self.energy,
                self.total_pdos(squeeze=True),
                "-",
                lw=1,
                color="black",
                alpha=0.7,
                label="PDOS",
            )
            ncol = 1
        if self.case == 2:
            ax.fill_between(
                self.energy,
                0,
                self.total_dos(squeeze=True)[0],
                lw=0,
                color="blue",
                alpha=0.3,
                label="DOS (up)",
            )
            ax.fill_between(
                self.energy,
                0,
                -self.total_dos(squeeze=True)[1],
                lw=0,
                color="red",
                alpha=0.3,
                label="DOS (down)",
            )

            ax.plot(
                self.energy,
                self.total_pdos(squeeze=True)[0],
                "-",
                lw=1,
                color="blue",
                alpha=0.7,
                label="PDOS (up)",
            )
            ax.plot(
                self.energy,
                -self.total_pdos(squeeze=True)[1],
                "-",
                lw=1,
                color="red",
                alpha=0.7,
                label="PDOS (down)",
            )
            ncol = 2
        if self.case == 3:
            ax.fill_between(
                self.energy,
                0,
                self.total_dos(squeeze=True),
                lw=0,
                color="grey",
                alpha=0.3,
                label="DOS",
            )
            ax.fill_between(
                self.energy,
                -self.total_dos(squeeze=True),
                0,
                lw=0,
                color="grey",
                alpha=0.3,
                label="-DOS",
            )

            ax.plot(
                self.energy,
                self.total_pdos(squeeze=True)[0],
                "-",
                lw=1,
                color="blue",
                alpha=0.7,
                label="PDOS (up)",
            )
            ax.plot(
                self.energy,
                -self.total_pdos(squeeze=True)[1],
                "-",
                lw=1,
                color="red",
                alpha=0.7,
                label="PDOS (down)",
            )
            ncol = 1

        ax.legend(loc="best", ncol=ncol, draggable=True)

        if interactive:
            plt.show()
        else:
            png_path = f"{output_name}.png"
            plt.savefig(png_path, dpi=600, bbox_inches="tight")
            if save_pickle:
                import pickle

                with open(f"{png_path}.pickle", "wb") as file:
                    pickle.dump(fig, file)
        plt.close()


class DOSIterator:
    def __init__(self, dos: DOSQE) -> None:
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

r"""
DOS
"""
import re
from os import walk
from os.path import join
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
from termcolor import cprint

from radtools.decorate.axes import plot_hlines, plot_vlines
from radtools.dos.pdos import PDOSQE

PATTERN = "\\.pdos_atm\\#[0-9]*\\([a-zA-Z0-9]*\\)_wfc\\#[0-9]*\\([spdf_0-9j\\.]*\\)"


class DOSQE:
    r"""
    Analyse the folder with the |QE|_ results and store all information.

    Parameters
    ----------
    seedname : str
        Seedname for the |QE|_ output files.
    input_folder : str
        Path to the folder with |QE|_ output files.
    energy_window : tuple, optional
        Energy limits, necessary for correct plotting.

    Attributes
    ----------
    seedname : str
        Seedname for the |QE|_ output files.
    energy : array
        Energy values of the DOS/PDOS.
    k_resolved : bool
    case : int
    casename : str
        Human readable name of the case base on :py:attr:`.case`.
    nkpoints : int
        Number of k points for k resolved DOS. 1 if DOS is not k resolved.
    nepoints : int
        Number of energy points.
    filenames : list
    atoms : list
    """

    def __init__(
        self, seedname: str, input_folder: str, energy_window=None, efermi=0
    ) -> None:
        self.seedname = seedname
        self._input_folder = input_folder

        # Case and whether DOS is k resolved
        self._k_resolved = None
        self._case = None

        # Detect number of energy and k points and read energy values
        self.nkpoints = None
        self.nepoints = None
        self.energy = None
        self.efermi = efermi
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
    def spin_pol(self):
        return self.case in [2, 3]

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

        Headers of "seedname".pdos_tot according to the
        |projwfc|_:

        * Collinear:
            E DOS(E) PDOS(E)

        * Collinear, spin-polarized:
            E DOSup(E) DOSdw(E) PDOSup(E) PDOSdw(E)

        * Non-collinear, non-spin-orbit:
            E DOS(E) PDOSup(E) PDOSdw(E)

        * Non-collinear, spin-orbit:
            E DOS(E) PDOS(E)
        """

        # Detect the case
        if self._case is None:
            # Read header
            with open(join(f"{self._input_folder}", f"{self.seedname}.pdos_tot")) as file:
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
            # Differentiate between case 1 and 4
            if self._case == 1:
                # pattern is intentionally different from PATTERN
                pattern = (
                    "\\.pdos_atm\\#[0-9]*\\([a-zA-Z]*\\)_wfc\\#[0-9]*\\([spdf\\.]*\\)"
                )
                filename = self.filenames[0]
                if not re.fullmatch(
                    f"{re.escape(self.seedname)}{pattern}",
                    filename,
                ):
                    self._case = 4

            if self._case is None:
                raise RuntimeError(
                    "Unable to detect case, analysed header:\n"
                    + f"{header}\n"
                    + f"of the file:\n"
                    + f"{join(f'{self._input_folder}', f'{self.seedname}.pdos_tot')}"
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
            with open(join(f"{self._input_folder}", f"{self.seedname}.pdos_tot")) as file:
                self._k_resolved = "ik" in file.readline().lower().split()

        return self._k_resolved

    def _extract_energy(self):
        # Load data
        dos = np.loadtxt(
            join(f"{self._input_folder}", f"{self.seedname}.pdos_tot"), skiprows=1
        ).T

        # Load energy
        if self.k_resolved:
            # ik E ...
            self.nkpoints = int(dos[0][-1])
            self.nepoints = len(dos[0]) // self.nkpoints
            self.energy = dos[1][0 : self.nepoints] - self.efermi
        else:
            # E ...
            self.nkpoints = 1
            self.nepoints = len(dos[0])
            self.energy = dos[0] - self.efermi

    @property
    def filenames(self):
        r"""
        List of filenames (without "seedname".pdos_tot).
        """
        if self._filenames is None:
            # Get list of files in the folder
            all_filenames = []
            for _, _, filenames in walk(self._input_folder):
                all_filenames.extend(filenames)
                break

            self._filenames = []
            # Search
            for filename in all_filenames:
                if re.fullmatch(
                    f"{re.escape(self.seedname)}{PATTERN}",
                    filename,
                ):
                    self._filenames.append(filename)

        return self._filenames

    def _decompose_filenames(self):
        r"""
        Decompose filenames and extract information about atoms and projectors.

        Notes
        -----
        atoms : dist
            Dictionary with the atom labels and their numbers.

            .. code-block:: python

                _atoms = {atom_symbol: [atom_number, ...], ...}
        wfcs : dict
            Dictionary of projectors wave functions and their numbers.

            .. code-block:: python

                _wfcs = {(atom_symbol, atom_number): {wfc_symbol: [wfc_number, ...], ...}, ...}
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
                self._wfcs[(atom_symbol, atom_number)] = {}
            if wfc_symbol not in self._wfcs[(atom_symbol, atom_number)]:
                self._wfcs[(atom_symbol, atom_number)][wfc_symbol] = []
            self._wfcs[(atom_symbol, atom_number)][wfc_symbol].append(wfc_number)

        # Sort entries
        for atom in self._atoms:
            self._atoms[atom] = list(set(self._atoms[atom]))
            self._atoms[atom].sort()
        for key in self._wfcs:
            for wfc_symbol in self._wfcs[key]:
                self._wfcs[key][wfc_symbol].sort()

    def __iter__(self):
        return DOSIterator(self)

    def __contains__(self, item):
        atom, atom_number, wfc, wfc_number = item
        return (
            atom in self.atoms
            and atom_number in self.atom_numbers(atom)
            and wfc in self.wfcs(atom, atom_number)
            and wfc_number in self.wfc_numbers(atom, wfc, atom_number)
        )

    def total_dos(self, squeeze=False):
        r"""
        Total density of states.

        Density of states, computed directly from plane-wave basis set.

        If DOS is k resolved, then sum over k points.

        Parameters
        ----------
        squeeze : bool, default False
            Whether to sum over k points. Ignored if DOS is not k resolved.
        Whether to fix updown in the noncollinear, non spin-orbit case.

        Returns
        -------
        total_dos : :numpy:`ndarray`
            Shape :math:`(2, n_e)` if DOS is collinear, spin-polarized;
            :math:`(n_e)` otherwise.
        """

        # Load data
        dos = np.loadtxt(
            join(self._input_folder, f"{self.seedname}.pdos_tot"), skiprows=1
        ).T

        if self.case == 2:
            # ik E DOS(up) DOS(down) PDOS(up) PDOS(down)
            if self.k_resolved:
                reshaped_dos = dos[2:4].reshape((2, self.nkpoints, self.nepoints))[
                    :, :, self.energy_window[0] : self.energy_window[1]
                ]
                # Squeeze if necessary
                if squeeze:
                    return np.sum(reshaped_dos, axis=1) / self.nkpoints
                return reshaped_dos

            # E DOS(up) DOS(down) PDOS(up) PDOS(down)
            return dos[1:3][:, self.energy_window[0] : self.energy_window[1]]

        # ik E DOS PDOS
        # or
        # ik E DOS PDOS(up) PDOS(down)
        if self.k_resolved:
            reshaped_dos = dos[2].reshape((self.nkpoints, self.nepoints))[
                :, self.energy_window[0] : self.energy_window[1]
            ]
            # Squeeze if necessary
            if squeeze:
                return np.sum(reshaped_dos, axis=0) / self.nkpoints
            return reshaped_dos

        # E DOS PDOS
        # or
        # E DOS PDOS(up) PDOS(down)
        return dos[1][self.energy_window[0] : self.energy_window[1]]

    def total_pdos(self, squeeze=False):
        r"""
        Total partial density of states.

        If DOS is k resolved, then sum over k points.

        Parameters
        ----------
        squeeze : bool, default False
            Whether to sum over k points. Ignored if DOS is not k resolved.

        Returns
        -------
        total_pdos : :numpy:`ndarray`
            Shape :math:`(2, n_e)` if PDOS is collinear, spin-polarized;
            :math:`(n_e)` otherwise.
        """

        # Load data
        dos = np.loadtxt(
            join(f"{self._input_folder}", f"{self.seedname}.pdos_tot"), skiprows=1
        ).T

        if self.case == 2:
            # ik E DOS(up) DOS(down) PDOS(up) PDOS(down)
            if self.k_resolved:
                reshaped_pdos = dos[4:6].reshape((2, self.nkpoints, self.nepoints))[
                    :, :, self.energy_window[0] : self.energy_window[1]
                ]
                # Squeeze if necessary
                if squeeze:
                    return np.sum(reshaped_pdos, axis=1) / self.nkpoints
                return reshaped_pdos

            # E DOS(up) DOS(down) PDOS(up) PDOS(down)
            return dos[3:5][:, self.energy_window[0] : self.energy_window[1]]

        elif self.case == 3:
            # ik E DOS PDOS(up) PDOS(down)
            if self.k_resolved:
                reshaped_pdos = dos[3:5].reshape((2, self.nkpoints, self.nepoints))[
                    :, :, self.energy_window[0] : self.energy_window[1]
                ]
                # Squeeze if necessary
                if squeeze:
                    total_pdos = np.sum(reshaped_pdos, axis=1) / self.nkpoints
                return reshaped_pdos

            # E DOS PDOS(up) PDOS(down)
            return dos[2:4][:, self.energy_window[0] : self.energy_window[1]]

        # ik E DOS PDOS
        if self.k_resolved:
            reshaped_pdos = dos[3].reshape((self.nkpoints, self.nepoints))[
                :, self.energy_window[0] : self.energy_window[1]
            ]
            # Squeeze if necessary
            if squeeze:
                return np.sum(reshaped_pdos, axis=0) / self.nkpoints
            return reshaped_pdos

        # E DOS PDOS
        return dos[2][self.energy_window[0] : self.energy_window[1]]

    @property
    def atoms(self):
        r"""
        List of atom's types.

        Returns
        -------
        atoms : list
            List of atom's types symbols

            .. code-block:: python

                atoms = [atom_symbol, ...]
        """

        return [atom for atom in self._atoms]

    def atom_numbers(self, atom: str):
        r"""
        List of atom's numbers for particular atom type.

        Parameters
        ----------
        atom : str
            Label of atom.

        Returns
        -------
        numbers : list

            .. code-block:: python

                numbers = [atom_number, ...]
        """

        return self._atoms[atom]

    def wfcs(self, atom: str, atom_number: int = None):
        r"""
        Return list of wave function symbols for particular atom.

        Parameters
        ----------
        atom : str
            Label of an atom.
        atom_number : int, optional
            Number of an atom. If ``None`` then return wfc symbols of first atom.

        Returns
        -------
        wfc_symbols : list
            List of wave function symbols

            .. code-block:: python

                wfc_symbols = [wfc_symbol, ...]

        """

        # If atom_number is not provided, then return wfc symbols of first atom
        if atom_number is None:
            atom_number = self.atom_numbers(atom)[0]
        return [symbol for symbol in self._wfcs[(atom, atom_number)]]

    def wfc_numbers(self, atom: str, wfc_symbol, atom_number: int = None):
        r"""
        Return list of wave function numbers for particular atom and wave function type.

        Parameters
        ----------
        atom : str
            Label of an atom.
        wfc_symbol : str
            Label of a wave function type.
        atom_number : int, optional
            Number of an atom. If ``None`` then return wfc numbers of first atom.

        Returns
        -------
        wfc_numbers : list
            List of wave function numbers

            .. code-block:: python

                wfc_numbers = [wfc_number, ...]

        """

        # If atom_number is not provided, then return wfc numbers of first atom
        if atom_number is None:
            atom_number = self.atom_numbers(atom)[0]
        return self._wfcs[(atom, atom_number)][wfc_symbol]

    def pdos(
        self, atom, wfc, wfc_numbers=None, atom_numbers=None, background_total=False
    ) -> PDOSQE:
        r"""
        Projected density of states for a particular atom.

        Parameters
        ----------
        atom : str
            Name of the atom type.
        wfc : str
            Name of the projector wave function.
        wfc_numbers : Iterable or int, optional
            Number of wave function projector.
        atom_numbers : Iterable or int, optional
            If ``None``, then PDOS summed over all atom numbers for ``atom``.

        Returns
        -------
        pdos : :py:class:`PDOSQE`
            Partial density of states, orbital-resolved.
        """

        # Check if atom is valid
        if atom not in self.atoms:
            raise ValueError(f"There is no {atom} in PDOS.")

        # Prepare atom_numbers
        if atom_numbers is None:
            atom_numbers = self.atom_numbers(atom)
        elif isinstance(atom_numbers, int):
            atom_numbers = [atom_numbers]
        elif not isinstance(atom_numbers, Iterable):
            raise ValueError(
                f"atom_numbers must be None, int or Iterable, "
                + f"got {type(atom_numbers)} :\n"
                + f"{atom_numbers}"
            )

        # Check if wfc is valid
        if wfc not in self.wfcs(atom, atom_numbers[0]):
            raise ValueError(f"There is no {wfc} projector of {atom} atom in PDOS.")

        # Load and sum data for all atom numbers of the atom and all wfc numbers
        for i, atom_number in enumerate(atom_numbers):
            # Check if atom number is valid
            if atom_number not in self.atom_numbers(atom):
                raise ValueError(f"There is no {atom}#{atom_number} atom in PDOS.")

            # Prepare wfc_numbers
            if wfc_numbers is None:
                wfc_numbers = self.wfc_numbers(atom, wfc, atom_number)
            elif isinstance(wfc_numbers, int):
                wfc_numbers = [wfc_numbers]
            elif not isinstance(wfc_numbers, Iterable):
                raise ValueError(
                    f"wfc_numbers must be None, int or Iterable,"
                    + f" got {type(wfc_numbers)} :\n"
                    + f"{wfc_numbers}"
                )

            for j, wfc_number in enumerate(wfc_numbers):
                # Check if wfc number is valid
                if wfc_number not in self.wfc_numbers(atom, wfc, atom_number):
                    raise ValueError(
                        f"There is no {wfc}#{wfc_number} wave-function for "
                        + f"{atom}#{atom_number} atom in PDOS."
                    )
                # Load and sum data
                path = join(
                    self._input_folder,
                    f"{self.seedname}.pdos_atm#{atom_number}({atom})_wfc#{wfc_number}({wfc})",
                )
                if i == 0 and j == 0:
                    pdos = np.loadtxt(path, skiprows=1).T
                else:
                    pdos += np.loadtxt(path, skiprows=1).T

        # Reshape data
        if self.spin_pol:
            if self.k_resolved:
                # ik E LDOS(up) LDOS(down) PDOS_1(up) PDOS_1(down) ... PDOS_2l+1(up) PDOS_2l+1(down)
                ldos = pdos[2:4].reshape(2, self.nkpoints, self.nepoints)[
                    :, :, self.energy_window[0] : self.energy_window[1]
                ]
                pdos = pdos[4:].reshape(
                    (pdos.shape[0] - 4) // 2, 2, self.nkpoints, self.nepoints
                )[:, :, :, self.energy_window[0] : self.energy_window[1]]
            else:
                # E LDOS(up) LDOS(down) PDOS_1(up) PDOS_1(down) ... PDOS_2l+1(up) PDOS_2l+1(down)
                ldos = pdos[1:3][:, self.energy_window[0] : self.energy_window[1]]
                pdos = pdos[3:].reshape((pdos.shape[0] - 3) // 2, 2, self.nepoints)[
                    :, :, self.energy_window[0] : self.energy_window[1]
                ]
        else:
            if self.k_resolved:
                # ik E LDOS PDOS_1 ... PDOS_2l+1
                ldos = pdos[2].reshape(self.nkpoints, self.nepoints)[
                    :, self.energy_window[0] : self.energy_window[1]
                ]
                pdos = pdos[3:].reshape(
                    pdos.shape[0] - 3, self.nkpoints, self.nepoints
                )[:, :, self.energy_window[0] : self.energy_window[1]]
            else:
                # E LDOS PDOS_1 ... PDOS_2l+1
                ldos = pdos[1][self.energy_window[0] : self.energy_window[1]]
                pdos = pdos[2:][:, self.energy_window[0] : self.energy_window[1]]

        if background_total:
            ldos = self.total_pdos()
        return PDOSQE(
            energy=self.energy,
            pdos=pdos,
            projectors_group=wfc,
            ldos=ldos,
            spin_pol=self.spin_pol,
        )

    def plot_pdos_tot(
        self,
        output_name,
        interactive=False,
        efermi=0,
        xlim=None,
        ylim=None,
        save_pickle=False,
        axes_labels_fontsize=18,
        legend_fontsize=12,
        title_fontsize=18,
    ):
        r"""
        Plot total DOS vs total PDOS.

        Parameters
        ----------
        output_name : str
            Name of the output file.
        interactive : bool, default False
            Whether to plot in interactive |matplotlib|_ window.
        efermi : float, default 0
            Fermi energy. Zero of energy scale is shifted to it.
        xlim : tuple, optional
            Limits for energy scale.
        ylim : tuple, optional
            Limits for dos scale.
        save_pickle : bool, default False
            Whether to save figure as a .pickle file.
            Helps for custom modification of particular figures.
        axes_label_fontsize : int, default 18
            Fontsize of the axes labels.
        legend_fontsize : int, default 12
            Fontsize of the legend.
        title_fontsize : int, default 18
            Title fontsize.
        """
        fig, ax = plt.subplots(figsize=(8, 4))

        #  x axis label
        if efermi == 0:
            ax.set_xlabel("E, eV", fontsize=axes_labels_fontsize)
        else:
            ax.set_xlabel("E - E$_{Fermi}$, eV", fontsize=axes_labels_fontsize)
        # y axis label
        ax.set_ylabel("DOS, states/eV", fontsize=axes_labels_fontsize)

        # x axis limits
        if xlim is None:
            xlim = (np.amin(self.energy), np.amax(self.energy))
        ax.set_xlim(*tuple(xlim))
        # y axis limits
        if ylim is not None:
            ax.set_ylim(*tuple(ylim))

        # Title
        ax.set_title(f"DOS vs PDOS", fontsize=title_fontsize)

        # Eye-guide lines
        if efermi != 0:
            plot_vlines(ax, 0)
        if self.spin_pol:
            plot_hlines(ax, 0)

        def fill(data, color, label):
            ax.fill_between(
                self.energy, 0, data, lw=0, color=color, alpha=0.3, label=label
            )

        def plot(data, color, label):
            ax.plot(self.energy, data, "-", lw=1, color=color, alpha=0.7, label=label)

        # E DOS(E) PDOS(E)
        if self.case in [1, 4]:
            fill(self.total_dos(squeeze=True), "grey", "DOS")
            plot(self.total_pdos(squeeze=True), "black", "PDOS")
            ncol = 1
        # E DOSup(E) DOSdw(E) PDOSup(E) PDOSdw(E)
        if self.case == 2:
            fill(self.total_dos(squeeze=True)[0], "blue", "DOS (up)")
            fill(-self.total_dos(squeeze=True)[1], "red", "DOS (down)")
            plot(self.total_pdos(squeeze=True)[0], "blue", "PDOS (up)")
            plot(-self.total_pdos(squeeze=True)[1], "red", "PDOS (down)")
            ncol = 2
        # E DOS(E) PDOSup(E) PDOSdw(E)
        if self.case == 3:
            fill(self.total_dos(squeeze=True), "grey", "DOS")
            fill(-self.total_dos(squeeze=True), "grey", "$-$DOS")
            plot(self.total_pdos(squeeze=True)[0], "blue", "PDOS (up)")
            plot(-self.total_pdos(squeeze=True)[1], "red", "PDOS (down)")
            ncol = 1

        if interactive:
            ax.legend(loc="best", ncol=ncol, draggable=True, fontsize=legend_fontsize)
        else:
            ax.legend(loc="best", ncol=ncol, fontsize=legend_fontsize)

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
        self._list = []
        for atom in dos.atoms:
            for atom_number in dos.atom_numbers(atom):
                for wfc in dos.wfcs(atom, atom_number):
                    for wfc_number in dos.wfc_numbers(atom, wfc, atom_number):
                        self._list.append([atom, atom_number, wfc, wfc_number])
        self._index = 0

    def __next__(self) -> str:
        if self._index < len(self._list):
            result = self._list[self._index]
            self._index += 1
            return result
        raise StopIteration

    def __iter__(self):
        return self


def detect_seednames(input_folder):
    r"""
    Analyze input folder, detects seednames for the dos output files.

    Parameters
    ----------
    input_folder : str
        Directory with DOS files.

    Returns
    -------
    seednames : list
        List of seednames found in ``input_folder``.
    """

    # Get list of files in the folder
    files = []
    for _, _, filenames in walk(input_folder):
        files.extend(filenames)
        break

    seednames = set()
    for file in files:
        if ".pdos_tot" in file and ".pdos_tot" == file[-9:]:
            seednames.add(file[:-9])
        elif re.match(f".*{PATTERN}$", file):
            seednames.add(re.split(f"{PATTERN}$", file)[0])
    seednames = list(seednames)

    return seednames


def prepare_custom_pdos(dos: DOSQE, custom, quiet=True):
    r"""
    Prepare custom PDOS. Based on the input line.

    Parameters
    ----------
    dos : :py:class:`.DOSQE`
        DOS input files wrapper.
    custom : list of str
        List of strings describing the custom PDOS.
        See :ref:`rad-plot-dos_custom-plots` for the string`s description.
    quiet : bool, default True
        Whether to print information about the runtime.

    Returns
    -------
    pdos : list of :numpy:`ndarray`
        PDOS array. Not the instance of the :py:class:`.PDOS`, but
        an array, which could be passed to the constructor of
        :py:class:`.PDOS` class
    """
    pdos = []

    # Process each custom specification
    for entry in custom:
        # Print info about current specification
        if not quiet:
            cprint(f'"{entry}":', "green")

        # Initialize pdos_element corresponding to the entry
        pdos_element = None

        # Remove spaces and closing brackets
        entry = entry.replace(" ", "").replace(")", "")

        # Split by semicolon into different atom types
        subentries = entry.split(";")

        # Process each atom specification
        for subentry in subentries:
            # Get the info about atoms
            atom_part = subentry.split("(")[0]

            # Get the info about atom type
            atom = atom_part.split("#")[0]

            # Get the info about atom numbers
            if "#" in atom_part:
                # Atom numbers provided directly
                atom_numbers = list(map(int, atom_part.split("#")[1:]))
                # Check if the custom numbers are correct
                if not set(atom_numbers).issubset(set(dos.atom_numbers(atom))):
                    raise ValueError(
                        f"For the atom {atom} the following numbers exist:\n"
                        + f"{dos.atom_numbers(atom)}\n"
                        + f"but the following numbers are provided:\n"
                        + f"{atom_numbers}"
                    )
            else:
                # Atom numbers are taken from the dos object
                atom_numbers = dos.atom_numbers(atom)

            # Print the summary of the atom with numbers
            if not quiet:
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

            # Get the info about wave function projectors, if any
            if "(" in subentry:
                # Get all wave function types
                wfc_parts = subentry.split("(")[1].split(",")
                wfcs = {}
                # Process each wave function type
                for wfc_part in wfc_parts:
                    # Get the name of the wave function type
                    wfc = wfc_part.split("#")[0]

                    # If wave function numbers are provided directly
                    if "#" in wfc_part:
                        wfc_numbers = list(map(int, wfc_part.split("#")[1:]))
                        if wfc in wfcs:
                            wfcs[wfc].extend(wfc_numbers)
                        else:
                            wfcs[wfc] = wfc_numbers
                    # If wave function numbers are not provided
                    else:
                        # It will result in the summ over all numbers
                        wfcs[wfc] = dos.wfc_numbers(atom, wfc, atom_numbers[0])

            # If no wave function projectors are provided, then sum over all projectors
            else:
                wfcs = dict(
                    [
                        (wfc, dos.wfc_numbers(atom, wfc, atom_numbers[0]))
                        for wfc in dos.wfcs(atom)
                    ]
                )

            if not quiet:
                print(
                    f"  * For each {atom} atom PDOS is summed among the following projections:",
                    end="\n      ",
                )
                for wfc in wfcs:
                    for number in wfcs[wfc]:
                        print(f"{wfc}#{number}", end=" ")
                print()

            for wfc in wfcs:
                tmp_pdos = dos.pdos(
                    atom=atom,
                    wfc=wfc,
                    wfc_numbers=wfcs[wfc],
                    atom_numbers=atom_numbers,
                ).ldos
                if pdos_element is None:
                    pdos_element = tmp_pdos
                else:
                    pdos_element += tmp_pdos

        pdos.append(pdos_element)

    return np.array(pdos)

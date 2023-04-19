#! /usr/local/bin/python3

import re
from argparse import ArgumentParser
from copy import deepcopy
from os import makedirs, walk
from os.path import abspath, join

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import __version__ as matplotlib_version

from rad_tools.routines import cprint


class PDOS:
    r"""
    Class which deals with the individual files of projected density of states.

    Parameters
    ----------
    case : int
        Case of the calculations.
    projector_family : str
        Type of projections
        (s, p, d, f + value of j in the case of spin-orbit).

    Attributes
    ----------
    case : int
    projector_family : str
        Type of projections
        (s, p, d, f + value of j in the case of spin-orbit).
    dos : array
        Projected density of states.
        Formats with according to the
        :projwfc:`QE <>`:

        * Collinear:
            LDOS(E) PDOS_1(E) ... PDOS_2l+1(E)

        * Collinear, spin-polarized:
            LDOS_up(E) LDOS_down(E) PDOS_1_up(E) PDOS_1_down(E) ... PDOS_2l+1_up(E) PDOS_2l+1_down(E)

        * Non-collinear, non-spin-orbit:
            LDOS_up(E) LDOS_down(E) PDOS_1_up(E) PDOS_1_down(E) ... PDOS_2l+1_up(E) PDOS_2l+1_down(E)

        * Non-collinear, spin-orbit:
            LDOS(E) PDOS_1(E) ... PDOS_2j+1(E)
    projectors : list
        List of projectors.
    ldos : array
    """

    _projectors = {
        "s": ["s"],
        "p": ["$p_z$", "$p_y$", "$p_x$"],
        "d": ["$d_{z^2}$", "$d_{zx}$", "$d_{zy}$", "$d_{x^2 - y^2}$", "$d_{xy}$"],
        "f": [
            "$f_{z^3}$",
            "$f_{yz^2}$",
            "$f_{xz^2}$",
            "$f_{z(x^2 - y^2)}$",
            "$f_{xyz}$",
            "$f_{y(3x^2 - y^2)}$",
            "$f_{x(x^2 - 3y^2)}$",
        ],
    }

    def __init__(self, case: int, projector_family: str, dos, projectors=None):
        self.case = case
        self.projector_family = projector_family
        self.dos = dos
        if self.projector_family in self._projectors:
            self.projectors = self._projectors[self.projector_family]
        elif re.fullmatch("[spdf]_[0-9]*\.[0-9]*"):  # TODO check re expression
            self.projectors = [
                f"{self.projector_family.split('_')[0]} ($m_j = {i+1}$)"
                for i in range(
                    0, int(2 * float(self.projector_family.split("_")[1]) + 1)
                )
            ]
        else:
            self.projectors = projectors

    def __add__(self, other):
        if isinstance(other, PDOS):
            dos = self.dos + other.dos
            return PDOS(case=self.case, projector_family=self.projector_family, dos=dos)
        else:
            raise TypeError

    def __iter__(self):
        return PDOSIterator(self)

    def __contains__(self, item):
        return item in self.projectors

    def __getitem__(self, key) -> np.ndarray:
        if self.case in [1, 4]:
            return self.dos[1 + self.projectors.index(key)]
        else:
            return self.dos[
                1 + 2 * self.projectors.index(key) : 3 + 2 * self.projectors.index(key)
            ]

    @property
    def ldos(self):
        r"""
        Local density of states.

        Returns
        -------
        ldos : array
        """
        if self.case in [1, 4]:
            return self.dos[0]
        else:
            return self.dos[0:2]


class PDOSIterator:
    def __init__(self, pdos: PDOS) -> None:
        self._projectors = pdos.projectors
        self._index = 0

    def __next__(self) -> str:
        if self._index < len(self._projectors):
            result = self._projectors[self._index]
            self._index += 1
            return result
        raise StopIteration

    def __iter__(self):
        return self


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

    @property
    def total_dos(self):
        r"""
        Total density of states.

        Density of states, computed directly from plane-wave basis set.
        """
        dos = np.loadtxt(
            join(self._input_path, f"{self.seedname}.pdos_tot"), skiprows=1
        ).T

        if self.case == 2:
            if self.k_resolved:
                return dos[2:4].reshape((2, self.nkpoints, self.nepoints))[
                    :, :, self.window[0] : self.window[1]
                ]
            else:
                return dos[1:3][:, self.window[0] : self.window[1]]
        if self.k_resolved:
            return dos[2].reshape((self.nkpoints, self.nepoints))[
                :, self.window[0] : self.window[1]
            ]
        else:
            return dos[1][self.window[0] : self.window[1]]

    @property
    def total_pdos(self):
        r"""
        Total density of states.

        Density of states, computed directly from plane-wave basis set.
        """
        dos = np.loadtxt(f"{self._input_path}/{self.seedname}.pdos_tot", skiprows=1).T

        if self.case == 2:
            if self.k_resolved:
                return dos[4:6].reshape((2, self.nkpoints, self.nepoints))[
                    :, :, self.window[0] : self.window[1]
                ]
            else:
                return dos[3:5][:, self.window[0] : self.window[1]]
        elif self.case == 3:
            if self.k_resolved:
                return dos[3:5].reshape((2, self.nkpoints, self.nepoints))[
                    :, :, self.window[0] : self.window[1]
                ]
            else:
                return dos[2:4][:, self.window[0] : self.window[1]]
        if self.k_resolved:
            return dos[3].reshape((self.nkpoints, self.nepoints))[
                :, self.window[0] : self.window[1]
            ]
        else:
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

    def pdos(self, atom, atom_number, wfc, wfc_number):
        r"""
        Projected density of states for a particular atom.

        If ``atom_number`` = 0 summ all pdos over the atoms of the same type.
        """

        if atom_number == 0:
            for i, atom_number in enumerate(self.atom_numbers(atom)):
                path = join(
                    self._input_path,
                    f"{self.seedname}.pdos_atm#{atom_number}({atom})_wfc#{wfc_number}({wfc})",
                )
                tmp_dos = np.loadtxt(path, skiprows=1).T[1:]
                if self.k_resolved:
                    tmp_dos = tmp_dos.reshape(
                        (tmp_dos.shape[0], self.nkpoints, self.nepoints)
                    )
                if i == 0:
                    dos = tmp_dos
                else:
                    dos += tmp_dos
            dos /= len(self.atom_numbers(atom))
        else:
            path = join(
                self._input_path,
                f"{self.seedname}.pdos_atm#{atom_number}({atom})_wfc#{wfc_number}({wfc})",
            )
            dos = np.loadtxt(path, skiprows=1).T[1:]
            if self.k_resolved:
                dos = dos.reshape((dos.shape[0], self.nkpoints, self.nepoints))
        if self.k_resolved:
            dos = dos[:, :, self.window[0] : self.window[1]]
        else:
            dos = dos[:, self.window[0] : self.window[1]]
        return PDOS(case=self.case, projector_family=wfc, dos=dos)

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
            total_dos = self.total_dos
            total_pdos = self.total_pdos
            if self.k_resolved:
                total_dos = np.sum(total_dos, axis=0)
                total_pdos = np.sum(total_pdos, axis=0)

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
            total_dos = self.total_dos
            total_pdos = self.total_pdos
            if self.k_resolved:
                total_dos = np.sum(total_dos, axis=1) / self.nkpoints
                total_pdos = np.sum(total_pdos, axis=1) / self.nkpoints

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


def plot_projected(
    energy,
    dos,
    labels,
    main_label,
    output_name,
    title,
    ylim=None,
    relative=False,
    normalize=False,
    updown=False,
):
    colours = [
        "#0000FF",
        "#FF0000",
        "#00FF00",
        "#FF00FF",
        "#00FFFF",
        "#3E3847",
        "#FFD600",
        "#366B35",
        "#FF6F00",
    ]
    nsubplots = len(labels)
    if relative:
        fig, ax = plt.subplots(figsize=(9, 4))
    else:
        fig, axs = plt.subplots(nsubplots, 1, figsize=(9, nsubplots * 2))
        if nsubplots == 1:
            axs = [axs]
    fig.subplots_adjust(hspace=0)

    def set_up_axis(ax, i):
        if normalize:
            ax.set_ylabel("PDOS / LDOS", fontsize=12)
        else:
            ax.set_ylabel("DOS, states/eV", fontsize=12)
        if i == nsubplots - 1:
            ax.set_xlabel("E, ev", fontsize=15)
        else:
            ax.axes.get_xaxis().set_visible(False)
        if ylim is not None:
            ax.set_ylim(*tuple(ylim))
        ax.set_xlim(np.amin(energy), np.amax(energy))
        ax.vlines(
            0,
            0,
            1,
            transform=ax.get_xaxis_transform(),
            color="grey",
            linewidths=0.5,
            linestyles="dashed",
        )
        if title is not None and i == 0:
            ax.set_title(title)

    if relative:
        set_up_axis(ax, nsubplots - 1)
        ax.hlines(
            0,
            0,
            1,
            transform=ax.get_yaxis_transform(),
            color="black",
            linewidths=1,
        )
        if normalize:
            if updown:
                dos[2] = np.where(dos[0] > 10e-8, dos[2] / dos[0], 0)
                dos[3] = np.where(dos[1] > 10e-8, dos[3] / dos[1], 0)
            else:
                dos[1] = np.where(dos[0] > 10e-8, dos[1] / dos[0], 0)

        for i in range(1, nsubplots):
            if not updown:
                if normalize:
                    dos[1 + i] = np.where(dos[0] > 10e-8, dos[1 + i] / dos[0], 0)
                dos[1 + i] += dos[i]
            else:
                if normalize:
                    dos[2 + 2 * i] = np.where(
                        dos[0] > 10e-8, dos[2 + 2 * i] / dos[0], 0
                    )
                    dos[3 + 2 * i] = np.where(
                        dos[1] > 10e-8, dos[3 + 2 * i] / dos[1], 0
                    )
                dos[2 + 2 * i] += dos[1 + 2 * (i - 1)]
                dos[3 + 2 * i] += dos[2 + 2 * (i - 1)]

    for i in range(0, nsubplots):
        if not relative:
            ax = axs[i]
            set_up_axis(ax, i)

        if not updown:
            if relative:
                colour = colours[i // len(colours)]
                ax.fill_between(
                    energy,
                    dos[i],
                    dos[1 + i],
                    lw=0,
                    color=colour,
                    alpha=0.3,
                    label=main_label,
                )
            else:
                colour = "black"
                ax.fill_between(
                    energy, 0, dos[0], lw=0, color=colour, alpha=0.3, label=main_label
                )
                ax.plot(
                    energy,
                    dos[1 + i],
                    "-",
                    lw=0.5,
                    color=colour,
                    alpha=0.8,
                    label=labels[i],
                )
        else:
            if relative:
                colour_up = colours[i // len(colours)]
                colour_down = colours[i // len(colours)]
                ax.fill_between(
                    energy,
                    dos[2 + 2 * (i - 1)],
                    dos[2 + 2 * i],
                    lw=0,
                    color=colour_up,
                    alpha=0.3,
                    label=f"{main_label} (up)",
                )
                ax.fill_between(
                    energy,
                    -dos[2 + 2 * (i - 1)],
                    -dos[2 + 2 * i],
                    lw=0,
                    color=colour_down,
                    alpha=0.3,
                    label=f"{main_label} (down)",
                )
            else:
                colour_up = "blue"
                colour_down = "red"
                ax.fill_between(
                    energy,
                    0,
                    dos[1],
                    lw=0,
                    color=colour_up,
                    alpha=0.3,
                    label=f"{main_label} (up)",
                )
                ax.fill_between(
                    energy,
                    0,
                    -dos[2],
                    lw=0,
                    color=colour_down,
                    alpha=0.3,
                    label=f"{main_label} (down)",
                )

                ax.plot(
                    energy,
                    dos[2 + 2 * i],
                    "-",
                    lw=0.5,
                    color=colour_up,
                    alpha=0.8,
                    label=f"{labels[i]} (up)",
                )
                ax.plot(
                    energy,
                    -dos[2 + 2 * i],
                    "-",
                    lw=0.5,
                    color=colour_down,
                    alpha=0.8,
                    label=f"{labels[i]} (down)",
                )

        ax.legend(loc=(1.025, 0.2), bbox_transform=ax.transAxes)

    plt.savefig(f"{output_name}.png", dpi=400, bbox_inches="tight")
    plt.close()


def detect_seednames(input_path):
    r"""
    Analyze input folder, detects seednames for the dos output files.

    Parameters
    ----------
    input_path : str
        Directory with DOS files.

    Returns
    -------
    seednames : list
        List of seednames found in ``input_path``.
    """

    pattern = ".pdos_atm#[0-9]*\\([a-zA-Z]*\\)_wfc#[0-9]*\\([spdf_0-9j.]*\\)"

    # Get list of files in the folder
    files = []
    for dirpath, dirnames, filenames in walk(input_path):
        files.extend(filenames)
        break

    seednames = set()
    for file in files:
        if ".pdos_tot" in file and ".pdos_tot" == file[-9:]:
            seednames.add(file[:-9])
        elif re.match(f".*{pattern}$", file):
            seednames.add(re.split(f"{pattern}$", file)[0])
    seednames = list(seednames)

    return seednames


def manager(
    input_path,
    seedname=None,
    output_path=".",
    energy_window=None,
    dos_window=None,
    efermi=0.0,
    separate=False,
    relative=False,
    normalize=False,
    verbose=False,
    interactive=False,
):
    makedirs(output_path, exist_ok=True)
    suffix = ""
    if relative:
        suffix += "r"
    if normalize:
        suffix += "n"

    if seedname is None:
        seednames = detect_seednames(input_path)
        print(f"Following DOS seednames (filpdos) are detected:")
        for item in seednames:
            cprint(f"   * {item}", colour="yellow")
    else:
        seednames = [seedname]

    for s_i, seedname in enumerate(seednames):
        output_root = join(output_path, f"{seedname}-{suffix}")
        # Preparations
        cprint(
            f"({s_i + 1}/{len(seednames)}) Start to work with {seedname} seedname",
            colour="yellow",
        )
        makedirs(output_root, exist_ok=True)

        dos = DOS(seedname, input_path, efermi=efermi, window=energy_window)
        print(f"{dos.case_name} case detected.")

        # Plot PDOS vs DOS
        dos.plot_pdos_tot(
            output_name=join(output_root, "pdos-vs-dos"),
            interactive=interactive,
            ylim=dos_window,
        )

        # Plot PDOS for each atom/wfc
        for atom in dos.atoms:
            for atom_number in dos.atom_numbers(atom):
                for wfc, wfc_number in dos.wfcs(atom, atom_number):
                    plot_projected(
                        dos.energy, dos.pdos(atom, atom_number, wfc, wfc_number)
                    )


def create_parser():
    parser = ArgumentParser(
        description="Script for visualisation of density of states."
    )
    parser.add_argument(
        "-ip",
        "--input-path",
        metavar="path",
        type=str,
        required=True,
        help="Relative or absulute path to the folder with dos files.",
    )
    parser.add_argument(
        "-s",
        "--seedname",
        metavar="filpdos",
        type=str,
        default=None,
        help="Prefix for output files containing PDOS(E). "
        + "As specified in the QE projwfc.x input file.",
    )
    parser.add_argument(
        "-op",
        "--output-path",
        metavar="path",
        type=str,
        default=".",
        help="Relative or absolute path to the folder for saving outputs.",
    )
    parser.add_argument(
        "-ew",
        "--energy-window",
        metavar=("min", "max"),
        type=float,
        nargs=2,
        default=None,
        help="Energy window for the plots. "
        + "By default whole range present in the files is plotted.",
    )
    parser.add_argument(
        "-dw",
        "--dos-window",
        metavar=("min", "max"),
        type=float,
        nargs=2,
        default=None,
        help="DOS window for the plots. "
        + "By default whole range present in the files is plotted.",
    )
    parser.add_argument(
        "-ef",
        "--efermi",
        metavar="energy",
        type=float,
        default=0.0,
        help="Fermi energy, zero by default.",
    )
    parser.add_argument(
        "-sep",
        "--separate",
        action="store_true",
        default=False,
        help="Whenever to plot projected DOS for each atom of the same type separately.",
    )
    parser.add_argument(
        "-r",
        "--relative",
        action="store_true",
        default=False,
        help="Whenever to use relative style.",
    )
    parser.add_argument(
        "-n",
        "--normalize",
        action="store_true",
        default=False,
        help="Whenever to use normalize relative style.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Verbose output, propagates to the called methods.",
    )
    parser.add_argument(
        "-i",
        "--interactive",
        action="store_true",
        default=False,
        help="Interactive plotting.",
    )

    return parser

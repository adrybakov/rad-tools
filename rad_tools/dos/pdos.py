r"""
PDOS
"""

import re
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np


class PDOS:
    r"""
    Partial density of states, projected on arbitrary projections.

    Supports k-resolved density of states.
    Support spin-polarised and spin-unpolarised cases.

    PDOS class is iterable (over ``projectors``) and
    supports item call (return PDOS for ``key`` projector).

    Operations of addition and substruction are defined.

    Parameters
    ----------
    energy : array
        Values of energy for the PDOS. Shape :math:`n_e` is assumed.
    pdos : array
        Array with the values of PDOS. Shape is assumed to be:

        * Spin-polarized, k-resolved: :math:`(n, 2, n_k, n_e)`
        * Spin-unpolarized, k-resolved: :math:`(n, n_k, n_e)`
        * Spin-polarized, non k-resolved: :math:`(n, 2, n_e)`
        * Spin-unpolarized, non k-resolved: :math:`(n, n_e)`

        where :math:`n` is the number of projections,
        :math:`n_k` is the number of k-points,
        :math:`n_e` is the number of energy points.
    projectors_group : str
        Name of the projectors group.
    projectors : list
        Names of the projectors.
        If ``projectors_group`` has the form "l" or "l_j",
        where l is "s", "p", "d", "f" and j is the total angular momentum,
        the projectors are assigned automatically,
        otherwise it is necessary to provide :math:`n` projectors manually.
        The names of projectors are directly used in the plots.
    ldos : array, default None
        Local density of states. Sum of partial density of states over all projectors.
        Computed based on ``pdos`` if not provided.
        Shape is assumed to be:

        * Spin-polarized, k-resolved: :math:`(2, n_k, n_e)`
        * Spin-unpolarized, k-resolved: :math:`(n_k, n_e)`
        * Spin-polarized, non k-resolved: :math:`(2, n_e)`
        * Spin-unpolarized, non k-resolved: :math:`(n_e)`

        where :math:`n_k` is the number of k-points,
        :math:`n_e` is the number of energy points.
    spin_pol : bool, default False
        Whenever PDOS is spin-polarized or not.

    Attributes
    ----------
    energy : array
        Values of energy for the PDOS. Has the shape :math:`n_e`.
    ldos : array
    pdos : array
    projectors_group : str
        Name of the projectors group.
    projectors : list
        Names of the projectors.
    spin_pol : bool, default False
        Whenever PDOS is spin-polarized or not.
    k_resolved : bool, default False
    """

    def __init__(
        self,
        energy,
        pdos,
        projectors_group: str,
        projectors,
        ldos=None,
        spin_pol=False,
    ):
        self.energy = energy
        self._pdos = None
        self._ldos = None
        self.projectors_group = projectors_group
        self.projectors = projectors
        self.spin_pol = spin_pol

        self.pdos = pdos
        self.ldos = ldos

    def __add__(self, other):
        if not isinstance(other, PDOS):
            raise TypeError(
                f"Addition is not supported between "
                + f"{type(self)} and {type(other)}"
            )
        if (
            self._pdos.shape == other._pdos.shape
            and self.spin_pol == other.spin_pol
            and self.projectors_group == other.projectors_group
            and set(self.projectors) == set(other.projectors)
        ):
            pdos = self._pdos + other._pdos
            ldos = self._ldos + other._ldos

            return PDOS(
                energy=self.energy,
                pdos=pdos,
                projectors_group=self.projectors_group,
                projectors=self.projectors,
                ldos=ldos,
                spin_pol=self.spin_pol,
            )
        else:
            raise ValueError(
                "There is a mismatch between self and other:\n"
                + f"    pdos shape: {self._pdos.shape} {other._pdos.shape}\n"
                + f"    spin_pol: {self.spin_pol} {other.spin_pol}\n"
                + f"    projectors_group: {self.projectors_group} {other.projectors_group}\n"
                + f"    projectors: {self.projectors} {other.projectors}\n"
            )

    def __sub__(self, other):
        if not isinstance(other, PDOS):
            raise TypeError(
                f"Subtraction is not supported between "
                + f"{type(self)} and {type(other)}"
            )
        if (
            self._pdos.shape == other._pdos.shape
            and self.spin_pol == other.spin_pol
            and self.projectors_group == other.projectors_group
            and set(self.projectors) == set(other.projectors)
        ):
            pdos = self._pdos - other._pdos
            ldos = self._ldos - other._ldos
            return PDOS(
                energy=self.energy,
                pdos=pdos,
                projectors_group=self.projectors_group,
                projectors=self.projectors,
                ldos=ldos,
                spin_pol=self.spin_pol,
            )
        else:
            raise ValueError(
                "There is a mismatch between self and other:\n"
                + f"    pdos shape: {self._pdos.shape} {other._pdos.shape}\n"
                + f"    spin_pol: {self.spin_pol} {other.spin_pol}\n"
                + f"    projectors_group: {self.projectors_group} {other.projectors_group}\n"
                + f"    projectors: {self.projectors} {other.projectors}\n"
            )

    def __iter__(self):
        return PDOSIterator(self)

    def __contains__(self, item):
        return item in self.projectors

    def __getitem__(self, key) -> np.ndarray:
        r"""
        Return pdos by the projector name.

        Parameters
        ----------
        key : str or int
            Projector`s name or index.
        Returns
        -------
        pdos : array
            Partial density of states.
        """
        if isinstance(key, str):
            key = self.projectors.index(key)

        return self.pdos[key]

    @property
    def ldos(self) -> np.ndarray:
        r"""
        Local density of states.

        Returns
        -------
        ldos : array
            Summed density of states along all projectors.
            Has the following shapes:

            * Spin-polarized, k-resolved: :math:`(2, n_k, n_e)`
            * Spin-unpolarized, k-resolved: :math:`(n_k, n_e)`
            * Spin-polarized, non k-resolved: :math:`(2, n_e)`
            * Spin-unpolarized, non k-resolved: :math:`(n_e)`

            where :math:`n_k` is the number of k-points,
            :math:`n_e` is the number of energy points.
        """

        return self._ldos

    @ldos.setter
    def ldos(self, new_ldos):
        if (
            self._pdos is not None
            and new_ldos is not None
            and np.array(new_ldos).shape != self.pdos.shape[1:]
        ):
            raise ValueError(
                f"New LDOS shape {new_ldos.shape} does not match "
                + f"PDOS shape ({self.pdos.shape[1:]})"
            )

        if new_ldos is not None:
            self._ldos = np.array(new_ldos)
        else:
            self._ldos = np.sum(self._pdos, axis=0)

        if len(self.energy) != self.ldos.shape[-1]:
            raise ValueError(
                f"LDOS does not match with energy: "
                + f"{len(self.energy)} energy points, {self.ldos.shape[-1]} LDOS points"
            )

    @property
    def pdos(self) -> np.ndarray:
        r"""
        Partial density of states.

        Returns
        -------
        pdos : array
            Has the following shapes:

            * Spin-polarized, k-resolved: :math:`(n, 2, n_k, n_e)`
            * Spin-unpolarized, k-resolved: :math:`(n, n_k, n_e)`
            * Spin-polarized, non k-resolved: :math:`(n, 2, n_e)`
            * Spin-unpolarized, non k-resolved: :math:`(n, n_e)`

            where :math:`n` is the number of projections,
            :math:`n_k` is the number of k-points,
            :math:`n_e` is the number of energy points.
        """

        return self._pdos

    @pdos.setter
    def pdos(self, new_pdos):
        new_pdos = np.array(new_pdos)
        if self._ldos is not None and new_pdos.shape[1:] != self.ldos.shape:
            raise ValueError(
                f"New PDOS shape {new_pdos.shape[1:]} does not match "
                + f"LDOS shape ({self.ldos.shape})"
            )
        self._pdos = new_pdos

        if len(self.energy) != self.pdos.shape[-1]:
            raise ValueError(
                f"PDOS does not match with energy: "
                + f"{len(self.energy)} energy points, {self.pdos.shape[-1]} PDOS points"
            )

        if len(self.projectors) != self.pdos.shape[0]:
            raise ValueError(
                f"PDOS does not match with projectors: "
                + f"{len(self.projectors)} projectors, {self.pdos.shape[0]} PDOS"
            )

    @property
    def k_resolved(self):
        r"""
        Check if pdos is k-resolved based on shape of ``self.pdos``.
        """

        return (
            self.spin_pol
            and len(self.pdos.shape) == 4
            or not self.spin_pol
            and len(self.pdos.shape) == 3
        )

    def squeeze(self):
        r"""
        Squeeze k-resolved PDOS.

        Note
        ----
        It modifies the instance on which called.
        See also :py:func:`PDOS.squeezed`
        """

        if self.k_resolved:
            self._pdos = np.sum(self._pdos, axis=1 + int(self.spin_pol))
            self._ldos = np.sum(self._ldos, axis=int(self.spin_pol))

    def squeezed(self):
        r"""
        Return new instance with squeezed PDOS.

        Calls :py:func:`PDOS.squeeze`.

        Returns
        -------
        pdos_squeezed : :py:class:`PDOS`
            Squeezed PDOS.
        """

        squeezed_pdos = deepcopy(self)
        squeezed_pdos.squeeze()
        return squeezed_pdos

    def normalize(self):
        r"""
        Normalize values of PDOS to 1 for each k and energy point.

        If :math:`x_i(E, k)` is  PDOS of the projector :math:`i`,
        then after this function does the following:

        .. math::
            x_i(E, k) \rightarrow \dfrac{x_i(E, k)}{\sum_{j=0}^{n} x_j(E, k)}

        where :math:`n` is the total number of projectors.
        Those sums are computed individually for spin-up and
        spin-down in the spin-polarized case.

        Note
        ----
        It modifies the instance on which called.
        See also :py:func:`PDOS.normalized`
        """

        for i in range(0, self.pdos.shape[0]):
            self._pdos[i] = np.where(self.ldos > 10e-8, self._pdos[i] / self.ldos, None)
        self._ldos = np.where(self.ldos > 10e-8, self._ldos / self.ldos, 0)

    def normalized(self):
        r"""
        Return new instance with normalized PDOS.

        Calls :py:func:`PDOS.normalize`.

        Returns
        -------
        normalized_pdos : :py:class:`PDOS`
            Normalized PDOS.
        """

        normalized_pdos = deepcopy(self)
        normalized_pdos.normalize()
        return normalized_pdos

    def dump_txt(self, ouptut_name):
        r"""
        Save PDOS as .txt file.

        First line is a header.
        """

        header = ""
        fmt = ""
        if self.k_resolved:
            header += "# ik  "
            fmt += "%5.0f "
        header += "E (eV)     "
        fmt += "%10.3f "

        if self.spin_pol:
            header += "ldosup(E) ldosdw(E) "
            fmt += "%9.3E %9.3E "
            for i in self.projectors:
                n = max(9, len(i) + 2)
                header += f"{i+'up':>{n}} {i+'dw':>{n}} "
                fmt += f"%{n}.3E %{n}.3E "
        else:
            header += "ldos(E)   "
            fmt += "%9.3E "
            for i in self.projectors:
                n = max(9, len(i))
                header += f"{i:>{n}} "
                fmt += f"%{n}.3E "

        data = []
        nepoints = self.pdos.shape[-1]
        if self.k_resolved:
            if self.spin_pol:
                nkpoints = self.pdos.shape[2]
                ldos = self.ldos.reshape(2, nkpoints * nepoints)
                pdos = self.pdos.reshape(len(self.projectors), 2, nkpoints * nepoints)
            else:
                nkpoints = self.pdos.shape[1]
                ldos = self.ldos.reshape(nkpoints * nepoints)
                pdos = self.pdos.reshape(len(self.projectors), nkpoints * nepoints)
            data.append(np.repeat(np.linspace(1, nkpoints, nkpoints), nepoints))
            data.append(np.tile(self.energy, nkpoints))
        else:
            data.append(self.energy)
            ldos = self.ldos
            pdos = self.pdos

        if self.spin_pol:
            data.append(ldos[0])
            data.append(ldos[1])
            for i in range(pdos.shape[0]):
                data.append(pdos[i][0])
                data.append(pdos[i][1])
        else:
            data.append(ldos)
            for i in range(pdos.shape[0]):
                data.append(pdos[i])
        np.savetxt(ouptut_name, np.array(data).T, fmt=fmt, header=header, comments="")


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


class PDOSQE(PDOS):
    r"""
    PDOS wrapper for Quantum Espresso pdos.

    Supports the order of projectors of :projwfc:`Quantum Espresso <>` (s,p,d,f) and
    the case of projection in the spin-orbit calculations.
    In the custom cases it is necessary to specify projectors manually.
    If ``projectors_group`` has the form "l" or "l_j",
    where l is "s", "p", "d", "f" and j is the total angular momentum,
    the projectors are assigned automatically,
    otherwise it is necessary to provide :math:`n` projectors manually.
    The names of projectors are directly used in the plots.
    If ``projectors_group`` is one of "s", "p", "d", "f", then the projectors are:

    * s : :math:`s`
    * p : :math:`p_z`, :math:`p_y`, :math:`p_x`
    * d : :math:`d_{z^2}`, :math:`d_{zx}`, :math:`d_{zy}`, :math:`d_{x^2 - y^2}`, :math:`d_{xy}`
    * f : :math:`f_{z^3}`, :math:`f_{yz^2}`, :math:`f_{xz^2}`, :math:`f_{z(x^2 - y^2)}`, :math:`f_{xyz}`, :math:`f_{y(3x^2 - y^2)}`, :math:`f_{x(x^2 - 3y^2)}`

    If ``projectors_group`` has the form "l_j", then the projectors are :math:`(1, ..., 2j+1)`

    """

    _pattern = "[spdf]_j[0-9.]*"
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

    def __init__(
        self,
        energy,
        pdos,
        projectors_group: str,
        projectors=None,
        ldos=None,
        spin_pol=False,
    ):
        if projectors is not None:
            pass
        elif projectors_group in self._projectors:
            projectors = self._projectors[projectors_group]
        elif re.fullmatch(self._pattern, projectors_group):
            l, j = projectors_group.split("_j")
            m_j = range(1, 1 + int(2 * float(j) + 1))
            projectors = [f"{l} ($m_j = {i}$)" for i in m_j]
        else:
            raise ValueError(
                "Projectors can not be assigned automatically, "
                + "you have to provide explicit list of projectors. "
                + f"Projectors group: {projectors_group}"
            )
        super().__init__(energy, pdos, projectors_group, projectors, ldos, spin_pol)


def plot_projected(
    pdos: PDOS,
    efermi=0.0,
    output_name="pdos",
    title=None,
    xlim=None,
    ylim=None,
    relative=False,
    normalize=False,
    interactive=False,
    save_pickle=False,
):
    r"""
    Plot PDOS.

    Parameters
    ----------
    pdos : PDOS
        PDOS for the plot.
    efermi : float, default 0
        Fermi energy.
    output_name : str, default "pdos"
        output_name for the plot file. Extension ".png" is added at the end.
    title : str, default None
        Title of the plot. Passed to the ``ax.set_title()``.
    xlim : tuple
        limits for the x (Energy) axis
    ylim : tuple
        limits for the y (PDOS) axis
    relative : bool, default False
        Relative plot style.
    normalize : bool, default False
        Whenever to norma;ize relative plot style.
    interactive : bool, default False
        Whenever to use interactive plotting mode.
    save_pickle : bool, default False
        Whether to save figure as a .pickle file.
        Helps for custom modification of particular figures.
    """

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
    n = len(pdos.projectors)
    pdos = pdos.squeezed()

    if relative:
        fig, ax = plt.subplots(figsize=(9, 4))
    else:
        fig, axs = plt.subplots(n, 1, figsize=(9, n * 2))
        if n == 1:
            axs = [axs]
    fig.subplots_adjust(hspace=0)

    def set_up_axis(ax, i):
        if normalize:
            ax.set_ylabel("PDOS / LDOS", fontsize=12)
        else:
            ax.set_ylabel("DOS, states/eV", fontsize=12)
        if i == n - 1:
            ax.set_xlabel("E, ev", fontsize=15)
        else:
            ax.axes.get_xaxis().set_visible(False)
        if ylim is not None:
            ax.set_ylim(*tuple(ylim))
        if xlim is not None:
            ax.set_xlim(*tuple(xlim))
        else:
            ax.set_xlim(np.amin(pdos.energy), np.amax(pdos.energy))
        ax.vlines(
            0,
            0,
            1,
            transform=ax.get_xaxis_transform(),
            color="grey",
            linewidths=0.5,
            linestyles="dashed",
        )
        if title is not None and (i == 0 or relative):
            ax.set_title(title)

    if normalize:
        pdos = pdos.normalized()

    if relative:
        set_up_axis(ax, n - 1)
        ax.hlines(
            0,
            0,
            1,
            transform=ax.get_yaxis_transform(),
            color="black",
            linewidths=1,
        )

    for i, projector in enumerate(pdos):
        if not relative:
            ax = axs[i]
            set_up_axis(ax, i)
        if pdos.spin_pol:
            if relative:
                ax.fill_between(
                    pdos.energy - efermi,
                    np.sum(pdos[:i], axis=0)[0],
                    np.sum(pdos[: i + 1], axis=0)[0],
                    lw=0,
                    color=colours[i % len(colours)],
                    alpha=0.3,
                    label=f"{projector} (up)",
                )
                ax.fill_between(
                    pdos.energy - efermi,
                    np.sum(pdos[:i], axis=0)[1],
                    -np.sum(pdos[: i + 1], axis=0)[1],
                    lw=0,
                    color=colours[i % len(colours)],
                    alpha=0.3,
                    label=f"{projector} (down)",
                )
            else:
                ax.fill_between(
                    pdos.energy - efermi,
                    0,
                    pdos.ldos[0],
                    lw=0,
                    color="blue",
                    alpha=0.3,
                    label=f"{pdos.projectors_group} (up)",
                )
                ax.fill_between(
                    pdos.energy - efermi,
                    0,
                    -pdos.ldos[1],
                    lw=0,
                    color="red",
                    alpha=0.3,
                    label=f"{pdos.projectors_group} (down)",
                )

                ax.plot(
                    pdos.energy - efermi,
                    pdos[projector][0],
                    "-",
                    lw=0.5,
                    color="blue",
                    alpha=0.8,
                    label=f"{projector} (up)",
                )
                ax.plot(
                    pdos.energy - efermi,
                    -pdos[projector][1],
                    "-",
                    lw=0.5,
                    color="red",
                    alpha=0.8,
                    label=f"{projector} (down)",
                )
        else:
            if relative:
                ax.fill_between(
                    pdos.energy - efermi,
                    np.sum(pdos[:i], axis=0),
                    np.sum(pdos[: i + 1], axis=0),
                    lw=0,
                    color=colours[i % len(colours)],
                    alpha=0.3,
                    label=projector,
                )
            else:
                ax.fill_between(
                    pdos.energy - efermi,
                    0,
                    pdos.ldos,
                    lw=0,
                    color="black",
                    alpha=0.3,
                    label=pdos.projectors_group,
                )
                ax.plot(
                    pdos.energy - efermi,
                    pdos[projector],
                    "-",
                    lw=0.5,
                    color="black",
                    alpha=0.8,
                    label=projector,
                )

        ax.legend(loc=(1.025, 0.2), bbox_transform=ax.transAxes)

    if interactive:
        plt.show()
    else:
        plt.savefig(f"{output_name}.png", dpi=400, bbox_inches="tight")
        if save_pickle:
            import pickle

            with open(f"{png_path}.pickle", "wb") as file:
                pickle.dump(fig, file)
    plt.close()

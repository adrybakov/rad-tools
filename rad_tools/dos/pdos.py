r"""
PDOS
"""

import re
import numpy as np
import matplotlib.pyplot as plt


class PDOS:
    r"""
    Partial density of states, projected on arbitrary projections.

    Supports the order of projectors of :projwfc:`Quantum Espresso <>` (s,p,d,f) and
    the case of projection in the spin-orbit calculations.
    In the custom cases it is necessary to specify projectors manually.

    Supports k-resolved density of states.
    Support spin-polarised and spin-unpolarised cases.

    PDOS class is iterable (over ``projectors``) and
    supports item call (return PDOS for ``key`` projector).

    Operations of addition and substruction are defined.

    Parameters
    ----------
    energy : array
        Values of energy for the PDOS. Sshape :math:`n_e` is assumed.
    pdos : array
        Array with the values of PDOS.
        For the k-resolved case shape :math:`(m, n_k, n_e)` is assumed.
        For the non k-resolved case  shape :math:`(m, n_e)` is assumed.
        If PDOS has :math:`n` projectors, then first index `math:`m`
        is assumed to indicate the folowing data:

        * Spin-polarized case (:math:`m = 2n + 2`):

          LDOS_up LDOS_down PDOS_1_up PDOS_1_down ... PDOS_n_up PDOS_n_down

        * Spin-unpolarized case (:math:`m = n + 1`):

          LDOS PDOS_1 ... PDOS_n

    projectors_group : str
        Name of the projectors group.
    projectors : list
        Names of the projectors.
        If ``projectors_group`` has the form "l" or "l_j",
        where l is "s", "p", "d", "f" and j is the total angular momentum,
        the projectors are assigned automatically,
        otherwise it is necessary to provide :math:`n` projectors manually.
        The names of projectors are directly used in the plots.
    spin_pol : bool, default False
        Whenever PDOS is spin-polarized or not.
    k_resolved : bool, default False
        Whenever PDOS is k_resolved.

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

    spin_pol : bool, default False
        Whenever PDOS is spin-polarized or not.
    k_resolved : bool, default False
        Whenever PDOS is k_resolved.
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

    _pattern = "[spdf]_[0-9.]*"

    def __init__(
        self,
        energy,
        pdos,
        projectors_group: str,
        projectors=None,
        spin_pol=False,
        k_resolved=False,
    ):
        self.energy = energy
        self._pdos = np.array(pdos)

        self.projectors_group = projectors_group
        if self.projectors_group in self._projectors:
            self.projectors = self._projectors[self.projectors_group]
        elif re.fullmatch(self._pattern, self.projectors_group):
            l, j = self.projectors_group.split("_")
            m_j = range(1, 1 + int(2 * float(j) + 1))
            self.projectors = [f"{l} ($m_j = {i}$)" for i in m_j]
        elif projectors is None:
            raise ValueError(
                "Projectors can not be assigned automatically, "
                + "you have to provide explicit list of projectors. "
                + f"Projectors group: {self.projectors_group}."
            )
        else:
            self.projectors = projectors

        self.spin_pol = spin_pol
        self.k_resolved = k_resolved

    def __add__(self, other):
        if not isinstance(other, PDOS):
            raise TypeError(
                f"Addition is not supported between " + "{type(self)} and {type(other)}"
            )
        if (
            self.k_resolved == other.k_resolved
            and self.spin_pol == other.spin_pol
            and self.projectors_group == other.projectors_group
            and set(self.projectors) == set(other.projectors)
        ):
            pdos = self._pdos + other._pdos
            return PDOS(
                energy=self.energy,
                pdos=pdos,
                projectors_group=self.projectors_group,
                projectors=self.projectors,
                spin_pol=self.spin_pol,
                k_resolved=self.k_resolved,
            )
        else:
            raise ValueError(
                "There is a mismatch between self and other:\n"
                + f"    k_resolved: {self.k_resolved} {other.k_resolved}\n"
                + f"    spin_pol: {self.spin_pol} {other.spin_pol}\n"
                + f"    projectors_group: {self.projectors_group} {other.projectors_group}\n"
                + f"    projectors: {self.projectors} {other.projectors}\n"
            )

    def __sub__(self, other):
        if not isinstance(other, PDOS):
            raise TypeError(
                f"Subtraction is not supported between "
                + "{type(self)} and {type(other)}"
            )
        if (
            self.k_resolved == other.k_resolved
            and self.spin_pol == other.spin_pol
            and self.projectors_group == other.projectors_group
            and set(self.projectors) == set(other.projectors)
        ):
            pdos = self._pdos - other._pdos
            return PDOS(
                energy=self.energy,
                pdos=pdos,
                projectors_group=self.projectors_group,
                projectors=self.projectors,
                spin_pol=self.spin_pol,
                k_resolved=self.k_resolved,
            )
        else:
            raise ValueError(
                "There is a mismatch between self and other:\n"
                + f"    k_resolved: {self.k_resolved} {other.k_resolved}\n"
                + f"    spin_pol: {self.spin_pol} {other.spin_pol}\n"
                + f"    projectors_group: {self.projectors_group} {other.projectors_group}\n"
                + f"    projectors: {self.projectors} {other.projectors}\n"
            )

    def __iter__(self):
        return PDOSIterator(self)

    def __contains__(self, item):
        return item in self.projectors

    def __getitem__(self, key) -> np.ndarray:
        if key == "zeros_placeholder":
            return np.zeros(self.pdos[0].shape, dtype=float)
        return self.pdos[self.projectors.index(key)]

    @property
    def ldos(self):
        r"""
        Local density of states.

        Parameters
        ----------
        squeeze : bool, default False
            Affects only k-resolved case.
            Sum LDOS among k-points and divide by the number of k-points.

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

        if self.spin_pol:
            return self._pdos[0:2]
        else:
            return self._pdos[0]

    @property
    def pdos(self, squeeze=False):
        r"""
        Partial density of states.

        Returns
        -------
        pdos : array
            Summed density of states along all projectors.
            Has the following shapes:

            * Spin-polarized, k-resolved: :math:`(n, 2, n_k, n_e)`
            * Spin-unpolarized, k-resolved: :math:`(n, n_k, n_e)`
            * Spin-polarized, non k-resolved: :math:`(n, 2, n_e)`
            * Spin-unpolarized, non k-resolved: :math:`(n, n_e)`

            where :math:`n` is the number of projections,
            :math:`n_k` is the number of k-points,
            :math:`n_e` is the number of energy points.
        """

        if self.spin_pol:
            shape = self._pdos.shape
            return self._pdos[2:].reshape((shape[0] // 2, 2, shape[1:]))
        else:
            return self._pdos[1:]

    def squeeze(self):
        r"""
        Squeeze k-resolved PDOS.
        """

        if self.k_resolved:
            self._pdos = np.sum(self._pdos, axis=1)

    def make_relative(self, normalize=False):
        r"""
        Recomputes to the relative values for the custom plots.

        If :math:`x_i(E)` is  PDOS of the projector :math:`i`,
        then after this function does the following:

        .. math::
            x_i(E) \rightarrow \sum_{j=0}^{i} x_j(E)

        With ``normalize`` = ``True`` it does the following:

        .. math::
            x_i(E) \rightarrow \dfrac{\sum_{j=0}^{i} x_j(E)}{\sum_{j=0}^{n} x_j(E)}

        where :math:`n` is the total number of projectors.
        Those sums are computed individually for spin-up and spin-down in the spin-polarized case.

        Parameters
        ----------
        normalize : bool, default False
            Whenever to normalize PDOS to 1 for each energy value.
        """
        if self.spin_pol:
            for i in range(2, len(self._pdos) // 2):
                self._pdos[2 + 2 * i] += self._pdos[1 + 2 * (i - 1)]
                self._pdos[3 + 2 * i] += self._pdos[2 + 2 * (i - 1)]
        else:
            for i in range(2, len(self._pdos)):
                self._pdos[i] += self._pdos[i - 1]

        if normalize:
            if self.spin_pol:
                for i in range(2, len(self._pdos) // 2):
                    self._pdos[2 + 2 * i] = np.where(
                        self.ldos[0] > 10e-8, self._pdos[2 + 2 * i] / self.ldos[0], 0
                    )
                    self._pdos[3 + 2 * i] = np.where(
                        self.ldos[1] > 10e-8, self._pdos[3 + 2 * i] / self.ldos[1], 0
                    )
            else:
                for i in range(2, len(self._pdos)):
                    self._pdos[i] += self._pdos[i - 1]
                    self._pdos[i] = np.where(
                        self.ldos > 10e-8, self._pdos[i] / self.ldos, 0
                    )


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
    if pdos.k_resolved:
        pdos.squeeze()

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

    if relative:
        pdos.make_relative(normalize=normalize)
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
        if i != 0:
            prev_projector = pdos.projectors[i - 1]
        else:
            prev_projector = "zeros_placeholder"
        if not relative:
            ax = axs[i]
            set_up_axis(ax, i)
        if pdos.spin_pol:
            if relative:
                ax.fill_between(
                    pdos.energy - efermi,
                    pdos[prev_projector][0],
                    pdos[projector][0],
                    lw=0,
                    color=colours[i // len(colours)],
                    alpha=0.3,
                    label=f"{projector} (up)",
                )
                ax.fill_between(
                    pdos.energy - efermi,
                    pdos[prev_projector][1],
                    pdos[projector][1],
                    lw=0,
                    color=colours[i // len(colours)],
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
                    pdos[prev_projector],
                    pdos[projector],
                    lw=0,
                    color=colours[i // len(colours)],
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
    plt.close()


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

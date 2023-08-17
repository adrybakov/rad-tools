import os

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

ROOT_DIR = "."


def plot_notation():
    colors = {
        f"{1}": to_rgb("#33CB86"),
        f"{2}": to_rgb("#50D42E"),
        f"{0.5}": to_rgb("#D5DD2B"),
        f"{1}" + R"$\cdot\vert S_i\vert \vert S_j\vert$": to_rgb("#E66027"),
        f"{2}" + R"$\cdot\vert S_i\vert \vert S_j\vert$": to_rgb("#EE2480"),
        f"{0.5}" + R"$\cdot\vert S_i\vert \vert S_j\vert$": to_rgb("#CB22F6"),
        f"{4}": to_rgb("#2A24FE"),
        f"{4}" + R"$\cdot\vert S_i\vert \vert S_j\vert$": to_rgb("#29C0FF"),
    }
    quiver_style = dict(
        angles="xy",
        scale_units="xy",
        scale=1,
        width=0.003,
        headlength=7,
        headwidth=4,
        headaxislength=6,
        color="grey",
    )
    text_style = dict(va="center", ha="center", fontsize=14)
    delta = 5
    ylim = (0, 31)
    xlim = (0.7 * delta, 6.9 * delta)
    fig, ax = plt.subplots(figsize=(xlim[1] / 2, ylim[1] / 2))
    ax.text(
        delta,
        31 / 2,
        "???",
        va="center",
        ha="right",
        fontsize=20,
        bbox=dict(
            boxstyle="round",
            ec=(0.1, 0.1, 0.1, 1),
            fc=(0.1, 0.1, 0.1, 0.5),
        ),
    )
    for dcount in range(2):
        dcount_pos = 31 / 4 * (1 + 2 * dcount)
        message = "Double\ncounting" if dcount else "No\ndouble\ncounting"
        ax.text(2 * delta, dcount_pos, message, **text_style)
        ax.quiver(delta + 0.5, 31 / 2, delta - 1.7, dcount_pos - 31 / 2, **quiver_style)
        for snorm in range(2):
            snorm_pos = 31 / 8 * (1 + 4 * dcount + 2 * snorm)
            message = "Spin\nnormalized" if snorm else "Spin\nnot\nnormalized"
            ax.text(3 * delta, snorm_pos, message, **text_style)
            ax.quiver(
                2 * delta + 1.2,
                dcount_pos,
                delta - 2.7,
                dcount_pos - snorm_pos,
                **quiver_style,
            )
            for onehalf in range(2):
                onehalf_pos = 31 / 16 * (1 + 8 * dcount + 4 * snorm + 2 * onehalf)
                message = (
                    "Factor $\dfrac{1}{2}$" if onehalf else "No\nfactor $\dfrac{1}{2}$"
                )
                ax.text(4 * delta, onehalf_pos, message, **text_style)
                ax.quiver(
                    3 * delta + 1.5,
                    snorm_pos,
                    delta - 2.7,
                    snorm_pos - onehalf_pos,
                    **quiver_style,
                )
                for two in range(2):
                    two_pos = (
                        31 / 32 * (1 + 16 * dcount + 8 * snorm + 4 * onehalf + 2 * two)
                    )
                    message = "Factor 2" if two else "No\nfactor 2"
                    ax.text(5 * delta, two_pos, message, **text_style)
                    ax.quiver(
                        4 * delta + 1.2,
                        onehalf_pos,
                        delta - 2.3,
                        onehalf_pos - two_pos,
                        **quiver_style,
                    )

                    ax.quiver(5 * delta + 1.2, two_pos, 2, 0, **quiver_style)
                    hamiltonian = ""
                    energy = 1.0
                    if two and not onehalf:
                        hamiltonian += "2"
                        energy *= 2
                    if onehalf and not two:
                        hamiltonian += "$\dfrac{1}{2}$"
                        energy *= 0.5
                    if dcount:
                        hamiltonian += "$\sum_{ij}$"
                        energy *= 2
                    else:
                        hamiltonian += "$\sum_{i<j}$"
                    if int(energy) != energy:
                        energy = f"{energy}"
                    else:
                        energy = f"{energy:.0f}"
                    if snorm:
                        hamiltonian += "$e_i\ J_{ij}\ e_j$"
                        conv_factor = f"{energy}"
                    else:
                        hamiltonian += "$S_i\ J_{ij}\ S_j$"
                        conv_factor = (
                            f"{energy}" + R"$\cdot\vert S_i\vert \vert S_j\vert$"
                        )
                    hamiltonian += f" = {conv_factor}"
                    ax.text(
                        5 * delta + 3.5,
                        two_pos,
                        hamiltonian,
                        va="center",
                        ha="left",
                        fontsize=15,
                        bbox=dict(
                            boxstyle="round",
                            ec=(*colors[conv_factor], 1),
                            fc=(*colors[conv_factor], 0.5),
                        ),
                    )

    ax.set_ylim(*ylim)
    ax.set_xlim(*xlim)
    ax.set_aspect(1)
    ax.axis("off")
    filename = os.path.join(ROOT_DIR, "docs", "source", "img", "notation-tree.png")
    plt.savefig(filename, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    plot_notation()

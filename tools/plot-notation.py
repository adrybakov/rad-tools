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
    delta = 6
    yrange = 25
    ylim = (0, yrange)
    xlim = (0.7 * delta, 5.7 * delta)
    factor_message = ["Factor $\dfrac{1}{2}$", "Factor 2", "No Factor"]
    fig, ax = plt.subplots(figsize=(xlim[1] / 2, ylim[1] / 2))
    ax.text(
        delta,
        yrange / 2,
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
        dcount_pos = yrange / 4 * (1 + 2 * dcount)
        message = "Double\ncounting" if dcount else "No\ndouble\ncounting"
        ax.text(2 * delta, dcount_pos, message, **text_style)
        ax.quiver(
            delta + 0.5,
            yrange / 2,
            delta - 1.7,
            dcount_pos - yrange / 2,
            **quiver_style,
        )
        for snorm in range(2):
            snorm_pos = yrange / 8 * (1 + 4 * dcount + 2 * snorm)
            message = "Spin\nnormalized" if snorm else "Spin\nnot\nnormalized"
            ax.text(3 * delta, snorm_pos, message, **text_style)
            ax.quiver(
                2 * delta + 1.2,
                dcount_pos,
                delta - 2.7,
                dcount_pos - snorm_pos,
                **quiver_style,
            )
            for factor in range(3):
                factor_pos = yrange / 24 * (1 + 12 * dcount + 6 * snorm + 2 * factor)

                message = (
                    "Factor $\dfrac{1}{2}$"
                    if factor
                    else 1
                    if 1
                    else "No\nfactor $\dfrac{1}{2}$"
                )
                ax.text(4 * delta, factor_pos, factor_message[factor], **text_style)
                ax.quiver(
                    3 * delta + 1.5,
                    snorm_pos,
                    delta - 2.7,
                    snorm_pos - factor_pos,
                    **quiver_style,
                )
                hamiltonian = ""
                energy = 1.0
                if factor == 0:
                    hamiltonian += "$\dfrac{1}{2}$"
                    energy *= 0.5
                if factor == 1:
                    hamiltonian += "2"
                    energy *= 2
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
                    conv_factor = f"{energy}" + R"$\cdot\vert S_i\vert \vert S_j\vert$"
                hamiltonian += f" = {conv_factor}"
                ax.text(
                    4 * delta + 3.5,
                    factor_pos,
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
                ax.quiver(
                    4 * delta + 1.4,
                    factor_pos,
                    1.7,
                    0,
                    **quiver_style,
                )

    ax.set_ylim(*ylim)
    ax.set_xlim(*xlim)
    ax.set_aspect(1)
    ax.axis("off")
    filename = os.path.join(ROOT_DIR, "docs", "source", "img", "notation-tree.png")
    plt.savefig(filename, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    plot_notation()

import os

import matplotlib.pyplot as plt

ROOT_DIR = "."


def class_structure():
    fig, ax = plt.subplots(figsize=(10, 5))
    bbox = dict(
        boxstyle="round",
        ec=(0, 129 / 255, 129 / 255, 1),
        fc=(0, 129 / 255, 129 / 255, 0.5),
    )
    text_style = dict(
        ha="center",
        va="center",
        bbox=bbox,
        size=15,
        backgroundcolor="#008181",
    )
    arrow_style = dict(
        angles="xy", scale_units="xy", scale=1, headlength=3, headaxislength=2.7
    )

    def legend(x, y, text, color):
        ax.quiver(x, y, 20, 0, color=color, **arrow_style)
        ax.text(x + 21, y, text, ha="left", va="center", fontsize=15)

    parent_color = "tab:green"
    build_with_color = "tab:blue"
    build_from_color = "purple"
    ax.text(15, 50, "Lattice", **text_style)
    ax.text(15, 80, "Kpoints", **text_style)
    ax.text(60, 50, "Crystal", **text_style)
    ax.text(60, 80, "Atom", **text_style)
    ax.text(114, 50, "SpinHamiltonian", **text_style)
    ax.text(114, 20, "ExchangeParameter", **text_style)
    ax.text(140, 80, "MagnonDispersion", **text_style)

    ax.quiver(26, 50, 23, 0, color=parent_color, **arrow_style)
    ax.quiver(71, 50, 21, 0, color=parent_color, **arrow_style)
    ax.quiver(114, 43, 0, -16, color=build_with_color, **arrow_style)
    ax.quiver(60, 57, 0, 16, color=build_with_color, **arrow_style)
    ax.quiver(15, 73, 0, -16, color=build_from_color, **arrow_style)
    ax.quiver(140, 73, -26, -16, color=build_from_color, **arrow_style)
    ax.set_xlim(0, 170)
    ax.set_ylim(-10, 87)
    ax.axis("off")
    legend(5, -5, "is a parent of", parent_color)
    legend(5, 5, "is build with", build_with_color)
    legend(5, 15, "is build from", build_from_color)
    filename = os.path.join(ROOT_DIR, "docs", "source", "img", "class-relation.png")
    plt.savefig(filename, dpi=600, bbox_inches="tight")
    print(f"File is save in {os.path.abspath(filename)}")


if __name__ == "__main__":
    class_structure()

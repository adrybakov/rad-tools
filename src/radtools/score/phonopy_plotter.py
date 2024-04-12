# RAD-tools - Sandbox (mainly condense matter plotting).
# Copyright (C) 2022-2024  Andrey Rybakov
#
# e-mail: anry@uv.es, web: rad-tools.org
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import yaml


def plot(file, kp=None, constant=1, on="output"):
    data_file = open(file, "r")
    data = yaml.load(data_file, Loader=yaml.FullLoader)
    data_file.close()
    X = []
    Y = []
    K_points = []
    miny = None
    maxy = None
    for i, point in enumerate(data["phonon"]):
        for j, band in enumerate(data["phonon"][i]["band"]):
            data["phonon"][i]["band"][j]["frequency"] *= constant
    for i, point in enumerate(data["phonon"]):
        print(i)
        X.append(point["distance"])
        Y.append(list(map(lambda x: x["frequency"], point["band"])))
        if (
            i == 0
            or i == len(data["phonon"]) - 1
            or data["phonon"][i]["q-position"] == data["phonon"][i + 1]["q-position"]
        ):
            K_points.append(point["distance"])
        if not miny:
            miny = min(Y[-1])
        if not maxy:
            maxy = max(Y[-1])
        miny = min(min(Y[-1]), miny)
        maxy = max(max(Y[-1]), maxy)

    print(K_points)

    fig, ax = plt.subplots()
    X = np.array(X)
    Y = np.array(Y)
    ax.plot(X, Y, color="red")
    ax.set_ylabel("Phonon frequency (THz)", fontsize=15)
    ax.set_xlim(X[0], X[-1])
    if kp is None or len(kp) != len(K_points):
        if kp is None:
            kp = []
        print(
            f"\033[93m"
            f"{len(K_points)} k-points expected, only {len(kp)} k-points are provided. "
            f"Please, check the input."
            f"\nProvided K-points: {kp}"
            f"\x1b[0m"
        )
    if kp and len(kp) == len(K_points):
        for i, k in enumerate(kp):
            if k == "G":
                kp[i] = "$\Gamma$"
        ax.xaxis.set_ticks(K_points, kp, fontsize=15)
        for k_point in K_points:
            ax.vlines(
                k_point,
                0,
                1,
                transform=ax.get_xaxis_transform(),
                color="black",
                linewidths=0.5,
            )
    plt.show()
    if constant != 1:
        output = yaml.dump(data)
        with open(f"{on}.yaml", "w", encoding="utf-8") as file:
            file.write(output)


def main():
    parser = ArgumentParser(description="Fit Exchange parameters map")

    parser.add_argument(
        "-file", type=str, default="band.yaml", help="Path to the band.yaml file"
    )
    parser.add_argument(
        "-kp", type=str, default=None, nargs="*", help="K-path, example: -kp G X R Y G"
    )
    parser.add_argument(
        "-constant",
        type=float,
        default=1,
        help="Constant for conversion. " "The values will be multiplied by it",
    )
    parser.add_argument("-on", type=str, default="output", help="output name file")
    args = parser.parse_args()
    if args.file:
        plot(args.file, kp=args.kp, constant=args.constant, on=args.on)


if __name__ == "__main__":
    main()

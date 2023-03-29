#! /usr/local/bin/python3

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
        with open(f"{on}.yaml", "w") as file:
            file.write(output)


if __name__ == "__main__":
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

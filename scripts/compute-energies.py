#! /usr/local/bin/python3

from argparse import ArgumentParser, RawTextHelpFormatter
from subprocess import run


def J_values(FM, AFM1, AFM2, AFM3):
    J2 = -(FM - AFM1) * 13605 / 8 / 8
    J3 = -(FM - AFM2) * 13605 / 2 / 2 / 8 - J2
    J1 = -(AFM2 - AFM3) * 13605 / 8 / 2 / 8 + J3

    return J1, J2, J3


def compute(filename):
    FM = {}
    AFM1 = {}
    AFM2 = {}
    AFM3 = {}
    file = open(filename, "r")
    for line in file:
        mode = int(line.split("mode")[0])
        conf = line.split("/")[1]
        negpos = line.split("crsbr.scf")[1].split(".out")[0]
        energy = float(line.split("=")[1].split()[0])
        if "neg" in negpos:
            negpos = -int(negpos.split("neg")[0])
        elif "pos" in negpos:
            negpos = int(negpos.split("pos")[0])

        if mode not in FM:
            FM[mode] = {}
        if mode not in AFM1:
            AFM1[mode] = {}
        if mode not in AFM2:
            AFM2[mode] = {}
        if mode not in AFM3:
            AFM3[mode] = {}

        if conf == "FM":
            FM[mode][negpos] = energy
        elif conf == "AF1":
            AFM1[mode][negpos] = energy
        elif conf == "AF2":
            AFM2[mode][negpos] = energy
        elif conf == "AF3":
            AFM3[mode][negpos] = energy

    modes_list = [key for key in FM]
    modes_list.sort()
    print("MODES:\n")
    print(f"    {modes_list}\n\n")
    print("DISPLACEMENTS:\n")
    for mode in modes_list:
        dis_list = [key for key in FM[mode]]
        dis_list.sort()
        print(f"    mode {mode}: {dis_list}")
    print("\n")

    for mode in modes_list:
        dis_list = [key for key in FM[mode]]
        dis_list.sort()
        for dis in dis_list:
            ISOK = True
            try:
                fm = FM[mode][dis]
            except KeyError:
                print(f"mode {mode}, displacement {dis/100:.2f} FM : problem with scf")
                ISOK = False
            try:
                afm1 = AFM1[mode][dis]
            except KeyError:
                print(
                    f"mode {mode}, displacement {dis/100:.2f} AFM1 : problem with scf"
                )
                ISOK = False
            try:
                afm2 = AFM2[mode][dis]
            except KeyError:
                print(
                    f"mode {mode}, displacement {dis/100:.2f} AFM2 : problem with scf"
                )
                ISOK = False
            try:
                afm3 = AFM3[mode][dis]
            except KeyError:
                print(
                    f"mode {mode}, displacement {dis/100:.2f} AFM3 : problem with scf"
                )
                ISOK = False
            if ISOK:
                J1, J2, J3 = J_values(fm, afm1, afm2, afm3)
                print(f"Mode {mode}, displacement {dis/100:.2f}")
                print(f"    J1 = {J1:.6f} J2 = {J2:.6f} J3 = {J3:.6f}")


if __name__ == "__main__":
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "filename",
        type=str,
        help="""File with grepped lines.
Example grep:

grep "!" *mode/*F*/{neg,pos}/crsbr.scf*out > summary.txt
    """,
    )
    args = parser.parse_args()

    compute(args.filename)

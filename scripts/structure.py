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

import math

import bpy


def cylinder_between(x1, y1, z1, x2, y2, z2, r, name="bond"):
    r = 1.5 * r

    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    dist = math.sqrt(dx**2 + dy**2 + dz**2)

    bpy.ops.mesh.primitive_cylinder_add(
        radius=r, depth=dist, location=(dx / 2 + x1, dy / 2 + y1, dz / 2 + z1)
    )
    activeObject = bpy.context.active_object
    for f in activeObject.data.polygons:
        f.use_smooth = True

    phi = math.atan2(dy, dx)
    theta = math.acos(dz / dist)

    activeObject.rotation_euler[1] = theta
    activeObject.rotation_euler[2] = phi
    if name is not None:
        activeObject.name = name
    return activeObject


def sphere(
    x,
    y,
    z,
    r,
    name=None,
):
    r = 2 * r

    bpy.ops.mesh.primitive_uv_sphere_add(radius=r, location=(x, y, z))
    activeObject = bpy.context.active_object
    for f in activeObject.data.polygons:
        f.use_smooth = True
    if name is not None:
        activeObject.name = name

    return activeObject


def distance(x1, y1, z1, x2, y2, z2, target):
    # print(math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2))
    if math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) <= target:
        return True
    else:
        return False


def all_atoms_and_bonds(filename):
    a = 3.506000
    b = 4.767020
    c = 23.570749
    with open(filename, "r") as file:
        data = file.readlines()
    data = [i.split() for i in data]
    radius = {
        "Cr": 166 / 53,
        "Br": 94 / 53,
        "S": 88 / 53 * 1.83,
        "Ni": 135 / 53 * 1.63,
        "P": 100 / 53 * 1.16,
        "Se": 115 / 53,
    }
    counter = {"Cr": 0, "Br": 0, "S": 0, "Ni": 0, "P": 0, "Se": 0}

    atoms = []
    for atom in data:
        x, y, z = 1 * float(atom[5]), 1 * float(atom[6]), 1 * float(atom[7])
        name = atom[0]
        atoms.append((name, x, y, z))

    for i in range(0, len(atoms) - 1):
        name1 = atoms[i][0]
        x1 = atoms[i][1]
        y1 = atoms[i][2]
        z1 = atoms[i][3]
        sphere(x1, y1, z1, r=radius[name1] / 10, name=f"{name1}.{counter[name1]}")
        counter[atom[0]] += 1
        for j in range(i + 1, len(atoms)):
            name2 = atoms[j][0]
            x2 = atoms[j][1]
            y2 = atoms[j][2]
            z2 = atoms[j][3]
            # if distance(x1, y1, z1, x2, y2, z2, 2.5):
            #     print(distance(x1, y1, z1, x2, y2, z2, 2.5))
            # if {name1, name2} == {"Ni", "S"} and distance(x1, y1, z1, x2, y2, z2, 2.5):
            #     cylinder_between(x1, y1, z1, x2, y2, z2, 0.05)
            # if {name1, name2} == {"Ni", "Se"} and distance(x1, y1, z1, x2, y2, z2, 2.6):
            #     cylinder_between(x1, y1, z1, x2, y2, z2, 0.05)
            # if {name1, name2} == {"P", "S"} and distance(x1, y1, z1, x2, y2, z2, 2.1):
            #     cylinder_between(x1, y1, z1, x2, y2, z2, 0.05)
            # if {name1, name2} == {"P", "Se"} and distance(x1, y1, z1, x2, y2, z2, 2.3):
            #     cylinder_between(x1, y1, z1, x2, y2, z2, 0.05)

            if {name1, name2} == {"Ni"} and distance(x1, y1, z1, x2, y2, z2, 0.91):
                cylinder_between(x1, y1, z1, x2, y2, z2, 0.05 * 2.58)
            if {name1, name2} == {"Ni", "S"} and distance(
                x1, y1, z1, x2, y2, z2, 3.299
            ):
                cylinder_between(x1, y1, z1, x2, y2, z2, 0.05 * 2.58)
            if {name1, name2} == {"P", "S"} and distance(x1, y1, z1, x2, y2, z2, 2.851):
                cylinder_between(x1, y1, z1, x2, y2, z2, 0.05 * 2.58)


atoms = all_atoms_and_bonds("/Users/rybakov.ad/Projects/Dalton/cover_art/final.txt")
